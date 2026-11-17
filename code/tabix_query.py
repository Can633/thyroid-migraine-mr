# -*- coding: utf-8 -*-
"""
Minimal tabix reader for remote bgzip files over HTTP Range requests.

Why: the GTEx v8 Thyroid eQTL all-pairs file is 4.1 GB and the link runs at
~0.5 MB/s, so a full download is not viable. The .tbi index is only 2.2 MB, and
tabix chunks let us fetch just the bytes covering our candidate loci.

Implements just enough of the tabix spec (SAMtools tbx format 1):
  parse_index()  -> per-reference bins and linear index
  reg2bins()     -> bins overlapping a query interval
  fetch()        -> HTTP Range -> bgzip members -> text lines in the interval

Only tested against tabix-built bgzf files with chromosome/position columns in
the first two fields (the eQTL Catalogue layout).
"""
import gzip
import io
import struct
import zlib
import urllib.request

MAX_BIN = 37450  # bin number used for the metadata pseudo-bin


def _read_index(path):
    with gzip.open(path, "rb") as fh:
        buf = fh.read()
    if buf[:4] != b"TBI\x01":
        raise ValueError("not a tabix index")
    off = 4
    n_ref, fmt, col_seq, col_beg, col_end, meta, skip, l_nm = struct.unpack_from("<8i", buf, off)
    off += 32
    names = buf[off:off + l_nm].split(b"\x00")
    names = [n.decode() for n in names if n]
    off += l_nm
    refs = []
    for _ in range(n_ref):
        n_bin = struct.unpack_from("<i", buf, off)[0]
        off += 4
        bins = {}
        for _ in range(n_bin):
            bin_id, n_chunk = struct.unpack_from("<Ii", buf, off)
            off += 8
            chunks = []
            for _ in range(n_chunk):
                b, e = struct.unpack_from("<QQ", buf, off)
                off += 16
                chunks.append((b, e))
            bins[bin_id] = chunks
        n_intv = struct.unpack_from("<i", buf, off)[0]
        off += 4
        linear = []
        for _ in range(n_intv):
            linear.append(struct.unpack_from("<Q", buf, off)[0])
            off += 8
        refs.append(dict(bins=bins, linear=linear))
    return dict(names=names, refs=refs, fmt=fmt, col_seq=col_seq,
                col_beg=col_beg, col_end=col_end, meta=meta, skip=skip)


def reg2bins(beg, end):
    """0-based half-open interval -> overlapping bin numbers."""
    end -= 1
    out = [0]
    for k in range(1 + (beg >> 26), 1 + (end >> 26) + 1):
        out.append(k)
    for k in range(9 + (beg >> 23), 9 + (end >> 23) + 1):
        out.append(k)
    for k in range(73 + (beg >> 20), 73 + (end >> 20) + 1):
        out.append(k)
    for k in range(585 + (beg >> 17), 585 + (end >> 17) + 1):
        out.append(k)
    for k in range(4681 + (beg >> 14), 4681 + (end >> 14) + 1):
        out.append(k)
    return out


def _bgzip_members(data, min_uoffset=0):
    """Decompress concatenated bgzip members; drop bytes before min_uoffset of
    the first member."""
    out = []
    pos = 0
    first = True
    while pos < len(data) - 1:
        d = zlib.decompressobj(31)
        try:
            chunk = d.decompress(data[pos:])
        except zlib.error:
            break
        if first and min_uoffset:
            chunk = chunk[min_uoffset:] if len(chunk) > min_uoffset else b""
            first = False
        else:
            first = False
        out.append(chunk)
        if not d.eof:
            break
        consumed = len(data) - pos - len(d.unused_data)
        if consumed <= 0:
            break
        pos += consumed
    return b"".join(out)


def fetch(url, index, chrom, start, end, timeout=120, max_bytes=200_000_000):
    """Return list of (pos, line) tuples for 1-based closed interval
    [start, end] on `chrom`."""
    if chrom not in index["names"]:
        alt = "chr" + chrom if not chrom.startswith("chr") else chrom[3:]
        if alt not in index["names"]:
            return []
        chrom = alt
    rid = index["names"].index(chrom)
    ref = index["refs"][rid]
    beg0, end0 = max(start - 1, 0), end
    bins = reg2bins(beg0, end0)
    chunks = []
    for b in bins:
        chunks.extend(ref["bins"].get(b, []))
    if not chunks:
        return []
    # linear index gives a lower bound on the starting offset
    li = ref["linear"][beg0 >> 14] if (beg0 >> 14) < len(ref["linear"]) else 0
    chunks = [(max(cb, li), ce) for cb, ce in chunks if ce > li]
    if not chunks:
        return []
    chunks.sort()
    merged = []
    for cb, ce in chunks:
        if merged and cb <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], ce))
        else:
            merged.append((cb, ce))

    out = []
    seq_i = index["col_seq"] - 1   # 0-based column of the sequence name
    beg_i = index["col_beg"] - 1   # 0-based column of the position
    for cb, ce in merged:
        c_off, u_off = cb >> 16, cb & 0xFFFF
        e_off = ce >> 16
        # fetch to the end of the last member; +64KB guard for member trailer
        req = urllib.request.Request(
            url, headers={"Range": f"bytes={c_off}-{e_off + 65536}",
                          "User-Agent": "tabix-remote/1.0"})
        with urllib.request.urlopen(req, timeout=timeout) as r:
            data = r.read()
        if len(data) > max_bytes:
            data = data[:max_bytes]
        text = _bgzip_members(data, min_uoffset=u_off).decode("utf-8", "replace")
        # drop a possibly truncated first line
        nl = text.find("\n")
        if nl > 0 and u_off and not text.startswith(("\n",)):
            head = text[:nl]
            # only drop if the fragment does not look like a complete record
            if len(head.split("\t")) < len(index.get("_ncol", [])) and head.count("\t") < beg_i:
                text = text[nl + 1:]
        for line in text.split("\n"):
            if not line:
                continue
            f = line.split("\t")
            if len(f) <= max(seq_i, beg_i):
                continue
            if f[seq_i] != chrom and f[seq_i] != chrom.replace("chr", ""):
                continue
            try:
                p = int(f[beg_i])
            except ValueError:
                continue
            if start <= p <= end:
                out.append((p, line))
    return out


def open_index(path):
    return _read_index(path)


if __name__ == "__main__":
    import sys
    idx = open_index("D:/GWAS/Migraine/eqtl/Thyroid.tsv.gz.tbi")
    print("refs:", idx["names"][:5], "... n=", len(idx["names"]))
    print("col_seq/beg/end:", idx["col_seq"], idx["col_beg"], idx["col_end"])
    url = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge/Thyroid.tsv.gz"
    hits = fetch(url, idx, "9", 133_255_000, 133_260_000)
    print("rows in chr9:133.255-133.260Mb ->", len(hits))
    for p, line in hits[:3]:
        print(" ", line[:160])
