"""Extract the 88 IV SNPs from each FinnGen outcome into small TSV files."""
import os, gzip
root = r"D:\GWAS\Migraine"
want = set(open(os.path.join(root, "iv_rsids.txt")).read().split())
print("want", len(want), "SNPs", flush=True)
for code in ["G6_MIGRAINE", "G6_MIGRAINE_WITH_AURA", "G6_MIGRAINE_NO_AURA"]:
    src = os.path.join(root, f"finngen_R12_{code}.gz")
    dst = os.path.join(root, f"outcome_{code}.tsv")
    n = 0
    with gzip.open(src, "rt") if False else open(src, "rb") as raw:
        import gzip
        with gzip.GzipFile(fileobj=raw) as f:
            header = f.readline().decode().rstrip("\n").split("\t")
            idx = {k: header.index(k) for k in ["rsids","ref","alt","pval","beta","sebeta","af_alt"]}
            with open(dst, "w", encoding="utf-8") as out:
                out.write("rsid\tref\talt\tpval_out\tbeta_out\tse_out\teaf_out\n")
                for line in f:
                    parts = line.decode().rstrip("\n").split("\t")
                    if len(parts) < len(header): continue
                    rs = parts[idx["rsids"]]
                    if rs in want:
                        out.write("\t".join([rs, parts[idx["ref"]], parts[idx["alt"]],
                                             parts[idx["pval"]], parts[idx["beta"]],
                                             parts[idx["sebeta"]], parts[idx["af_alt"]]]) + "\n")
                        n += 1
    print(code, "->", n, "rows", flush=True)
print("DONE", flush=True)
