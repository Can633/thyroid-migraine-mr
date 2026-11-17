"""Parallel chunked downloader for FinnGen R13 migraine GWAS (range requests + resume)."""
import os
import time
import threading
import urllib.request

BASE = "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90473001-GCST90474000/GCST90473326"
OUT = r"D:\GWAS\Migraine"
FILES = ["GCST90473326.tsv.gz"]
NCHUNK = 16
CHUNK = 16 * 1024 * 1024  # 16 MiB


def total_size(url):
    req = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(req, timeout=60) as r:
        return int(r.headers["Content-Length"])


def fetch_range(url, start, end, dest, idx, progress, lock):
    part = dest + f".part{idx}"
    have = os.path.getsize(part) if os.path.exists(part) else 0
    s = start + have
    if s > end:
        return
    for attempt in range(6):
        try:
            req = urllib.request.Request(url, headers={"Range": f"bytes={s}-{end}"})
            with urllib.request.urlopen(req, timeout=180) as r, open(part, "ab") as f:
                while True:
                    b = r.read(1 << 20)
                    if not b:
                        break
                    f.write(b)
                    with lock:
                        progress[idx] = f.tell()
            return
        except Exception as e:
            have = os.path.getsize(part) if os.path.exists(part) else 0
            s = start + have
            if s > end:
                return
            print(f"  chunk {idx} retry {attempt+1}: {e}", flush=True)
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"chunk {idx} failed")


def download(fname):
    url = f"{BASE}/{fname}"
    dest = os.path.join(OUT, fname)
    size = total_size(url)
    if os.path.exists(dest) and os.path.getsize(dest) == size:
        print(f"SKIP {fname} (complete, {size} bytes)", flush=True)
        return
    if os.path.exists(dest):
        print(f"REMOVE partial {fname} ({os.path.getsize(dest)} bytes)", flush=True)
        os.remove(dest)

    print(f"DOWNLOAD {fname} size={size} bytes", flush=True)
    chunks = [(i, i * CHUNK, min((i + 1) * CHUNK - 1, size - 1))
              for i in range((size + CHUNK - 1) // CHUNK)]
    progress = {i: 0 for i, _, _ in chunks}
    lock = threading.Lock()
    threads = []
    sem = threading.Semaphore(NCHUNK)

    def run(i, s, e):
        with sem:
            fetch_range(url, s, e, dest, i, progress, lock)

    t0 = time.time()
    for i, s, e in chunks:
        t = threading.Thread(target=run, args=(i, s, e), daemon=True)
        t.start()
        threads.append(t)
    while any(t.is_alive() for t in threads):
        time.sleep(30)
        done = sum(min(progress[i], e - s + 1) for i, s, e in chunks)
        el = time.time() - t0
        print(f"  progress {fname}: {done*100//size}% ({done//1024//1024} MiB, "
              f"{done/max(el,1)/1024:.0f} KiB/s)", flush=True)
    for t in threads:
        t.join()

    with open(dest, "wb") as out:
        for i, _, _ in chunks:
            part = dest + f".part{i}"
            with open(part, "rb") as f:
                while True:
                    b = f.read(1 << 20)
                    if not b:
                        break
                    out.write(b)
            # 分片保留，避免大批量删除触发安全策略；下次运行会因目标文件已完整而跳过
    print(f"DONE {fname} -> {os.path.getsize(dest)} bytes in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    for f in FILES:
        download(f)
    print("ALL FINISHED", flush=True)
