"""Parallel chunked downloader for FinnGen migraine GWAS files (supports range requests, resume)."""
import os, sys, time, threading, urllib.request

BASE = "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release"
OUT = r"D:\GWAS\Migraine"
FILES = [
    "finngen_R12_G6_MIGRAINE.gz",
    "finngen_R12_G6_MIGRAINE.gz.tbi",
    "finngen_R12_G6_MIGRAINE_WITH_AURA.gz",
    "finngen_R12_G6_MIGRAINE_WITH_AURA.gz.tbi",
    "finngen_R12_G6_MIGRAINE_NO_AURA.gz",
    "finngen_R12_G6_MIGRAINE_NO_AURA.gz.tbi",
]
NCHUNK = 8
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
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        print(f"SKIP {fname} ({os.path.getsize(dest)} bytes)", flush=True)
        return
    size = total_size(url)
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
    for i, s, e in chunks:
        t = threading.Thread(target=run, args=(i, s, e), daemon=True)
        t.start()
        threads.append(t)
    while any(t.is_alive() for t in threads):
        time.sleep(20)
        done = sum(min(progress[i], e - s + 1) for i, s, e in chunks)
        print(f"  progress {fname}: {done*100//size}% ({done//1024//1024} MiB)", flush=True)
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
            os.remove(part)
    print(f"DONE {fname} -> {os.path.getsize(dest)} bytes", flush=True)

if __name__ == "__main__":
    for f in FILES:
        download(f)
    print("ALL FINISHED", flush=True)
