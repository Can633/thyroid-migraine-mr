"""Finish FinnGen R12 downloads without deleting any file (safe-delete shim friendly).
Complete 16MiB .partN files are reused; partial ones keep their downloaded prefix
and only the missing tail is fetched in 1MiB micro-chunks via many connections."""
import os, time, threading, urllib.request

BASE = "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release"
OUT = r"D:\GWAS\Migraine"
BIG = 16 * 1024 * 1024
MICRO = 1 * 1024 * 1024
MAXWORK = 32
FILES = ["finngen_R12_G6_MIGRAINE.gz", "finngen_R12_G6_MIGRAINE.gz.tbi",
         "finngen_R12_G6_MIGRAINE_WITH_AURA.gz", "finngen_R12_G6_MIGRAINE_WITH_AURA.gz.tbi",
         "finngen_R12_G6_MIGRAINE_NO_AURA.gz", "finngen_R12_G6_MIGRAINE_NO_AURA.gz.tbi"]

def size_of(url):
    with urllib.request.urlopen(urllib.request.Request(url, method="HEAD"), timeout=60) as r:
        return int(r.headers["Content-Length"])

def fetch(url, s, e, dest, lock, prog, key):
    for attempt in range(8):
        try:
            req = urllib.request.Request(url, headers={"Range": f"bytes={s}-{e}"})
            with urllib.request.urlopen(req, timeout=120) as r, open(dest, "wb") as f:
                while True:
                    b = r.read(1 << 18)
                    if not b: break
                    f.write(b)
                    with lock: prog[key] = prog.get(key, 0) + len(b)
            return True
        except Exception as ex:
            if attempt == 7:
                print(f"FAIL {dest}: {ex}", flush=True); return False
            time.sleep(1.5 * (attempt + 1))
    return False

def process(fname):
    url = f"{BASE}/{fname}"
    dest = os.path.join(OUT, fname)
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        print(f"SKIP {fname} ({os.path.getsize(dest)} bytes)", flush=True); return
    size = size_of(url)
    print(f"FILE {fname} size={size}", flush=True)
    nb = (size + BIG - 1) // BIG
    todo = []
    for bi in range(nb):
        bs, be = bi*BIG, min((bi+1)*BIG, size)-1
        bigpart = dest + f".part{bi}"
        have = os.path.getsize(bigpart) if os.path.exists(bigpart) else 0
        if have >= be - bs + 1:
            continue  # complete big chunk
        # fetch only missing tail [bs+have, be]
        for x in range(bs+have, be+1, MICRO):
            mp = dest + f".m{bi}_{x}"
            if os.path.exists(mp) and os.path.getsize(mp) == min(x+MICRO-1, be)-x+1:
                continue
            todo.append((x, min(x+MICRO-1, be), bi))
    print(f"  {len(todo)} micro-chunks to fetch", flush=True)
    lock = threading.Lock(); prog = {}; sem = threading.Semaphore(MAXWORK); threads = []
    def run(s, e, bi):
        with sem:
            fetch(url, s, e, dest + f".m{bi}_{s}", lock, prog, fname)
    for s, e, bi in todo:
        t = threading.Thread(target=run, args=(s, e, bi), daemon=True)
        t.start(); threads.append(t)
    while any(t.is_alive() for t in threads):
        time.sleep(15)
        print(f"  {fname}: fetched {prog.get(fname,0)//1024//1024} MiB this session", flush=True)
    for t in threads: t.join()
    # assemble (leave all parts in place; just write final file)
    got = 0
    with open(dest, "wb") as out:
        for bi in range(nb):
            bs, be = bi*BIG, min((bi+1)*BIG, size)-1
            bigpart = dest + f".part{bi}"
            if os.path.exists(bigpart):
                with open(bigpart, "rb") as f:
                    while True:
                        b = f.read(1<<20)
                        if not b: break
                        out.write(b); got += len(b)
            for x in range(bs, be+1, MICRO):
                mp = dest + f".m{bi}_{x}"
                if os.path.exists(mp):
                    with open(mp, "rb") as f:
                        while True:
                            b = f.read(1<<20)
                            if not b: break
                            out.write(b); got += len(b)
    print(f"DONE {fname}: {got} bytes (expected {size})", flush=True)

for f in FILES:
    process(f)
print("ALL FINISHED", flush=True)
