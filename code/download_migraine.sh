#!/usr/bin/env bash
set -u
BASE="https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release"
OUT="D:/GWAS/Migraine"
mkdir -p "$OUT"
for code in G6_MIGRAINE G6_MIGRAINE_WITH_AURA G6_MIGRAINE_NO_AURA; do
  for ext in gz gz.tbi; do
    f="finngen_R12_${code}.${ext}"
    url="$BASE/$f"
    dest="$OUT/$f"
    if [ -s "$dest" ]; then echo "SKIP existing $f"; continue; fi
    echo "DOWNLOAD $f"
    curl -fL --retry 5 --retry-delay 3 -C - --max-time 1800 "$url" -o "$dest" || echo "FAIL $f"
    echo "DONE $f $(stat -c%s "$dest" 2>/dev/null || echo 0)"
  done
done
echo "ALL FINISHED"
