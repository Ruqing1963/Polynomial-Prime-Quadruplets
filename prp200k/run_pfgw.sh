#!/bin/sh
# One PFGW instance per work_* directory (each needs its own pfgw.ini / pfgw.log).
# PRPs found are appended to work_*/pfgw.log.
PFGW=${PFGW:-pfgw64}
cd "$(dirname "$0")"
for d in work_*; do
  (cd "$d" && nohup nice "$PFGW" -f0 -l"pfgw_out.txt" candidates.txt > /dev/null 2>&1 &)
done
echo "started; check work_*/pfgw.log for PRP results"
