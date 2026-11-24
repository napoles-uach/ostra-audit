#!/bin/bash
# job.sh species xc tag   (tag: libxc | bundled | patched | conv)
X=/home/claude/src_xcfun
case $3 in
  libxc) be=libxc; pre="";; bundled) be=xcfun; pre="";;
  patched) be=xcfun; pre=$X/patched/b/lib/libxcfun.so.2;; conv) be=xcfun; pre=$X/conv/b/lib/libxcfun.so.2;;
esac
save=""; [ "$3" = patched ] && save=1
out=$(SAVEDM=$save LD_PRELOAD=$pre OMP_NUM_THREADS=1 timeout 900 python3 run1.py "$1" "$2" "$be" def2-tzvp 2>/dev/null)
[ -z "$out" ] && out="{\"sp\": \"$1\", \"xc\": \"$2\", \"E\": null, \"err\": \"crash/timeout\"}"
echo "$out" | python3 -c "import json,sys;d=json.loads(sys.stdin.read());d['tag']='$3';print(json.dumps(d))"
