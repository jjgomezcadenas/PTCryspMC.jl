#!/usr/bin/env zsh
# Collect every production's control-plot PNG into one flat gallery dir (prod/_gallery/) so the runs
# can be opened side by side. Copies prod/<tag>/<tag>_lors_det.png -> prod/_gallery/<tag>_lors_det.png.
# Idempotent: mirrors the current set (drops stale copies, refreshes the rest). run_prod.sh calls it
# after a batch, or run it by hand any time.
#
# Usage (from anywhere):
#   scripts/run/gallery.sh
# Env: PRODDIR (base, default prod).

set -u
cd ${0:A:h:h:h}                       # repo root (in scripts/run/)

PRODDIR=${PRODDIR:-prod}
GAL=$PRODDIR/_gallery
mkdir -p $GAL

pngs=($PRODDIR/*/*_lors_det.png(N))   # (N) = no error if none match
rm -f $GAL/*.png(N)                    # drop stale copies so the gallery mirrors the current runs
if (( ${#pngs} == 0 )); then
  print "no PNGs found under $PRODDIR/*/  (run scripts/run/run_prod.sh first)"
  exit 0
fi
for p in $pngs; do cp -f $p $GAL/; done
print "gallery: ${#pngs} PNG(s) -> $GAL/"
