#!/usr/bin/env zsh
# Make all coincidence plots IN PARALLEL: the per-config 9-panel summaries
# (plot_coincidences.py, one per config) plus the cross-config comparison
# (plot_matrix.py, once). Reads the existing output/<tag>/coincidences_*.csv, so run
# scripts/run/run_matrix.sh first to generate the data. Plotting is decoupled from the
# transport — re-run this freely after tweaking a plotter, no re-simulation needed.
#
# Usage (from anywhere):
#   scripts/run/plot_all.sh                                  # every runs/*.toml
#   scripts/run/plot_all.sh sphere_water_csi cylinder_water_bgo
#
# Logs: /tmp/ptcplot_<tag>.log (+ /tmp/ptcplot_matrix.log). Outputs:
#   output/<tag>/coincidences_{truth,det}.png   and   output/matrix_summary.png

set -u
cd ${0:A:h:h:h}                       # repo root

if (( $# )); then
  configs=()
  for a in "$@"; do configs+=(runs/${a:t:r}.toml); done
else
  configs=(runs/*.toml(N))
fi
(( ${#configs} )) || { print "no configs found in runs/"; exit 1 }

print "plotting ${#configs} config(s) + matrix in parallel: ${(j: :)${(@)configs:t:r}}"
start=$SECONDS

for cfg in $configs; do
  python3 py/plot_coincidences.py --config $cfg >/tmp/ptcplot_${cfg:t:r}.log 2>&1 &
done
python3 py/plot_matrix.py $configs >/tmp/ptcplot_matrix.log 2>&1 &
wait

print "\nplots done in $((SECONDS - start)) s wall\n"
for cfg in $configs; do
  tag=${cfg:t:r}
  printf "%-24s " $tag
  grep -h "wrote" /tmp/ptcplot_$tag.log 2>/dev/null || print "ERROR -> /tmp/ptcplot_$tag.log"
done
printf "%-24s " matrix
grep -h "wrote" /tmp/ptcplot_matrix.log 2>/dev/null || print "ERROR -> /tmp/ptcplot_matrix.log"
