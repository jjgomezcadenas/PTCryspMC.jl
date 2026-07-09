#!/usr/bin/env zsh
# Generate shards of an API config SEQUENTIALLY (each run already saturates the cores, so running
# shards in parallel would only thrash). For each realization R_LO..R_HI: run the full chain
# (overwriting prod/<tag>/), publish the shard into PtCryspProds, run light per-shard gates, then
# PRUNE the heavy intermediates. Stops on the first failure. A master check at the end.
# See dev/data_generation_strategy.md for the shard/master model.
#
#   PRODS=~/Projects/PtCryspProds THREADS=18 scripts/run/run_shards.sh uniform_headep_bgo_api 1 9

set -u
cd ${0:A:h:h:h}

CFG=runs/${1:t:r}.toml
LO=${2:?usage: run_shards.sh <config> <R_lo> <R_hi>}
HI=${3:?usage: run_shards.sh <config> <R_lo> <R_hi>}
PRODS=${PRODS:-$HOME/Projects/PtCryspProds}
THREADS=${THREADS:-$(sysctl -n hw.ncpu)}
NCHUNKS=${NCHUNKS:-144}
[[ -f $CFG ]] || { print "missing config: $CFG"; exit 1 }
tag=${CFG:t:r}; dir=prod/$tag

print "shards $LO..$HI of $tag  →  $PRODS   (threads=$THREADS, nchunks=$NCHUNKS)"
prev_nev=""
for N in $(seq $LO $HI); do
  log=/tmp/shard_${tag}_${N}.log
  s3=$(printf '%03d' $N)
  print "\n=== shard $s3 ==="
  if ! {
    julia -t $THREADS --project=. scripts/simulate_source_mt.jl --config $CFG --realization $N --nchunks $NCHUNKS --format hdf5 &&
    julia -t $THREADS --project=. scripts/build_true_coincidences_from_singles.jl --config $CFG &&
    julia --project=. scripts/build_randoms_from_singles.jl --config $CFG &&
    julia --project=. scripts/reco_lors.jl --config $CFG &&
    julia --project=. scripts/run/publish_prod.jl --config $CFG --prods-root $PRODS --force
  } >$log 2>&1; then
    print "  CHAIN FAILED — tail of $log:"; tail -6 $log; exit 1
  fi

  # --- light per-shard gates (from the logs) ---
  nev=$(grep -oE '[0-9]+ events, seed' $log | grep -oE '^[0-9]+' | head -1)
  ratio=$(grep -oE 'measured/analytic = [0-9.]+' $log | grep -oE '[0-9.]+$' | head -1)
  ok=1
  grep -q "shard $s3" $log || { print "  GATE FAIL: publish did not report shard $s3 (realization mismatch)"; ok=0 }
  awk -v r="${ratio:-0}" 'BEGIN{exit !(r>0.7 && r<1.3)}' || { print "  GATE FAIL: randoms ratio ${ratio:-?} off 1"; ok=0 }
  [[ -n $prev_nev && $nev == $prev_nev ]] && { print "  GATE FAIL: nevents==prev ($nev) — realization not varying"; ok=0 }
  (( ok )) || { print "  gates failed for shard $s3 — see $log"; exit 1 }
  print "  OK: nevents=$nev, randoms measured/analytic=$ratio, published lors_shard$s3.h5"
  prev_nev=$nev

  # --- prune ALL heavy .h5: the shard is exported, so prod/ keeps no .h5 (only config.toml).
  # Guard: delete the source lors_det.h5 only once THIS run's published shard is confirmed on
  # disk. The leaf comes from the publish log ("  leaf: <scenario>/<scanner>/<crystal>/<bd>/"),
  # NOT from a tree glob — a glob matches every master in the tree and picks the alphabetically
  # first (the frozen reference), silently validating the wrong leaf once several masters exist.
  leafrel=$(grep -oE 'leaf: [^ ]+' $log | head -1); leafrel=${leafrel#leaf: }; leafrel=${leafrel%/}
  published=$PRODS/$leafrel/lors_shard$s3.h5
  if [[ -n $leafrel && -f $published ]]; then
    rm -f $dir/singles.h5 $dir/lors_truth.h5 $dir/randoms.h5 $dir/lors_det.h5
  else
    print "  WARN: published shard $s3 not found — keeping lors_det.h5"
    rm -f $dir/singles.h5 $dir/lors_truth.h5 $dir/randoms.h5
  fi
done

# --- master check: all shards present, self-describing, with distinct realizations ---
# $leafrel survives the loop (same leaf every shard) — the check runs on THIS run's leaf.
print "\n=== master check ==="
leaf=$PRODS/$leafrel
[[ -n $leafrel && -d $leaf ]] || { print "  no published leaf (no 'leaf:' line in the publish log)"; exit 1 }
julia --project=. -e '
using HDF5
d = ARGS[1]
fs = sort(filter(f -> startswith(f, "lors_shard") && endswith(f, ".h5"), readdir(d)))
reals = Int[]; nevs = Int[]
for f in fs
    h5open(joinpath(d, f)) do h
        push!(reals, Int(read(attributes(h)["realization"]))); push!(nevs, Int(read(attributes(h)["nevents"])))
    end
end
println("  leaf: ", d)
println("  shards: ", length(fs), "  realizations: ", reals)
println("  nevents: ", extrema(nevs), "  pooled ΣM = ", sum(nevs))
length(unique(reals)) == length(reals) || error("duplicate realizations — shards NOT independent")
println("  OK: ", length(unique(reals)), " distinct realizations, all self-describing.")
' "$leaf" 2>&1 | grep -vE "Precompiling|precompiled"
