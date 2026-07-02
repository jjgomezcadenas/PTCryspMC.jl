#!/usr/bin/env zsh
# Prune a production run down to its deliverable: keep lors_det.h5 + config.toml, delete the
# regenerable intermediates (singles.h5, lors_truth.h5, randoms.h5). Useful for ≥1e8 runs where the
# intermediates are several GB each.
#
# Everything deleted here regenerates bit-identically from config.toml (pinned + propagated
# seed/nchunks/nevents) via run_prod.sh — in ~3 min for a 1e8 run. config.toml is the regeneration
# recipe and is NEVER deleted. lors_truth.h5 is the DT/τ-study input (examine_dt.jl); prune it only
# once τ is fixed (it is: 3 ns for CsI & BGO).
#
# SAFETY: a dir is pruned only if BOTH lors_det.h5 and config.toml are present — i.e. a complete run.
# An incomplete/failed run (no lors_det.h5) is skipped, so pruning never loses an unfinished result.
#
# Usage (from anywhere):
#   scripts/run/prune_prod.sh                                 # every prod/*/
#   scripts/run/prune_prod.sh sphere_air_bgo_10MBq            # named run(s)
#   DRYRUN=1 scripts/run/prune_prod.sh                        # show what would be deleted, delete nothing
#
# Env: PRODDIR (base, default prod), DRYRUN (set → preview only).

set -u
cd ${0:A:h:h:h}                       # repo root (in scripts/run/)

PRODDIR=${PRODDIR:-prod}
KEEP=(lors_det.h5 config.toml)
PRUNE=(singles.h5 lors_truth.h5 randoms.h5)

# Collect target dirs: the args (tag or prod/tag both accepted), else all prod/*/.
if (( $# )); then
  dirs=()
  for a in "$@"; do dirs+=($PRODDIR/${a:t}); done
else
  dirs=($PRODDIR/*(N/))               # (N/) = directories only, no error if none
fi
(( ${#dirs} )) || { print "no run dirs found under $PRODDIR/"; exit 1 }

[[ -n ${DRYRUN:-} ]] && print "DRY RUN — nothing will be deleted\n"
total=0
for d in $dirs; do
  tag=${d:t}
  [[ $tag == _* ]] && continue        # skip meta dirs (e.g. _gallery), not runs
  if [[ ! -d $d ]]; then print "skip $tag: no such dir ($d)"; continue; fi
  if [[ ! -f $d/lors_det.h5 ]]; then print "skip $tag: no lors_det.h5 (incomplete run — not pruning)"; continue; fi
  if [[ ! -f $d/config.toml ]]; then print "skip $tag: no config.toml (cannot regenerate — not pruning)"; continue; fi

  freed=0; victims=()
  for f in $PRUNE; do
    [[ -f $d/$f ]] || continue
    bytes=$(stat -f%z $d/$f)          # macOS stat
    freed=$((freed + bytes)); victims+=($f)
  done
  if (( ${#victims} == 0 )); then print "$tag: already pruned"; continue; fi

  mb=$((freed / 1048576))
  if [[ -n ${DRYRUN:-} ]]; then
    print "$tag: would delete ${(j:, :)victims}  (${mb} MB)"
  else
    for f in $victims; do rm -f $d/$f; done
    print "$tag: deleted ${(j:, :)victims}  (freed ${mb} MB)"
  fi
  total=$((total + freed))
done
print "\ntotal $([[ -n ${DRYRUN:-} ]] && print reclaimable || print freed): $((total / 1048576)) MB"
