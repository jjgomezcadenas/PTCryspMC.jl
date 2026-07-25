# AGENTS.md — PTCryspMC.jl: PET detector simulation

Working rules for any coding-agent session on this repo. `CLAUDE.md` holds the
repo orientation (purpose, modes, guides); read both before working.

## Working style

**Follow the user's explicit instructions and scope.** Never take an action the
user has not asked for. A question about whether an action is needed is not
authorization to perform it: answer the question and wait for an explicit
instruction before searching, reading, editing, running commands, or otherwise
acting.

**Never create, switch, or otherwise change Git branches with an uncommitted
worktree.** Before every branch operation, inspect the current branch and full
worktree status. Existing changes must first be committed on their current
branch, or the user must explicitly direct how they are to be handled.

**Produce every calculation and plot with tracked repository code.** Every
numerical result, data transformation, table and figure used in an analysis or
document must be reproducible from version-controlled code in this repository.
Terminal one-liners and REPL sessions are for diagnostics only; promote any
result used in the work to a tracked script before reporting it.

Ask questions plainly, in prose — no multiple-choice option menus. No defensive
writing in docs or comments: state what something is, not what it is not.
Terminal responses use bold for emphasis and plain text for paths and
identifiers (no backticks or links — they render blue in JJ's terminal).

## Active feature branch: `xsection-weighted-lors`

This branch carries the detector side of the cross-section uncertainty work.
The goal: each recorded coincidence remembers which source event produced it,
so that after the detector simulation has run once, the effect of exchanging
one production cross-section curve for another is computed by re-weighting the
recorded events instead of re-simulating them. Concretely, the branch will

1. add a second parent-event column to the coincidence schema, so accidental
   (random) coincidences keep both of the source events that formed them;
2. add a source mode that reads the shared source bank written by `ptcryspg4`,
   with a one-to-one mapping between the simulation event number and the bank's
   stable source-event identifier;
3. add a weights evaluator that joins the recorded coincidences with per-event
   weight tables and produces weighted activity profiles and distal-edge
   positions per cross-section replica;
4. keep the existing unweighted API path unchanged.

The governing plan is `docs/shared_plan.tex` in the `ptcryspg4` repo (its
feature branch there is `xsection-source-bank`). Detector-side implementation
details are recorded in `dev/xsection_weighted_lors_plan.md` as they are
decided.
