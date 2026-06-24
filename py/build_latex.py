#!/usr/bin/env python3
"""Build the LaTeX manuals in docs/ to PDF, then remove the auxiliary clutter, leaving only the
essentials (.tex, .pdf, the figures, a .bib if present). Runnable from the command line.

    python3 py/build_latex.py              # build every docs/*.tex, then clean
    python3 py/build_latex.py phys app     # build named docs (a stem, PTCryspMC_<stem>, or a file)
    python3 py/build_latex.py --figures    # regenerate the py/fig_*.py PNGs first, then build
    python3 py/build_latex.py --clean       # only remove the aux clutter (no build)

Needs `latexmk` (TeX Live) on PATH. Exits non-zero if any document fails to build.
"""
import argparse
import glob
import os
import shutil
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOCS = os.path.join(REPO, "docs")

# Regeneratable LaTeX by-products. Suffix match (str.endswith) handles the compound
# ".synctex.gz" / ".fdb_latexmk" / ".run.xml" as well.
AUX_SUFFIXES = (".aux", ".log", ".out", ".toc", ".lof", ".lot", ".fls", ".fdb_latexmk",
                ".synctex.gz", ".bbl", ".blg", ".nav", ".snm", ".vrb", ".idx", ".ilg",
                ".ind", ".bcf", ".run.xml", ".xdv")


def resolve_texs(names):
    """Map CLI names to docs/*.tex paths; with no names, every docs/*.tex."""
    if not names:
        return sorted(glob.glob(os.path.join(DOCS, "*.tex")))
    out = []
    for n in names:
        base = os.path.basename(n)
        for cand in (base, base + ".tex", f"PTCryspMC_{base}.tex"):
            p = os.path.join(DOCS, cand)
            if os.path.isfile(p):
                out.append(p)
                break
        else:
            sys.exit(f"no such .tex for '{n}' in docs/")
    return out


def clean(verbose=False):
    """Remove every aux by-product from docs/, keeping .tex/.pdf/.bib/figures."""
    removed = 0
    for f in sorted(os.listdir(DOCS)):
        if f.endswith(AUX_SUFFIXES):
            os.remove(os.path.join(DOCS, f))
            removed += 1
            if verbose:
                print(f"  rm {f}")
    return removed


def build(tex):
    """Compile one .tex with latexmk (run in docs/ so figures/ resolves). Returns (ok, log)."""
    name = os.path.basename(tex)
    r = subprocess.run(["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error", name],
                       cwd=DOCS, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    ok = r.returncode == 0 and os.path.isfile(os.path.join(DOCS, name[:-4] + ".pdf"))
    return ok, r.stdout


def main():
    ap = argparse.ArgumentParser(description="Build the LaTeX docs and clean the aux clutter.")
    ap.add_argument("docs", nargs="*", help="doc names/stems (default: all docs/*.tex)")
    ap.add_argument("--figures", action="store_true", help="regenerate the py/fig_*.py PNGs first")
    ap.add_argument("--clean", action="store_true", help="only clean aux files, do not build")
    a = ap.parse_args()

    if a.clean:
        print(f"cleaned {clean(verbose=True)} auxiliary file(s) in docs/")
        return

    if shutil.which("latexmk") is None:
        sys.exit("latexmk not found on PATH (install TeX Live).")

    if a.figures:
        for fig in sorted(glob.glob(os.path.join(REPO, "py", "fig_*.py"))):
            print(f"figure: {os.path.basename(fig)}")
            subprocess.run([sys.executable, fig], check=False)

    ok_all = True
    for tex in resolve_texs(a.docs):
        name = os.path.basename(tex)
        print(f"building {name} ...", end=" ", flush=True)
        ok, log = build(tex)
        if ok:
            kb = os.path.getsize(os.path.join(DOCS, name[:-4] + ".pdf")) // 1024
            print(f"OK ({kb} KB)")
        else:
            ok_all = False
            print("FAILED")
            for line in log.splitlines():
                if line.startswith("!") or "Error" in line or "Undefined" in line:
                    print("    " + line)

    n = clean()
    kept = ", ".join(sorted(os.listdir(DOCS)))
    print(f"cleaned {n} aux file(s); docs/ now holds: {kept}")
    sys.exit(0 if ok_all else 1)


if __name__ == "__main__":
    main()
