#!/usr/bin/env python3
"""Optional annotator hook - skip if index missing."""
from __future__ import annotations

import argparse
import pathlib
import subprocess
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sv", required=True)
    ap.add_argument("--pcawg-index", default="")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    
    if not args.pcawg_index or not pathlib.Path(args.pcawg_index).exists():
        print("[annotate] PCAWG index not found; skipping annotation.", flush=True)
        # pass-through: copy input to output if paths differ
        if args.sv != args.out:
            pathlib.Path(args.out).write_text(pathlib.Path(args.sv).read_text())
        return
    
    # call your annotator if present; otherwise skip
    if not pathlib.Path("annotate_sv_catalogs.py").exists():
        print("[annotate] annotate_sv_catalogs.py not present; skipping.", flush=True)
        if args.sv != args.out:
            pathlib.Path(args.out).write_text(pathlib.Path(args.sv).read_text())
        return
    
    cmd = [
        sys.executable, "annotate_sv_catalogs.py",
        "--sv", args.sv, "--pcawg", args.pcawg_index,
        "--out", args.out, "--junction-window", "500", "--interval-ovl", "0.90",
    ]
    print("[annotate] Running:", " ".join(cmd), flush=True)
    subprocess.check_call(cmd)


if __name__ == "__main__":
    main()

