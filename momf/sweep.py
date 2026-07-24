#!/usr/bin/env python3
"""Accuracy of fixed-precision MoM against the exact rational MoM.

Generates a family of closed multiclass models of growing population,
runs bin/mom (exact rational) and bin/momf (fixed precision, optionally
with Wilkinson/Moler iterative refinement) on each, and reports the
relative error of the floating-point run.

Metrics
  dlogG : |log G_fp - log G_exact|, i.e. the relative error of G itself
  eX    : max relative error over class throughputs
  eQ    : max relative error over per-station-per-class queue lengths

Usage
  python3 sweep.py --base <model.qn> --scale 1 2 4 8 ... [-- extra momf args]
"""
import argparse
import math
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MOM = os.path.join(ROOT, "bin", "mom")
MOMF = os.path.join(ROOT, "bin", "momf")


def read_qn(path):
    toks = []
    with open(path) as fh:
        for line in fh:
            line = line.split("#")[0].strip()
            if line:
                toks.append(line.split())
    R = int(toks[0][0])
    N = [int(x) for x in toks[1][:R]]
    Z = [x for x in toks[2][:R]]
    M = int(toks[3][0])
    rows = [toks[4 + i][: R + 1] for i in range(M)]
    return R, N, Z, M, rows


def write_qn(path, R, N, Z, M, rows):
    with open(path, "w") as fh:
        fh.write("%d\n" % R)
        fh.write(" ".join(str(x) for x in N) + "\n")
        fh.write(" ".join(str(x) for x in Z) + "\n")
        fh.write("%d\n" % M)
        for row in rows:
            fh.write(" ".join(str(x) for x in row) + "\n")


def run(binary, model, flag, extra=()):
    cmd = [binary, model, flag] + list(extra)
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    if p.returncode != 0:
        return None
    vals = []
    for line in p.stdout.strip().splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            vals.extend(float(x) for x in line.split())
        except ValueError:
            return None
    return vals


def log_int(n):
    """ln(n) for an arbitrary-size positive integer"""
    b = n.bit_length()
    shift = max(0, b - 60)
    return math.log(n >> shift) + shift * math.log(2.0)


def exact_logG(model):
    """log G from the exact solver's num/den output.

    bin/mom -l converts G to a double before taking the log, so it
    overflows to inf for any model with more than a few hundred jobs.
    That is a defect of the output path, not of the computation: the
    rational G is exact.  Reading -e and taking the log on the integers
    recovers it."""
    p = subprocess.run([MOM, model, "-e"], capture_output=True, text=True,
                       timeout=1800)
    if p.returncode != 0:
        return None
    lines = [l.strip() for l in p.stdout.strip().splitlines() if l.strip()]
    if len(lines) < 2:
        return None
    try:
        num, den = int(lines[0]), int(lines[1])
    except ValueError:
        return None
    if num <= 0 or den <= 0:
        return None
    return log_int(num) - log_int(den)


def relerr(a, b):
    """max relative error of a against reference b"""
    worst = 0.0
    for x, y in zip(a, b):
        if y == 0.0:
            worst = max(worst, abs(x))
        else:
            worst = max(worst, abs(x - y) / abs(y))
    return worst


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True)
    ap.add_argument("--scale", nargs="+", type=int, required=True)
    ap.add_argument("--momf-args", default="")
    ap.add_argument("--label", default="fp")
    args = ap.parse_args()
    args.momf_args = args.momf_args.split()

    R, N0, Z, M, rows = read_qn(args.base)
    tmpdir = tempfile.mkdtemp(prefix="momsweep")

    print("model=%s  R=%d M=%d  momf args: %s"
          % (os.path.basename(args.base), R, M, args.momf_args or "(none)"))
    print("%8s %10s %14s %14s %14s" % ("Ntot", "logG", "dlogG", "eX", "eQ"))

    for s in args.scale:
        N = [n * s for n in N0]
        mf = os.path.join(tmpdir, "m_%d.qn" % s)
        write_qn(mf, R, N, Z, M, rows)

        ex_lg = exact_logG(mf)
        ex_t = run(MOM, mf, "-t")
        ex_q = run(MOM, mf, "-q")
        if ex_lg is None or ex_t is None or ex_q is None:
            print("%8d  exact solver failed (singular?)" % sum(N))
            continue

        fp_l = run(MOMF, mf, "-l", args.momf_args)
        fp_t = run(MOMF, mf, "-t", args.momf_args)
        fp_q = run(MOMF, mf, "-q", args.momf_args)
        if fp_l is None or fp_t is None or fp_q is None:
            print("%8d  fp solver failed" % sum(N))
            continue

        dlog = abs(fp_l[0] - ex_lg)
        eX = relerr(fp_t, ex_t)
        eQ = relerr(fp_q, ex_q)
        print("%8d %10.4f %14.3e %14.3e %14.3e"
              % (sum(N), ex_lg, dlog, eX, eQ))
        sys.stdout.flush()


if __name__ == "__main__":
    main()
