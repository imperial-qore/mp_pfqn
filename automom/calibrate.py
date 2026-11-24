#!/usr/bin/env python3
"""Calibrate automom's cost model.

Runs every candidate solver on a grid of random non-degenerate models, fits
    log t = c0 + c1 log W + c2 log b
per solver (W = the solver's structural work and b = the bit-size estimate,
both from `bin/automom --features`), writes automom/calib.h and reports how
often the predicted-fastest solver is the measured-fastest one.

    python3 automom/calibrate.py --workdir /scratch/bitbucket/automom-cal -p 4
    make -C automom          # rebuild with the new constants

Each run is single-threaded; a chain (solver, R, population profile) stops
growing M at the first run that exceeds --tmax seconds.

--frac calibrates the models with fractional multiplicities instead: the
same grid with multiplicities drawn from FRAC_MI, only the solvers that
accept them (mom, comom, gcomom, gmom, gmomb), constants written to
automom/calib_frac.h (CALIB_FRAC, used by automom when qn_mden > 1).
"""
import argparse, json, math, os, random, subprocess, sys, time, zlib
from concurrent.futures import ThreadPoolExecutor

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, "bin")
# order = automom/main.c enum (calib.h is indexed by it); gmomb = gmom --bmax R-1
SOLVERS = ["ca", "mom", "comom", "gcomom", "gmom", "binmom", "bincomom", "mommod", "comommod", "gmomb"]
MS = [1, 2, 3, 4, 6, 8, 12, 16, 24]
PROFILES = {"u2": lambda R: [2] * R, "u8": lambda R: [8] * R,
            "u25": lambda R: [25] * R, "tail200": lambda R: [2] * (R - 1) + [200]}
FRAC_SOLVERS = ["mom", "comom", "gcomom", "gmom", "gmomb"]
FRAC_MI = ["1/2", "3/2", "2/3", "5/4", "7/3", "1", "2"]

def write_model(path, R, M, N, seed, frac=False):
    rnd = random.Random(seed)
    zpos = rnd.random() < 0.5                      # all-zero or all-positive think times
    if frac:                                       # at least one fractional station
        mi = [rnd.choice(FRAC_MI) for _ in range(M)]
        mi[rnd.randrange(M)] = rnd.choice(FRAC_MI[:5])
    else:
        mi = [1] * M if rnd.random() < 0.7 else [rnd.choice([1, 2]) for _ in range(M)]
    with open(path, "w") as f:
        f.write(f"{R}\n{' '.join(map(str, N))}\n")
        f.write(" ".join(str(rnd.randint(1, 10**6) if zpos else 0) for _ in range(R)) + "\n")
        f.write(f"{M}\n")
        for k in range(M):
            f.write(f"{mi[k]} " + " ".join(str(rnd.randint(1, 10**6)) for _ in range(R)) + "\n")

def features(path):
    out = subprocess.run([os.path.join(BIN, "automom"), "--features", path], capture_output=True, text=True).stdout.split("\n")
    head = out[0].split()
    logb = float(head[head.index("logb") + 1])
    W = {}
    for line in out[1:]:
        p = line.split()
        if len(p) == 3: W[p[0]] = (float(p[1]), p[2] == "1")
    return logb, W

def timed(solver, path, tmax):
    t = time.time()
    try:
        # through automom --solver, which adds gmom's -X and gmomb's --bmax
        p = subprocess.run([os.path.join(BIN, "automom"), "--solver", solver, "-t", path],
                           capture_output=True, text=True, timeout=tmax)
        return (time.time() - t, p.returncode)
    except subprocess.TimeoutExpired:
        return (tmax, "timeout")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", required=True)
    ap.add_argument("-p", "--parallel", type=int, default=4)
    ap.add_argument("--tmax", type=float, default=20.0)
    ap.add_argument("--maxR", type=int, default=5)
    ap.add_argument("--reuse", action="store_true", help="refit from workdir/runs.json without running")
    ap.add_argument("--frac", action="store_true", help="fractional multiplicities, writes calib_frac.h")
    a = ap.parse_args()
    solvers = FRAC_SOLVERS if a.frac else SOLVERS
    os.makedirs(a.workdir, exist_ok=True)
    runs_path = os.path.join(a.workdir, "runs.json")

    if not a.reuse:
        models = {}
        for R in range(1, a.maxR + 1):
            for pn, pf in PROFILES.items():
                for M in MS:
                    path = os.path.join(a.workdir, f"R{R}_{pn}_M{M}.qn")
                    write_model(path, R, M, pf(R), zlib.crc32(f"{R}/{pn}/{M}".encode()), a.frac)   # hash() is salted per process
                    models[(R, pn, M)] = (path,) + features(path)
        runs = []
        def chain(solver, R, pn):
            for M in MS:
                path, logb, W = models[(R, pn, M)]
                w, ok = W[solver]
                if not ok or not math.isfinite(w): continue
                t, rc = timed(solver, path, a.tmax)
                runs.append(dict(solver=solver, R=R, prof=pn, M=M, path=path, W=w, logb=logb, t=t, rc=rc))
                print(f"{solver:9s} R={R} {pn:7s} M={M:2d} t={t:7.3f} rc={rc}", flush=True)
                if rc == "timeout": break
        jobs = [(s, R, pn) for s in solvers for R in range(1, a.maxR + 1) for pn in PROFILES]
        random.Random(1).shuffle(jobs)                 # spread heavy chains over the pool
        with ThreadPoolExecutor(a.parallel) as ex:
            list(ex.map(lambda j: chain(*j), jobs))
        json.dump(runs, open(runs_path, "w"))
    runs = json.load(open(runs_path))

    coef = {}
    for s in solvers:
        pts = [r for r in runs if r["solver"] == s and r["rc"] == 0 and r["t"] >= 0.02]
        if len(pts) < 5:
            coef[s] = None; print(f"{s}: only {len(pts)} usable points, default kept"); continue
        X = np.array([[1.0, math.log(max(r["W"], 1.0)), r["logb"]] for r in pts])
        y = np.array([math.log(r["t"]) for r in pts])
        c, *_ = np.linalg.lstsq(X, y, rcond=None)
        err = np.abs(X @ c - y)
        coef[s] = c
        print(f"{s:9s} n={len(pts):3d}  c0={c[0]:8.3f} c1={c[1]:6.3f} c2={c[2]:6.3f}  "
              f"median |log err|={np.median(err):.2f} (x{math.exp(np.median(err)):.1f})")

    default = {"mommod": (-19.0, 1.0, 0.0), "comommod": (-19.0, 1.0, 0.0)}
    if a.frac:
        # same enum order as CALIB; solvers that never see a fractional model
        # keep their integer constants (not used: automom filters them out)
        with open(os.path.join(ROOT, "automom", "calib_frac.h"), "w") as f:
            f.write("/* Generated by automom/calibrate.py --frac: constants for models with a\n"
                    " * fractional multiplicity, same form and order as calib.h. */\n")
            f.write("#ifndef AUTOMOM_CALIB_FRAC_H\n#define AUTOMOM_CALIB_FRAC_H\n#include \"calib.h\"\nstatic const calib_t CALIB_FRAC[] = {\n")
            for i, s in enumerate(SOLVERS):
                c = coef.get(s)
                if c is None:
                    f.write(f"\t{{ 0, 0, 0 }},  /* {s}: from calib.h */\n")
                else:
                    f.write(f"\t{{ {c[0]:.6f}, {c[1]:.6f}, {c[2]:.6f} }},  /* {s} */\n")
            f.write("};\n#endif\n")
    else:
      with open(os.path.join(ROOT, "automom", "calib.h"), "w") as f:
        f.write("/* Generated by automom/calibrate.py.\n * Predicted seconds = exp(c0) * W^c1 * b^c2 (see main.c). */\n")
        f.write("#ifndef AUTOMOM_CALIB_H\n#define AUTOMOM_CALIB_H\ntypedef struct { double c0, c1, c2; } calib_t;\nstatic const calib_t CALIB[] = {\n")
        for s in SOLVERS:
            c = coef[s] if coef[s] is not None else default.get(s, (-16.0, 1.0, 1.0))
            f.write(f"\t{{ {c[0]:.6f}, {c[1]:.6f}, {c[2]:.6f} }},  /* {s} */\n")
        f.write("};\n#endif\n")

    # selection quality on the grid: chosen (predicted fastest) vs measured fastest
    by_model = {}
    for r in runs: by_model.setdefault(r["path"], []).append(r)
    hit = tot = 0; regret = []
    for path, rs in by_model.items():
        done = [r for r in rs if r["rc"] == 0]
        if len(done) < 2: continue
        def pred(r):
            c = coef[r["solver"]]
            return math.inf if c is None else c[0] + c[1] * math.log(max(r["W"], 1.0)) + c[2] * r["logb"]
        best = min(rs, key=lambda r: r["t"] if r["rc"] == 0 else math.inf)
        cand = [r for r in rs if r["rc"] in (0, "timeout")]
        pick = min(cand, key=pred)
        tot += 1; hit += pick["solver"] == best["solver"]
        regret.append(pick["t"] / max(best["t"], 1e-3))
    if tot:
        regret.sort()
        print(f"selection: fastest picked on {hit}/{tot} grid models; time/best median {regret[len(regret)//2]:.2f}, "
              f"90th pct {regret[int(0.9*(len(regret)-1))]:.2f}, worst {regret[-1]:.1f}")

if __name__ == "__main__":
    main()
