#!/usr/bin/env python3
"""AMVA-seeded, certified rational reconstruction of a closed CoMoM system.

Two positive-demand, distinct single-server queues; two closed classes; p>=2.
The exact reference is computed AFTER reconstruction and is never used by it.
Only the single-class boundary is exact input to the reconstruction.
No third-party Python packages are required.
"""

import argparse
import json
import math
import subprocess
import tempfile
import time
from dataclasses import dataclass
from decimal import Decimal, localcontext
from fractions import Fraction as F
from pathlib import Path


@dataclass(frozen=True)
class Model:
    name: str
    p: int
    q: int
    demands: tuple = ((1, 2), (3, 1))
    delays: tuple = (3, 2)


CASES = [Model("small", 3, 3), Model("medium", 6, 6),
         Model("larger", 12, 12),
         Model("accurate_amva", 6, 6, ((1, 2), (3, 1)), (100, 80)),
         Model("low_delay", 8, 8, ((10, 5), (5, 9)), (0, 0))]
CASES.append(Model("close_class1", 3, 6, ((1000001, 5), (1000000, 9)), (1, 1)))


def amva(model, p, q):
    """Damped Bard-Schweitzer fixed point; no exact metrics are used."""
    pop = (p, q)
    demand, delay = model.demands, model.delays
    Q = [[float(pop[r]) / 3 for r in range(2)] for _ in range(2)]
    for iteration in range(20000):
        R = [[0.0, 0.0] for _ in range(2)]
        for k in range(2):
            total = sum(Q[k])
            for r in range(2):
                if pop[r]:
                    R[k][r] = demand[k][r] * (1 + total - Q[k][r] / pop[r])
        X = [pop[r] / (delay[r] + sum(R[k][r] for k in range(2)))
             if pop[r] else 0.0 for r in range(2)]
        new = [[X[r] * R[k][r] for r in range(2)] for k in range(2)]
        change = max(abs(new[k][r] - Q[k][r]) / (1 + abs(new[k][r]))
                     for k in range(2) for r in range(2))
        if change < 1e-13:
            return X, new, iteration + 1
        Q = [[(Q[k][r] + new[k][r]) / 2 for r in range(2)] for k in range(2)]
    raise RuntimeError("AMVA did not converge")


def single_class(model, population, replica):
    """Exact boundary at class-2 population zero, positive convolution."""
    v = [F(model.delays[0] ** n, math.factorial(n))
         for n in range(population + 1)]
    stations = [0, 1] + ([replica - 1] if replica else [])
    for k in stations:
        demand = model.demands[k][0]
        # In-place convolution with 1/(1-demand*z).
        for n in range(1, population + 1):
            v[n] += demand * v[n - 1]
    return v[population]


def add(row, column, value):
    row[column] = row.get(column, F(0)) + value
    if not row[column]:
        del row[column]


class System:
    """Common normalization G(target)=1; s=1/G(target) is an unknown."""
    def __init__(self, model, start=0, closed=True):
        self.model, self.closed = model, closed
        first = 1 if closed else start
        self.keys = [(n, t, k) for n in range(first, model.q + 1)
                     for t in range(3) for k in range(3)]
        if closed:
            self.keys.append("s")
        self.index = {key: i for i, key in enumerate(self.keys)}
        self.size = len(self.keys)
        self.rows, self.rhs = [], []
        self.boundary = {(t, k): single_class(model, model.p - t, k)
                         for t in range(3) for k in range(3)} if closed else {}
        for n in range(1 if closed else start + 1, model.q + 1):
            # CE and class-1 PC at offsets t=0,1; class-2 PC at t=0,1,2.
            for t in range(2):
                for k in (1, 2):
                    self.equation([(n,t,k,1), (n,t,0,-1),
                                   (n,t+1,k,-model.demands[k-1][0]),
                                   (n-1,t,k,-model.demands[k-1][1])])
                self.equation([(n,t,0,model.p-t),
                               (n,t+1,0,-model.delays[0]),
                               (n,t+1,1,-model.demands[0][0]),
                               (n,t+1,2,-model.demands[1][0])])
            for t in range(3):
                self.equation([(n,t,0,n), (n-1,t,0,-model.delays[1]),
                               (n-1,t,1,-model.demands[0][1]),
                               (n-1,t,2,-model.demands[1][1])])
        self.equation([(model.q,0,0,1)], F(1))

    def component(self, n, t, k, value=1):
        if self.closed and n == 0:
            return {self.index["s"]: F(value) * self.boundary[t,k]}
        return {self.index[n,t,k]: F(value)}

    def equation(self, terms, rhs=F(0)):
        row = {}
        for n,t,k,value in terms:
            if value:
                for column, coefficient in self.component(n,t,k,value).items():
                    add(row, column, coefficient)
        self.rows.append(row)
        self.rhs.append(rhs)

    def metrics(self):
        m = self.model
        rows = {"X1": self.component(m.q,1,0),
                "X2": self.component(m.q-1,0,0)}
        for k in (1, 2):
            rows[f"Q{k}1"] = self.component(m.q,1,k,m.demands[k-1][0])
            rows[f"Q{k}2"] = self.component(m.q-1,0,k,m.demands[k-1][1])
        return rows

    def residual(self, y):
        return [b - dot(row,y) for row,b in zip(self.rows,self.rhs)]


def dot(row, y):
    return sum((coefficient * y[j] for j,coefficient in row.items()), F(0))


def echelon(rows):
    """Exact sparse row echelon form, used ONLY for uniqueness certificates."""
    pivots = {}
    for original in rows:
        row = dict(original)
        while row:
            j = min(row)
            if j not in pivots:
                value = row[j]
                pivots[j] = {k: v / value for k,v in row.items()}
                break
            value = row[j]
            for k,v in pivots[j].items():
                add(row,k,-value*v)
    return pivots


def rank_audit(system):
    pivots = echelon(system.rows)
    rank = len(pivots)
    extended = len(echelon(system.rows + list(system.metrics().values())))
    audit = {"unknowns": system.size, "equations": len(system.rows),
             "rank": rank, "nullity": system.size-rank,
             "metric_ambiguity_dimension": extended-rank}
    if extended > rank:
        for free in range(system.size):
            if free in pivots:
                continue
            witness = [F(0)] * system.size
            witness[free] = F(1)
            for pivot in sorted(pivots, reverse=True):
                witness[pivot] = -sum((v*witness[k] for k,v in pivots[pivot].items()
                                       if k != pivot), F(0))
            assert all(dot(row,witness) == 0 for row in system.rows)
            changes = {name:dot(row,witness) for name,row in system.metrics().items()}
            if any(changes.values()):
                audit["nullspace_witness"] = {str(j):str(v) for j,v in enumerate(witness) if v}
                audit["witness_metric_changes"] = {name:str(v) for name,v in changes.items()}
                break
        assert "nullspace_witness" in audit
    return audit


class CorrectionLU:
    """Floating correction solver, equilibrated using only H, never the answer.

    Default arithmetic is binary64. Decimal precision is an explicit fallback,
    not an exact linear solve. The same certification checks apply to both.
    """
    def __init__(self, system, digits=None):
        self.digits = digits
        with localcontext() as context:
            context.prec = digits or 28
            self.factor(system)

    def convert(self, value):
        if self.digits is None:
            return float(value)
        value = F(value)
        return Decimal(value.numerator) / Decimal(value.denominator)

    def total(self, values):
        return math.fsum(values) if self.digits is None else sum(values, Decimal(0))

    def factor(self, system):
        self.size = n = system.size
        assert len(system.rows) == n
        zero,one = self.convert(0),self.convert(1)
        colmax = [zero] * n
        for row in system.rows:
            for j,v in row.items():
                colmax[j] = max(colmax[j], abs(self.convert(v)))
        self.columns = [one / v if v else one for v in colmax]
        self.rowscale = []
        lu = []
        for row in system.rows:
            dense = [zero] * n
            for j,v in row.items():
                dense[j] = self.convert(v) * self.columns[j]
            scale = one / max(map(abs,dense))
            self.rowscale.append(scale)
            lu.append([v*scale for v in dense])
        self.permutation = list(range(n))
        for j in range(n):
            pivot = max(range(j,n), key=lambda i: abs(lu[i][j]))
            if not lu[pivot][j]:
                raise ArithmeticError("Floating-point LU has a zero pivot")
            lu[pivot],lu[j] = lu[j],lu[pivot]
            self.permutation[pivot],self.permutation[j] = self.permutation[j],self.permutation[pivot]
            for i in range(j+1,n):
                lu[i][j] /= lu[j][j]
                factor = lu[i][j]
                if factor:
                    for k in range(j+1,n):
                        lu[i][k] -= factor*lu[j][k]
        self.lu = lu

    def solve(self, residual):
        with localcontext() as context:
            context.prec = self.digits or 28
            return self.backsolve(residual)

    def backsolve(self, residual):
        rhs = [self.convert(r)*scale for r,scale in zip(residual,self.rowscale)]
        magnitude = max(map(abs,rhs))
        if magnitude == 0:
            raise ArithmeticError("Nonzero exact residual underflows in correction arithmetic")
        b = [rhs[i]/magnitude for i in self.permutation]
        for i in range(self.size):
            b[i] -= self.total(self.lu[i][j]*b[j] for j in range(i))
        for i in reversed(range(self.size)):
            b[i] = (b[i]-self.total(self.lu[i][j]*b[j]
                                   for j in range(i+1,self.size))) / self.lu[i][i]
        return [b[i]*magnitude*self.columns[i] for i in range(self.size)]


def seeds(system):
    m = system.model
    cache = {}
    def approx(p,q):
        if (p,q) not in cache:
            cache[p,q] = amva(m,p,q)
        return cache[p,q]
    X,Q,_ = approx(m.p,m.q)
    target_metrics = {"X1": X[0], "X2": X[1]}
    target_metrics.update({f"Q{k+1}{r+1}":Q[k][r] for k in range(2) for r in range(2)})
    target = [F(0)] * system.size
    target[system.index[m.q,0,0]] = F(1)
    for name,row in system.metrics().items():
        assert len(row) == 1
        j,weight = next(iter(row.items()))
        target[j] = F(target_metrics[name]) / weight

    base_ratio = {m.q: 1.0}
    for n in reversed(range(m.q)):
        base_ratio[n] = base_ratio[n+1] * approx(m.p,n+1)[0][1]
    full = []
    for key in system.keys:
        if key == "s":
            value = base_ratio[0] / float(system.boundary[0,0])
        else:
            n,t,k = key
            value = base_ratio[n]
            for population in range(m.p,m.p-t,-1):
                value *= approx(population,n)[0][0]
            if k:
                value *= 1 + sum(approx(m.p-t,n)[1][k-1])
        full.append(F(value))
    # Preserve the supplied target AMVA metrics exactly in all warm starts.
    for row in system.metrics().values():
        j = next(iter(row))
        full[j] = target[j]
    return {"zero": [F(0)]*system.size,
            "target_amva_only": target, "nearby_amva": full}, target_metrics, len(cache)


def reconstruct(system, initial, solver, max_iterations):
    """Exact rational accumulation and residuals; speculative CF then exact check."""
    y = list(initial)
    history = []
    started = time.perf_counter()
    for iteration in range(max_iterations+1):
        residual = system.residual(y)
        norm = max(map(abs,residual))
        history.append(float(norm))
        if norm == 0:
            return y, {"status":"certified", "correction_solves":iteration,
                       "reconstruction_denominator_cap_bits":None,
                       "residual_max_history":history,
                       "seconds":time.perf_counter()-started}
        # No reference-derived denominator bound or error tolerance is used.
        for bits in (8,16,32,64,128,256,512):
            candidate = [v.limit_denominator(2**bits) for v in y]
            if all(v == 0 for v in system.residual(candidate)):
                return candidate, {"status":"certified", "correction_solves":iteration,
                                   "reconstruction_denominator_cap_bits":bits,
                                   "residual_max_history":history,
                                   "seconds":time.perf_counter()-started}
        if iteration == max_iterations:
            break
        try:
            correction = solver.solve(residual)
            y = [v+F(delta) for v,delta in zip(y,correction)]
        except (ArithmeticError, OverflowError, ValueError) as error:
            return y, {"status":"failed", "reason":str(error),
                          "correction_solves":iteration,
                          "residual_max_history":history}
    return y, {"status":"not_certified", "correction_solves":max_iterations,
                  "residual_max_history":history,
                  "seconds":time.perf_counter()-started}


def convolution(model, replica=0):
    """Independent full-population positive convolution, validation only."""
    p,q = model.p,model.q
    values = [[F(model.delays[0]**i,math.factorial(i))
               * F(model.delays[1]**j,math.factorial(j))
               for j in range(q+1)] for i in range(p+1)]
    for k in [0,1] + ([replica-1] if replica else []):
        local = [[F(math.comb(i+j,i) * model.demands[k][0]**i
                    * model.demands[k][1]**j)
                  for j in range(q+1)] for i in range(p+1)]
        updated = [[F(0) for j in range(q+1)] for i in range(p+1)]
        for i in range(p+1):
            for j in range(q+1):
                updated[i][j] = sum((local[a][b]*values[i-a][j-b]
                                     for a in range(i+1) for b in range(j+1)),F(0))
        values = updated
    return values


def check_repository_ca(model, references):
    repo = Path(__file__).resolve().parents[2]
    binary = repo / "bin" / "ca"
    if not binary.exists():
        raise FileNotFoundError(f"Required independent reference binary missing: {binary}")
    checks = [(model.p,model.q,0), (model.p-1,model.q,0), (model.p,model.q-1,0)]
    checks += [(model.p-1,model.q,k) for k in (1,2)]
    checks += [(model.p,model.q-1,k) for k in (1,2)]
    with tempfile.TemporaryDirectory(prefix="amva-exact-poc-") as temporary:
        path = Path(temporary)/"model.qn"
        for p,q,replica in checks:
            stations = [0,1]+([replica-1] if replica else [])
            body = ["2",f"{p} {q}"," ".join(map(str,model.delays)),str(len(stations))]
            body += ["1 " + " ".join(map(str,model.demands[k])) for k in stations]
            path.write_text("\n".join(body)+"\n")
            process = subprocess.run([str(binary),str(path),"-e"],capture_output=True,
                                     text=True,timeout=30,check=True)
            num,den = process.stdout.split()
            assert F(int(num),int(den)) == references[replica][p][q]
    return len(checks)


def run_case(model, max_iterations, check_ca):
    system = System(model)
    # Local systems carry free preceding bases. Extra nearby levels alone do
    # not constitute an exact boundary condition.
    local = rank_audit(System(model,start=model.q-1,closed=False))
    window = rank_audit(System(model,start=max(0,model.q-3),closed=False))
    closed = rank_audit(system)
    assert local["metric_ambiguity_dimension"] > 0
    assert window["metric_ambiguity_dimension"] > 0
    assert closed["rank"] == system.size
    initial,amva_metrics,amva_calls = seeds(system)
    solver = CorrectionLU(system)
    outputs,report = {},{}
    for name,guess in initial.items():
        outputs[name],report[name] = reconstruct(system,guess,solver,max_iterations)

    # Retry stalled runs using a predetermined precision, still without any
    # multiclass reference. The failed binary64 attempts remain in the report.
    retries, retry_outputs = {}, {}
    failed = [name for name,r in report.items() if r["status"] != "certified"]
    if failed:
        precise_solver = CorrectionLU(system,digits=80)
        for name in failed:
            retry_outputs[name],retries[name] = reconstruct(system,outputs[name],precise_solver,max_iterations)
            retries[name]["correction_decimal_digits"] = 80

    # The exact multiclass oracle is deliberately first called here, after
    # every reconstruction attempt has completed.
    references = [convolution(model,k) for k in range(3)]
    normalizer = references[0][model.p][model.q]
    exact = [1/normalizer if key == "s"
             else references[key[2]][model.p-key[1]][key[0]]/normalizer
             for key in system.keys]
    assert all(r == 0 for r in system.residual(exact))
    metrics = {name:dot(row,exact) for name,row in system.metrics().items()}
    for attempts,summaries in ((outputs,report),(retry_outputs,retries)):
        for name,value in attempts.items():
            summary = summaries[name]
            summary["final_metric_max_relative_error"] = float(max(
                abs(dot(system.metrics()[metric],value)/reference-1)
                for metric,reference in metrics.items()))
            if summary["status"] == "certified":
                assert value == exact
                summary["full_basis_matches_convolution"] = True
                summary["max_denominator_bits"] = max(v.denominator.bit_length() for v in value)
                summary["metrics"] = {metric:str(dot(row,value))
                                      for metric,row in system.metrics().items()}
    initial_error = max(abs(F(amva_metrics[name])/value-1) for name,value in metrics.items())
    return {"name":model.name, "populations":[model.p,model.q],
            "demands":model.demands,"delays":model.delays,
            "rank_audit":{"last_step":local,"three_step_window":window,"closed":closed},
            "amva_target_max_relative_error":float(initial_error),
            "amva_target_metrics":amva_metrics,
            "nearby_amva_model_evaluations":amva_calls,
            "runs":report,
            "higher_precision_retry":retries,
            "reference_max_denominator_bits":max(v.denominator.bit_length() for v in exact),
            "repository_ca_constants_checked":check_repository_ca(model,references) if check_ca else 0}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=[m.name for m in CASES], action="append")
    parser.add_argument("--iterations",type=int,default=16)
    parser.add_argument("--check-ca",action="store_true")
    parser.add_argument("--output",type=Path)
    args = parser.parse_args()
    selected = [m for m in CASES if not args.case or m.name in args.case]
    results = []
    for model in selected:
        result = run_case(model,args.iterations,args.check_ca)
        results.append(result)
        status = ", ".join(f"{name}={r['status']}({r['correction_solves']})"
                           for name,r in result["runs"].items())
        print(f"{model.name}: AMVA error={result['amva_target_max_relative_error']:.3e}; {status}",flush=True)
        if result["higher_precision_retry"]:
            retry = ", ".join(f"{name}={r['status']}({r['correction_solves']})"
                              for name,r in result["higher_precision_retry"].items())
            print(f"  80-digit correction retry: {retry}",flush=True)
    all_certified = all(
        run["status"] == "certified"
        or case["higher_precision_retry"].get(name,{}).get("status") == "certified"
        for case in results for name,run in case["runs"].items())
    document = {"method":"closed normalized two-class CoMoM, exact residuals, floating corrections, rational reconstruction",
                "all_cases_eventually_certified":all_certified,
                "all_cases_binary64_certified":all(run["status"] == "certified"
                                                  for case in results for run in case["runs"].values()),
                "repository_ca_constants_checked":sum(case["repository_ca_constants_checked"] for case in results),
                "cases":results}
    if args.output:
        args.output.write_text(json.dumps(document,indent=2)+"\n")
    else:
        print(json.dumps(document,indent=2))
    if not all_certified:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
