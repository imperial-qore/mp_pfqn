# oimom: MoM-family analysis with order-independent stations

`oimom` computes the exact normalizing constant, throughputs and mean queue lengths of a closed
multiclass product-form network made of a delay, load-independent single-server queues (with
multiplicity `mi`) and **order-independent (OI) stations whose rate depends only on the set of
classes present**: compatibility pools (LINE's `sn_compat_rate`), balanced-fairness bandwidth
sharing, redundancy with cancel-on-completion, any table `mu_S`.

## Model

An OI station has visit weights `v_r` and a total rate `mu_S` for each nonempty set `S` of classes
present. Its balance function is

    Phi(n) = (1/mu_supp(n)) sum_{r: n_r>0} v_r Phi(n - e_r),   Phi(0) = 1,

the convention of LINE's `pfqn_ncoi` and `Pfqn_clwoi`. With `mu_S = 1` for every `S` the station is a
load-independent queue with demands `v`.

## How it works

Unrolling that recursion over the supports gives the station's generating function as a sum over
chains `{} < S_1 < ... < S_k` (each step adds one class `r_i`):

    F(z) = sum_chains prod_i (v_{r_i}/mu_{S_i}) z_{r_i} / (1 - theta_{S_i} . z),
    theta_{S,r} = v_r/mu_S for r in S, 0 otherwise.

Each factor `1/(1 - theta_S . z)` is a load-independent queue. So `G(N)` is a weighted sum, over
one chain per OI station, of normalizing constants of ordinary networks in which each OI station is
replaced by one **virtual queue** per support of its chain, at population `N` minus one job of each
class of the chain's last support. Each of those networks is solved exactly by a zero-aware CoMoM
recursion (`zc.c`, below), or, when that fails, by `bin/comom`'s tiered route, linked in-process
(`inner.c` includes `comom/main.c` with `main` renamed, as `bincomom` does). The constants they
export give the rest:

    G(N-1_r)  = sum_c w_c G_c(P_c - 1_r)
    Q_kr G(N) = sum_c w_c m_k L_kr G_c^{+k}(P_c - 1_r)                         LI queue k
    Q_jr G(N) = sum_c w_c ( [r in S_j] G_c(P_c)
                            + sum_{T in chain_j} theta_{T,r} G_c^{+T}(P_c - 1_r) )   OI station j

The first term of the last line counts the one class-r job that the chain's shift places at the
station; the second is the derivative of the virtual queues, a replica in MoM terms.

Details:

- Derived demands are rational; each class is scaled by the lcm of its denominators, and every
  exported constant at population `n` is divided by `prod_r D_r^{n_r}` afterwards.
- Classes with no jobs left and stations with zero demand on every remaining class are dropped
  from a derived network; a network with no stations left is a pure delay, solved in closed form.
- A chain choice that leaves a class with jobs and nowhere to go contributes `G_c = 0`, but can
  still contribute to `G(N-1_r)` when that class is `r` and holds one job; that case is solved
  separately at `P_c - 1_r`.
- `--conv` solves the network instead by convolving the station tables over the population
  lattice (the oracle); `--validate` runs both and requires `G`, every `X_r` and every `Q_kr` to
  agree exactly.

## Zero-aware inner recursion (`zc.c`)

The virtual queues have zero demand outside their support, and comom's square class system is
singular as soon as a station has a zero demand on a class processed earlier. `zc.c` runs the same
class-oriented recursion (basis `G(N_E - c, t)`, `G^{+k}(N_E - c, t)` for `|c| <= M` over the
earlier classes `E`; class-`r` PC for `G`, CEs and earlier-class PCs for `G^{+k}`) without
requiring the system to be square or nonsingular:

- Rows are chosen by exact elimination (the system is consistent, possibly overdetermined), and the
  null space of the remaining matrix is computed exactly: an unknown is used only if every null
  vector vanishes on it. The square block of pivot rows and pivot columns is put in block upper
  triangular form once per class: a perfect matching of rows to columns, then Tarjan's strongly
  connected components of the dependency graph, sinks first. Only the diagonal blocks are inverted. A
  step is a block back-substitution: subtract the sparse terms on the columns of blocks already
  solved, then apply the block inverse as integer dot products over a common denominator. The forward
  elimination stays sparse and costs little; the dense inverse of the whole block, which this replaces,
  filled in and was most of the run.
- Two defects cause the null space, both at the deepest level `|c| = M`. One is "orphans":
  `G^{+k}(N_E - d)` with `d` supported on classes where `k` has zero demand, which no equation
  touches. The other is stations proportional on two or more earlier classes, which partial zeros
  make common (`(2,0)` and `(2/7,0)`). Values below the deepest level stay determined. A deep `G`
  that the class-`r` PC cannot form (it needs an undetermined `G^{+k}`) becomes an unknown of the
  next system, determined through the `Z_s G` term of the earlier-class PCs when it is.
- Classes that visit no station are factored out; class orders are tried, the last successful one
  first.
- Only a value needed for the output (`G(N)`, `G(N-1_r)`, `G^{+k}(N-1_r)` with `L_kr > 0`) being
  undetermined fails a class order. When every order fails, the network goes to the tiered route.

A deep `G` that stays undetermined is still a plain normalizing constant: the network restricted to
the classes processed so far, at population `(N_E - d, t)`. `OIMOM_ZC_SUBRUNS=b` lets `zc.c`
compute up to `b` of them per class order by a recursive run at that population (memoized), which
makes the recursion complete in principle. It is off by default: on the models below each sub-run
costs a full solve, and they cut fallbacks less than they cost (`b3_16`: 4.2 s at 0, 9.8 s at 4,
the same two fallbacks).

`OIMOM_ZCHECK=1` checks every usable equation, not only the pivot rows, at every step.
`OIMOM_DEBUG=1` reports the route of each inner solve; `OIMOM_DEBUG=2` also the rank and
undetermined count of each factorization and the networks that fall back.

## Single basis (`--single`)

The chain decomposition runs one CoMoM recursion per combination of chains, `16^K` at R = 3. Over a
common denominator the whole model is one load-independent network instead. Each OI station is

    F_j(z) = P_j(z) / D_j(z),   D_j = prod_S (1 - theta_S . z),

`S` over the nonempty supports of visited classes, `P_j` a polynomial of degree `|D_j|` (every chain
term is a monomial times the factors of the supports it skips). `P_j` is computed as `D_j F_j` on the
box `n_r <= min(N_r, sum_j |D_j|)`, `F_j` from the balance recursion, and every factor of `D_j` that
divides `P_j` is cancelled (the quotient has degree below `|D_j|`). With `G'` the network of the LI queues
plus one queue per remaining support of every OI station (identical rows merged into `mi`), and
`p = prod_j P_j`:

    G(N)       = sum_a p_a G'(N-a)
    G(N-1_r)   = sum_a p_a G'(N-a-1_r)
    Q_kr G(N)  = m_k L_kr sum_a p_a G'^{+k}(N-a-1_r)                                   LI queue k
    Q_jr G(N)  = sum_a q^{jr}_a G'(N-a) + sum_{S in D_j} theta_{S,r} sum_a p_a G'^{+S}(N-a-1_r)

with `q^{jr} = (z_r d/dz_r P_j) prod_{i != j} P_i`: the derivative of the numerator, then that of the
denominator. One `zc.c` run on `G'` gives every shifted value (`zc_solve_shift`): the basis level is
raised to the largest shift over the earlier classes when it exceeds `M`, and the last class keeps its
steps `N_r - lev .. N_r`. An order is dropped as soon as one of its systems leaves every unknown
undetermined. If a needed value stays undetermined in every class order, oimom says so and uses the
chain decomposition.

Exact on all nineteen models: `--validate` on the ten small ones, and output identical to the chain
decomposition on the nine large ones and the five benchmarks. One large model (`r2k2_15`, R=2, two OI
stations) falls back. Time (one core, `-d`):

| model | `--single` | chains |
|---|---|---|
| R=3, 1 OI + 2 LI, N=(12,12,12) | 4.5 s | 1.1 s |
| R=3, 1 OI + 2 LI, N=(24,24,24) | 4.7 s | 3.7 s |
| R=3, 1 OI + 2 LI, N=(48,48,48) | 5.1 s | 25.8 s |
| R=3, 2 OI + 1 LI, N=(8,8,8) | 18.6 s | 37.8 s |
| R=3, 2 OI + 1 LI, N=(7,7,7) (`r3k2_7`) | 8.0 s | 9.4 s |
| R=3, 2 OI, zero think times, N=(8,8,8) (`z0_8`) | 3.3 s | 20.0 s |
| R=3, `oi_partial`, N=9 | 9.4 s | 7.3 s |

The single network is larger (8 queues at R = 3 with one OI station, 15 with two), and its class
systems have about 360 unknowns. Their block triangular form still has one diagonal block of about
200, and its exact inverse is most of the run. That cost does not depend on `N`, and the steps are
cheap. So `--single` wins at large populations and with several OI stations, and loses at small `N`
with one. Three of the four class orders fail on the R=3 benchmarks and cost about half the run.

## Usage

    oimom [-l|-e|-g|-t|-q|-d] [--conv] [--single] [--validate] [--comom-inner] [-j n] model.qnoi

Same output flags as the other solvers. `-q` prints one row per station: the LI queues first, then
the OI stations in file order. `--single` solves the model as one network (above). `--comom-inner` skips `zc.c` and sends every derived network to the
tiered CoMoM route (the previous behaviour). `-j` sets the threads of the interpolation tier of that
route (see Cost).

File format: a `.qn` file, then `OI K` and one row per OI station with `R` visit weights and the
`2^R - 1` rates `mu_S`, `S = 1 .. 2^R-1` as a bitmask (bit `r-1` set when class `r` is present).
Every value may be a rational `a/b`. Example (`models/oi_li.qnoi`, an OI station with every rate 1,
equal to a second LI queue with demands 4 3):

    2
    5 4
    3 2
    1
    1 2 5
    OI 1
    4 3 1 1 1

## Validation

`models/*.qnoi`: ten models, R = 2 and 3, one or two OI stations, with and without LI queues and
think times, `mi > 1`, a class with zero population, a class not visiting an OI station, and a
network served only by OI stations with zero think times (the degenerate-chain path). On every one
`--validate` passes: `G`, `X` and `Q` equal the lattice convolution exactly. `oi_li.qnoi` also
matches `bin/ca` byte for byte for `-e`, `-t`, `-q` on the equivalent LI model.

These models have `N <= M` in places, where the deepest corners of the basis are impossible
(`c_s > N_s`) and the defects above do not show. Nine copies with populations 7 to 15 exercise
them. `--validate` with `OIMOM_ZCHECK=1` passes on all nineteen. `zc.c` solves every derived network
of the ten small models, and of the large ones all but 1 of 16 (R=3, K=1), 2 of 25 (R=2, K=2),
20 of 80 (`oi_partial`) and 62 of 192 (`oi_z0_onlyoi`: zero think times, two OI stations).

## Cost

There are `sum_k R!/(R-k)!` chains per OI station (5 at R = 2, 16 at R = 3, 65 at R = 4), and one
inner solve per combination of chains, so the count is that number raised to the number of OI
stations. Each inner solve is polynomial in `N`, against roughly `prod_r (N_r+1)^2` for convolving
an OI balance function.

Through bin/comom, the system of every derived network with virtual queues is singular and goes to
exact interpolation (`Ntot+1` modular solves), which dominated the run before `zc.c`. Measured
(one core, `-e`, output identical in every column):

| model | oimom | `--comom-inner` | convolution (`--conv`) |
|---|---|---|---|
| R=2, 1 OI + 2 LI, N=(40,40) | 0.04 s | 0.08 s | 5.3 s |
| R=3, 1 OI + 2 LI, N=(12,12,12) | 1.1 s | 4.5 s | 4.4 s |
| R=3, 1 OI + 2 LI, N=(16,16,16) | 1.6 s | 9.4 s | 21.7 s |
| R=3, 1 OI + 2 LI, N=(24,24,24) | 3.7 s | 29.2 s | not run |
| R=3, 2 OI + 1 LI, N=(8,8,8) | 38 s | 122 s | 0.3 s |

In the R=3, one-OI rows 15 of the 16 inner solves go through `zc.c` and one falls back. Most of the
remaining time is that fallback and the class orders tried on it. With two OI stations every inner
solve goes through `zc.c`; the derived networks have up to 7 stations (about 250 unknowns). The
block triangular factorization (above) took the `zc.c` share of these rows down by 1.2x to 3x. oimom wins at large populations with
one OI station. With two or more OI stations at small `N`, the 16^K inner solves make convolution
the better choice.

## Scope and next steps

- Support-only rates only. OI rates that saturate at a count `tau > 1` decompose the same way, with
  one pole per saturation region and shifts up to `tau`, but are not implemented. Rates that never
  saturate have no finite decomposition.
- Closed networks only.
- Some derived networks still fall back (zero think times with several stations confined to two
  classes). The recursive sub-runs for the deep corners make them solvable but cost more than
  the fallback. A cheaper source for those normalizing constants would close the gap.
- `zc.c` inverts each diagonal block of the block triangular form densely. On the larger networks one
  block of about 200 unknowns remains; a sparse LU of that block would cut both its factorization and
  the per-step product.
- The chain combinations are independent and could run in parallel, but the inner CoMoM keeps
  global state, so they run one after another.
