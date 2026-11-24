# oica: exact convolution for networks with order-independent stations

`oica` implements, literally, the exact algorithm of the OI/pass-and-swap notes (overleaf
`26_oipas`, Sec. 3.1), as an independent check of `oimom`:

- `G` by the shift-vector recursion (eq. oiconv): for an OI station `m`,
  `g(m,n,t) = g(m-1,n,0) + sum_{r: n_r>0} v_{m,r}/mu_m(t+e_r) g(m,n-e_r,t+e_r)`; for a
  load-independent one the plain convolution step; `g(0,n) = prod sigma_r^{n_r}/n_r!`. The shifted
  table of an OI station lives on the pairs `n + t <= N`.
- `X_r = G(N-e_r)/G(N)` (eq. tput).
- `Q_{i,r} = F_{i,r}(N)/G(N)`, `F_{i,r}(N) = sum_{0<n<=N} f_{i,r}(n) g(M,N-n,0)`, with `f` from
  `(f * Phi_i)(n) = (1+n_r) Phi_i(n)`, `f(0) = 1`, and `Phi_i` the station's balance function
  (with visit ratios). The `n = 0` term excluded from `F` is exactly the extra `G` that the `1+n_r`
  adds, so `F = Q G`.

Exact rationals throughout. The input is oimom's `.qnoi` format. The OI rate is a function of the
macrostate (property P1) and is evaluated only in `oi_rate()`; it currently reads the rates given
per set of classes present, and a general `mu(n)` needs only a new input section and that function.

    oica [-l|-e|-g|-t|-q|-d] [--check] model.qnoi

`--check` verifies Little's law on the delay, `sum_i Q_{i,r} + sigma_r X_r = N_r`, exactly.

## Validation

On the ten `oimom/models/*.qnoi` and nine copies with populations 7 to 15, `oica -d` prints the
same exact report as `oimom -d` (G, every X and every Q), and `--check` holds on all of them. This
also checks the notes' formulas for `G`, `X` and `Q`.

## Cost

About `prod_r (N_r+1)(N_r+2)/2` table entries per OI station for `G`, and `prod_r (N_r+1)^2` per
station and class for the `Q` deconvolution. R=3, N=(12,12,12): 1.3 s for `G`, 3.8 s with `Q`.
