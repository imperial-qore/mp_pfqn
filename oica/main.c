/* oica - exact convolution algorithm for closed networks with order-independent
 * (OI) stations, as written in the OI/pass-and-swap notes (26_oipas, Sec. 3.1):
 *
 *  normalizing constant, eq. (oiconv): stations m = 1..K are OI, K+1..M
 *  load-independent (LI), station 0 the delay; with t_m the shift vector of
 *  OI station m,
 *    g(m,n,t) = g(m-1,n,0) + sum_{r: n_r>0} v_{m,r}/mu_m(t+e_r) g(m,n-e_r,t+e_r)   m <= K
 *    g(m,n,0) = g(m-1,n,0) + sum_{r: n_r>0} theta_{m,r} g(m,n-e_r,0)               m > K
 *    g(m,0,t) = 1,  g(0,n,.) = prod_r sigma_r^{n_r}/n_r!
 *  and G = g(M,N,0);
 *  throughputs, eq. (tput): X_r = G(N-e_r)/G(N);
 *  mean queue lengths: Q_{i,r} = F_{i,r}(N)/G(N), F_{i,r}(N) = sum_{0<n<=N}
 *  f_{i,r}(n) g(M,N-n,0), with f_{i,r} obtained by deconvolving
 *  (f_{i,r} * Phi_i)(n) = (1+n_r) Phi_i(n): f(0) = 1,
 *    f(n) = n_r Phi_i(n) - sum_{0 != k != n, k <= n} f(k) Phi_i(n-k),
 *  Phi_i the balance function of station i,
 *    Phi_i(0) = 1, Phi_i(n) = (1/mu_i(n)) sum_{r: n_r>0} v_{i,r} Phi_i(n-e_r)
 *  (mu_i = 1 and v = theta for an LI station).
 *
 * Exact rationals throughout.  The OI rate is a function of the macrostate
 * (property P1); oi_rate() is the only place that knows its form.  The input
 * is oimom's .qnoi format, rates mu_S given per set S of classes present.
 * A station with multiplicity mi is mi identical copies (convolved mi times,
 * its queue length reported as the total over the copies).  --check also
 * verifies Little's law on the delay, sum_i Q_{i,r} + sigma_r X_r = N_r. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <gmp.h>

typedef struct {
	int R, M, K;
	int* N;
	mpq_t* Z;
	mpq_t** L;   /* M x R */
	int* mi;
	mpq_t** v;   /* K x R */
	mpq_t** mu;  /* K x 2^R */
} model;

static double cputime(void)
{
	struct rusage u; getrusage(RUSAGE_SELF, &u);
	return u.ru_utime.tv_sec + u.ru_stime.tv_sec + 1e-6 * (u.ru_utime.tv_usec + u.ru_stime.tv_usec);
}

static int read_q(FILE* f, mpq_t q)
{
	if (gmp_fscanf(f, "%Qd", q) != 1) return 0;
	mpq_canonicalize(q);
	return 1;
}

static model* read_model(const char* path)
{
	FILE* f = fopen(path, "r");
	if (!f) { fprintf(stderr, "oica: cannot open %s\n", path); return NULL; }
	model* m = (model*) calloc(1, sizeof(model));
	int r, k, j, S;
	char kw[64];
	if (fscanf(f, "%d", &m->R) != 1 || m->R < 1 || m->R > 16) goto bad;
	m->N = (int*) calloc(m->R, sizeof(int));
	for (r = 0; r < m->R; r++) if (fscanf(f, "%d", &m->N[r]) != 1 || m->N[r] < 0) goto bad;
	m->Z = (mpq_t*) malloc(m->R * sizeof(mpq_t));
	for (r = 0; r < m->R; r++) { mpq_init(m->Z[r]); if (!read_q(f, m->Z[r])) goto bad; }
	if (fscanf(f, "%d", &m->M) != 1 || m->M < 0) goto bad;
	m->L = (mpq_t**) malloc((m->M ? m->M : 1) * sizeof(mpq_t*));
	m->mi = (int*) calloc(m->M ? m->M : 1, sizeof(int));
	for (k = 0; k < m->M; k++) {
		if (fscanf(f, "%d", &m->mi[k]) != 1 || m->mi[k] < 1) goto bad;
		m->L[k] = (mpq_t*) malloc(m->R * sizeof(mpq_t));
		for (r = 0; r < m->R; r++) { mpq_init(m->L[k][r]); if (!read_q(f, m->L[k][r])) goto bad; }
	}
	if (fscanf(f, "%63s", kw) == 1) {
		if (strcmp(kw, "OI") || fscanf(f, "%d", &m->K) != 1 || m->K < 0) goto bad;
		m->v = (mpq_t**) malloc((m->K ? m->K : 1) * sizeof(mpq_t*));
		m->mu = (mpq_t**) malloc((m->K ? m->K : 1) * sizeof(mpq_t*));
		for (j = 0; j < m->K; j++) {
			m->v[j] = (mpq_t*) malloc(m->R * sizeof(mpq_t));
			for (r = 0; r < m->R; r++) { mpq_init(m->v[j][r]); if (!read_q(f, m->v[j][r])) goto bad; }
			m->mu[j] = (mpq_t*) malloc((1 << m->R) * sizeof(mpq_t));
			mpq_init(m->mu[j][0]);
			for (S = 1; S < (1 << m->R); S++) { mpq_init(m->mu[j][S]); if (!read_q(f, m->mu[j][S])) goto bad; }
		}
	}
	fclose(f);
	for (j = 0; j < m->K; j++) {
		int vis = 0;
		for (r = 0; r < m->R; r++) if (mpq_sgn(m->v[j][r]) > 0) vis |= 1 << r;
		for (S = 1; S < (1 << m->R); S++)
			if ((S & ~vis) == 0 && mpq_sgn(m->mu[j][S]) <= 0) {
				fprintf(stderr, "oica: OI station %d needs a positive rate on support %d\n", m->M + j + 1, S);
				return NULL;
			}
	}
	return m;
bad:
	fprintf(stderr, "oica: malformed model file %s\n", path);
	fclose(f);
	return NULL;
}

/* total rate of OI station j in macrostate n (property P1: a function of the
 * macrostate only).  Here: the rate of the set of classes present. */
static void oi_rate(const model* m, int j, const int* n, mpq_t out)
{
	int r, S = 0;
	for (r = 0; r < m->R; r++) if (n[r] > 0) S |= 1 << r;
	mpq_set(out, m->mu[j][S]);
}

/* ---- population lattice 0 <= n <= N ---- */
typedef struct { int R; const int* N; long size; long* stride; } lat;

static lat lat_make(int R, const int* N)
{
	lat L; int r;
	L.R = R; L.N = N; L.stride = (long*) malloc(R * sizeof(long)); L.size = 1;
	for (r = 0; r < R; r++) { L.stride[r] = L.size; L.size *= N[r] + 1; }
	return L;
}
static void lat_vec(const lat* L, long i, int* n) { int r; for (r = 0; r < L->R; r++) { n[r] = i % (L->N[r] + 1); i /= L->N[r] + 1; } }
static long lat_idx(const lat* L, const int* n) { long i = 0; int r; for (r = 0; r < L->R; r++) i += n[r] * L->stride[r]; return i; }

static mpq_t* qarr(long n) { long i; mpq_t* a = (mpq_t*) malloc((n ? n : 1) * sizeof(mpq_t)); for (i = 0; i < n; i++) mpq_init(a[i]); return a; }
static void qfree(mpq_t* a, long n) { long i; for (i = 0; i < n; i++) mpq_clear(a[i]); free(a); }

/* OI station j on top of g = g(m-1,.,0): returns g(m,.,0).  The shifted
 * table h(n,t) lives on pairs n + t <= N, indexed per class by the
 * triangular index of (n_r, t_r). */
static void oi_step(const model* m, int j, const lat* L, mpq_t* g)
{
	int R = m->R, r, s;
	long* base = (long*) malloc(R * sizeof(long));
	long size = 1;
	for (r = 0; r < R; r++) { base[r] = size; size *= (long) (m->N[r] + 1) * (m->N[r] + 2) / 2; }
	/* tri(n,t) = u(u+1)/2 + t with u = n + t */
	#define TRI(nn, tt) ((long) ((nn) + (tt)) * ((nn) + (tt) + 1) / 2 + (tt))
	mpq_t* h = qarr(size);
	int n[32], t[32], u[32], tt[32];
	mpq_t w, rate; mpq_init(w); mpq_init(rate);
	/* h(n,t) needs h(n-e_r, t+e_r), which has the same u = n + t: for each
	 * u, take n <= u in increasing mixed-radix index */
	long nu = L->size, iu;
	for (iu = 0; iu < nu; iu++) {
			lat_vec(L, iu, u);
			long in, nn = 1;
			for (r = 0; r < R; r++) nn *= u[r] + 1;
			for (in = 0; in < nn; in++) {
				long x = in;
				for (r = 0; r < R; r++) { n[r] = x % (u[r] + 1); x /= u[r] + 1; t[r] = u[r] - n[r]; }
				long id = 0;
				for (r = 0; r < R; r++) id += TRI(n[r], t[r]) * base[r];
				mpq_set(h[id], g[lat_idx(L, n)]);       /* g(m-1, n, 0) */
				for (r = 0; r < R; r++) {
					if (n[r] == 0 || !mpq_sgn(m->v[j][r])) continue;
					memcpy(tt, t, R * sizeof(int)); tt[r]++;
					oi_rate(m, j, tt, rate);
					long id2 = 0;
					for (s = 0; s < R; s++) id2 += TRI(n[s] - (s == r), tt[s]) * base[s];
					mpq_div(w, m->v[j][r], rate);
					mpq_mul(w, w, h[id2]);
					mpq_add(h[id], h[id], w);
				}
			}
		}
	/* g(m, n, 0) */
	long i;
	for (i = 0; i < L->size; i++) {
		lat_vec(L, i, n);
		long id = 0;
		for (r = 0; r < R; r++) id += TRI(n[r], 0) * base[r];
		mpq_set(g[i], h[id]);
	}
	#undef TRI
	mpq_clear(w); mpq_clear(rate);
	qfree(h, size); free(base);
}

/* LI station with demands th on top of g: in place, by increasing index */
static void li_step(const model* m, mpq_t* th, const lat* L, mpq_t* g)
{
	int n[32], r; long i; mpq_t w; mpq_init(w);
	for (i = 0; i < L->size; i++) {
		lat_vec(L, i, n);
		for (r = 0; r < m->R; r++) {
			if (n[r] == 0 || !mpq_sgn(th[r])) continue;
			mpq_mul(w, th[r], g[i - L->stride[r]]);
			mpq_add(g[i], g[i], w);
		}
	}
	mpq_clear(w);
}

/* balance function of station s (s < M: LI queue, one copy; else OI station s-M) */
static mpq_t* balance(const model* m, int s, const lat* L)
{
	mpq_t* P = qarr(L->size);
	int n[32], r; long i; mpq_t w, rate; mpq_init(w); mpq_init(rate);
	mpq_set_ui(P[0], 1, 1);
	for (i = 1; i < L->size; i++) {
		lat_vec(L, i, n);
		for (r = 0; r < m->R; r++) {
			if (n[r] == 0) continue;
			mpq_mul(w, s < m->M ? m->L[s][r] : m->v[s - m->M][r], P[i - L->stride[r]]);
			mpq_add(P[i], P[i], w);
		}
		if (s >= m->M) { oi_rate(m, s - m->M, n, rate); mpq_div(P[i], P[i], rate); }
	}
	mpq_clear(w); mpq_clear(rate);
	return P;
}

/* convolution power: the balance function of mi identical copies */
static mpq_t* conv(const lat* L, mpq_t* a, mpq_t* b)
{
	mpq_t* c = qarr(L->size);
	int n[32], k[32], d[32], r; long i, j; mpq_t w; mpq_init(w);
	for (i = 0; i < L->size; i++) {
		lat_vec(L, i, n);
		for (j = 0; j < L->size; j++) {
			lat_vec(L, j, k);
			for (r = 0; r < L->R; r++) { if (k[r] > n[r]) break; d[r] = n[r] - k[r]; }
			if (r < L->R) continue;
			mpq_mul(w, a[j], b[lat_idx(L, d)]);
			mpq_add(c[i], c[i], w);
		}
	}
	mpq_clear(w);
	return c;
}

/* F_{i,r}(N) = sum_{0<n<=N} f(n) g(N-n), f from (f*Phi)(n) = (1+n_r) Phi(n) */
static void qlen(const lat* L, mpq_t* Phi, mpq_t* g, int r, mpq_t F)
{
	mpq_t* f = qarr(L->size);
	int n[32], k[32], d[32], s; long i, j; mpq_t w; mpq_init(w);
	mpq_set_ui(f[0], 1, 1);
	for (i = 1; i < L->size; i++) {
		lat_vec(L, i, n);
		mpq_set_si(w, n[r], 1); mpq_mul(f[i], w, Phi[i]);
		for (j = 1; j < i; j++) {                      /* 0 != k != n, k <= n */
			lat_vec(L, j, k);
			for (s = 0; s < L->R; s++) { if (k[s] > n[s]) break; d[s] = n[s] - k[s]; }
			if (s < L->R) continue;
			mpq_mul(w, f[j], Phi[lat_idx(L, d)]);
			mpq_sub(f[i], f[i], w);
		}
	}
	int N[32]; lat_vec(L, L->size - 1, N);
	mpq_set_ui(F, 0, 1);
	for (i = 1; i < L->size; i++) {
		lat_vec(L, i, n);
		for (s = 0; s < L->R; s++) d[s] = N[s] - n[s];
		mpq_mul(w, f[i], g[lat_idx(L, d)]);
		mpq_add(F, F, w);
	}
	mpq_clear(w); qfree(f, L->size);
}

static double qlog(const mpq_t q)
{
	long en, ed;
	double a = mpz_get_d_2exp(&en, mpq_numref(q)), b = mpz_get_d_2exp(&ed, mpq_denref(q));
	return log(a / b) + (en - ed) * log(2.0);
}
static void put_d(const mpq_t q) { printf("%.15e", mpq_get_d(q)); }

static void usage(const char* p)
{
	printf("USAGE: %s [-l|-e|-g|-t|-q|-d] [--check] model.qnoi\n", p);
	printf("  -l, --log    : log of the normalizing constant\n");
	printf("  -e, --ex     : exact normalizing constant, numerator and denominator\n");
	printf("  -g, --nc     : normalizing constant as double\n");
	printf("  -t, --tput   : throughputs, one per row\n");
	printf("  -q, --qlen   : mean queue lengths, one row per station (LI queues, then OI stations)\n");
	printf("  -d, --exact  : full report in exact rationals\n");
	printf("  --check      : verify Little's law on the delay exactly (exit 2 if it fails)\n");
	printf("Model: oimom's .qnoi format (a .qn file, then 'OI K' and K rows v_1..v_R mu_1..mu_{2^R-1}).\n");
	printf("Exact OI convolution of the 26_oipas notes, Sec. 3.1; cost about prod_r (N_r+1)(N_r+2)/2\n");
	printf("per OI station for G, prod_r (N_r+1)^2 per station and class for Q.\n");
}

int main(int argc, char** argv)
{
	int lg = 0, ex = 0, nc = 0, tp = 0, ql = 0, dbg = 0, chk = 0, i, r, k, j;
	const char* file = NULL;
	for (i = 1; i < argc; i++) {
		if (!strcmp(argv[i], "-l") || !strcmp(argv[i], "--log")) lg = 1;
		else if (!strcmp(argv[i], "-e") || !strcmp(argv[i], "--ex")) ex = 1;
		else if (!strcmp(argv[i], "-g") || !strcmp(argv[i], "--nc")) nc = 1;
		else if (!strcmp(argv[i], "-t") || !strcmp(argv[i], "--tput")) tp = 1;
		else if (!strcmp(argv[i], "-q") || !strcmp(argv[i], "--qlen")) ql = 1;
		else if (!strcmp(argv[i], "-d") || !strcmp(argv[i], "--exact")) dbg = 1;
		else if (!strcmp(argv[i], "--check")) chk = 1;
		else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) { usage(argv[0]); return 0; }
		else if (argv[i][0] == '-') { fprintf(stderr, "oica: unknown option %s\n", argv[i]); return 1; }
		else file = argv[i];
	}
	if (!file) { usage(argv[0]); return 1; }
	model* m = read_model(file);
	if (!m) return 1;
	double t0 = cputime();
	int R = m->R, S = m->M + m->K;
	lat L = lat_make(R, m->N);
	int n[32];
	long ii;
	/* station 0: the delay */
	mpq_t* g = qarr(L.size);
	mpq_t w; mpq_init(w);
	mpz_t fz; mpz_init(fz);
	for (ii = 0; ii < L.size; ii++) {
		lat_vec(&L, ii, n);
		mpq_set_ui(g[ii], 1, 1);
		for (r = 0; r < R; r++) {
			mpz_pow_ui(mpq_numref(w), mpq_numref(m->Z[r]), n[r]);
			mpz_pow_ui(mpq_denref(w), mpq_denref(m->Z[r]), n[r]);
			mpz_fac_ui(fz, n[r]); mpz_mul(mpq_denref(w), mpq_denref(w), fz);
			mpq_canonicalize(w);
			mpq_mul(g[ii], g[ii], w);
		}
	}
	/* OI stations m = 1..K, then LI stations (mi copies each) */
	for (j = 0; j < m->K; j++) oi_step(m, j, &L, g);
	for (k = 0; k < m->M; k++) for (i = 0; i < m->mi[k]; i++) li_step(m, m->L[k], &L, g);
	mpq_t G; mpq_init(G); mpq_set(G, g[L.size - 1]);
	if (!mpq_sgn(G)) { fprintf(stderr, "oica: G = 0 (a class with jobs has nowhere to go)\n"); return 1; }
	mpq_t* X = qarr(R);
	for (r = 0; r < R; r++)
		if (m->N[r] > 0) { lat_vec(&L, L.size - 1, n); n[r]--; mpq_div(X[r], g[lat_idx(&L, n)], G); }
	mpq_t** Q = (mpq_t**) malloc(S * sizeof(mpq_t*));
	int needQ = ql || chk || !(lg || ex || nc || tp);
	for (k = 0; k < S; k++) {
		Q[k] = qarr(R);
		if (!needQ) continue;
		mpq_t* Phi = balance(m, k, &L);
		if (k < m->M && m->mi[k] > 1) {   /* mi copies: convolution power */
			mpq_t* one = Phi; mpq_t* acc = Phi;
			for (i = 1; i < m->mi[k]; i++) { mpq_t* nx = conv(&L, acc, one); if (acc != one) qfree(acc, L.size); acc = nx; }
			qfree(one, L.size); Phi = acc;
		}
		for (r = 0; r < R; r++) { qlen(&L, Phi, g, r, Q[k][r]); mpq_div(Q[k][r], Q[k][r], G); }
		qfree(Phi, L.size);
	}
	int bad = 0;
	if (chk) {
		mpq_t tot; mpq_init(tot);
		for (r = 0; r < R; r++) {
			mpq_mul(tot, m->Z[r], X[r]);
			for (k = 0; k < S; k++) mpq_add(tot, tot, Q[k][r]);
			mpq_set_si(w, m->N[r], 1);
			if (!mpq_equal(tot, w)) { bad = 1; gmp_fprintf(stderr, "oica: Little's law fails for class %d: %Qd != %d\n", r + 1, tot, m->N[r]); }
		}
		mpq_clear(tot);
		fprintf(stderr, "oica: Little's law on the delay %s\n", bad ? "FAILS" : "holds exactly for every class");
	}
	if (lg) printf("%.15e\n", qlog(G));
	else if (ex) gmp_printf("%Zd\n%Zd\n", mpq_numref(G), mpq_denref(G));
	else if (nc) { put_d(G); printf("\n"); }
	else if (tp) for (r = 0; r < R; r++) { put_d(X[r]); printf("\n"); }
	else if (ql) for (k = 0; k < S; k++) { for (r = 0; r < R; r++) { put_d(Q[k][r]); if (r < R - 1) printf(" "); } printf("\n"); }
	else {
		printf("Model: R=%d, %d LI queue(s), %d OI station(s)\n", R, m->M, m->K);
		printf("\n========== Performance Metrics ==========\n");
		if (dbg) gmp_printf("G = %Qd\n", G); else { printf("G = "); put_d(G); printf("\nlog(G) = %.15e\n", qlog(G)); }
		printf("\nX (throughputs):\n");
		for (r = 0; r < R; r++) { if (dbg) gmp_printf("X[%d] = %Qd\n", r + 1, X[r]); else { printf("X[%d] = ", r + 1); put_d(X[r]); printf("\n"); } }
		printf("\nQ (mean queue lengths):\n");
		mpq_t tot; mpq_init(tot);
		for (k = 0; k < S; k++) {
			printf("Q[%d]%s =", k + 1, k >= m->M ? " (OI)" : "");
			mpq_set_ui(tot, 0, 1);
			for (r = 0; r < R; r++) {
				mpq_add(tot, tot, Q[k][r]);
				if (dbg) gmp_printf("\t%Qd", Q[k][r]); else { printf("\t"); put_d(Q[k][r]); }
			}
			if (dbg) gmp_printf("\t(total: %Qd)\n", tot); else { printf("\t(total: "); put_d(tot); printf(")\n"); }
		}
		mpq_clear(tot);
		printf("=========================================\n");
		printf("Elapsed time (OI convolution): %.6f s\n", cputime() - t0);
	}
	return bad ? 2 : 0;
}
