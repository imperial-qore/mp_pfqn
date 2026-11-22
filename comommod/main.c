/*
 * comommod - multi-modular CoMoM.
 *
 * The exact rational CoMoM solver (bin/comom) carries normalizing
 * constants whose length grows as Theta(N log N) bits, so the cost of
 * every arithmetic operation in the class-oriented recursion grows with
 * the population already processed.  Here the whole recursion is instead
 * run independently in Z/p_i for a sequence of word-size primes, and the
 * rational answers are recovered at the end by incremental CRT plus
 * rational reconstruction.
 *
 * Consequences, as for mommod:
 *   - no operand ever grows: every inner-loop operation is a single
 *     machine-word multiply-and-reduce;
 *   - the primes are completely independent, so the growth in total work
 *     is absorbed by parallelism rather than by wall clock;
 *   - memory per image is O(cardG*M) words rather than
 *     O(cardG*M * N log N) bits.
 *
 * Correctness caveats handled here:
 *   - p must exceed max_r N_r so that the division by n_r in the
 *     recursion is invertible;
 *   - a prime for which a BTF diagonal block is singular falls back to a
 *     direct LU of A11, and a prime for which A11 itself is singular
 *     ("unlucky prime") is discarded and another is drawn.  A matrix
 *     singular over Q is singular modulo every prime and shows up as a
 *     run of consecutive failures; the model is then solved exactly by
 *     the tiers of util/singular.h;
 *   - the sub-system solved on models with unreachable basis populations
 *     is rank-revealing, so a prime at which its rank drops is discarded
 *     rather than allowed to pick a different particular solution;
 *   - the number of primes is not fixed in advance.  Primes are added
 *     until rational reconstruction of every exported quantity succeeds
 *     and repeats unchanged.
 *
 * Unlike bin/comom this solver never returns a perturbed approximation:
 * a singular system is solved exactly through other class orders, merged
 * stations or classes, or per-prime interpolation of the exact limit of
 * a perturbation, as in bin/mommod.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include <time.h>
#include <zpla.h>
#include "comommod.h"
#include "singular.h"

/* Number of consecutive prime counts that must reproduce the same
 * reconstructed rational for G before the full set is reconstructed. */
#define STABLE_ROUNDS 1
#ifdef _OPENMP
#include <omp.h>
#endif

/* globals expected by profiling.h and the shared sources */
int AORSCTR;
int MULCTR;
int DIVCTR;
double AORSTIME;
double MULTIME;
double DIVTIME;
double t0,t1;
struct rusage ruse;

qnmodel* qnm;
bool INTERACTIVE,RANDGEN,CANON,ZSCALE,DEBUG,VERBOSE;
int VERBOSE_MOD = 0;

/* ---- reconstruction --------------------------------------------------
 *
 * Plain Wang splits the modulus evenly between numerator and denominator,
 * so it needs 2B bits of modulus to recover a B-bit value even when the
 * denominator is only a few bits wide, which is the case for every
 * quantity exported here: the CoMoM recursion divides only by the class
 * populations.  Bounding the denominator recovers the value from about B
 * bits and so from roughly half as many primes.
 *
 * A bounded reconstruction that succeeds is not necessarily correct,
 * so it is accepted only when it also matches every witness prime: a
 * prime that is solved but deliberately held out of the CRT accumulator,
 * whose residue therefore carries information the reconstruction cannot
 * have used.  A wrong candidate passes one witness with probability about
 * 1/p.  This is off by default; the default path is plain Wang, whose
 * size condition 2|num|den < M guarantees uniqueness on its own. */
static mpz_t DENBOUND;
static int   DENBITS_SET = 0;   /* --denbits given explicitly */
static int   USE_BOUNDED = 0;
static int   NWITNESS    = 0;
static zp_t  WPRIME[4];
static zp_t* WRES[4];

static int witness_ok(const mpz_t num, const mpz_t den, int idx)
{
	int w;
	for (w=0; w<NWITNESS; w++)
	{
		zp_t p = WPRIME[w];
		zp_t n = (zp_t) mpz_fdiv_ui((mpz_ptr)num, (unsigned long) p);
		zp_t d = (zp_t) mpz_fdiv_ui((mpz_ptr)den, (unsigned long) p);
		zp_t rhs = (zp_t)(((__uint128_t)WRES[w][idx] * (__uint128_t)d) % (__uint128_t)p);
		if (rhs != n) return 0;
	}
	return 1;
}

/* Denominator bound used by -B.
 *
 * mommod can bound the denominators of its exported constants by a fixed
 * few dozen bits, because mdecrease divides only by populations, station
 * multiplicities and think times.  CoMoM has no such bound: every
 * population step divides by n_r and every basis solve divides by the
 * pivots of A11, so the denominator grows linearly with the numerator.
 * Measured on the models here it runs at roughly a fifth of the length of
 * G, so a fixed 64-bit bound never fires and -B degenerates into plain
 * Wang.  Unless --denbits says otherwise the bound is therefore taken as
 * a third of the modulus accumulated so far: enough room for a
 * denominator of that shape, while still leaving the numerator two thirds
 * of the modulus instead of the half plain Wang leaves it.  The split is
 * a heuristic, not a theorem, which is exactly what the witness primes
 * are there to check. */
static void set_denbound(const mpz_t M)
{
	if (DENBITS_SET) return;
	mpz_ui_pow_ui(DENBOUND, 2, (unsigned long) (mpz_sizeinbase(M,2)/3 + 1));
}

static int reconstruct(mpz_t num, mpz_t den, const mpz_t c, const mpz_t M, int idx)
{
	if (USE_BOUNDED && NWITNESS > 0)
	{
		set_denbound(M);
		if (zp_ratrecon_bounded(num, den, c, M, DENBOUND)
		    && witness_ok(num, den, idx))
			return 1;
	}
	return zp_ratrecon(num, den, c, M);
}

/* Draw the next prime below cand, descending.  Returns 0 once the
 * sequence would fall to or below the floor, i.e. once the primes of the
 * requested size are exhausted.  Without this guard the search runs past
 * 2 into negative candidates, whose absolute values re-enter the sequence
 * as primes already used; a repeated prime makes M mod p zero, so
 * zp_crt_step adds no information while M keeps growing and the
 * reconstruction never converges. */
static int next_prime_below(mpz_t cand, unsigned long floorval, zp_t* out)
{
	for (;;)
	{
		mpz_sub_ui(cand, cand, 1);
		if (mpz_cmp_ui(cand, floorval) <= 0) return 0;
		if (mpz_probab_prime_p(cand, 25) != 0) { *out = (zp_t) mpz_get_ui(cand); return 1; }
	}
}

/* ------------------------------------------------------------------ */
static double to_double(mpz_t num, mpz_t den)
{
	/* Same conversion path as bin/comom (mpq -> mpf -> double) so that the
	 * two solvers' printed output is byte-identical and can be diffed. */
	double d;
	mpq_t q; mpf_t f;
	mpq_init(q); mpf_init(f);
	mpz_set(mpq_numref(q), num);
	mpz_set(mpq_denref(q), den);
	mpq_canonicalize(q);
	mpf_set_q(f, q);
	d = mpf_get_d(f);
	mpf_clear(f); mpq_clear(q);
	return d;
}

static double log_mpq(mpz_t num, mpz_t den)
{
	/* Take the same route as bin/comom (mpq -> mpf -> double, then log) so
	 * that the printed value is byte-identical to every other solver in
	 * the tree; the two-term form below rounds twice and can land one ulp
	 * away.  It is kept as the fallback for a G that no double can hold,
	 * where comom itself prints inf: mpz_get_d_2exp splits off the binary
	 * exponent, so the logarithm survives an overflowing ratio. */
	double d = to_double(num, den);
	if (d > 0.0 && !isinf(d)) return log(d);
	{
		signed long en, ed;
		double dn = mpz_get_d_2exp(&en, num);
		double dd = mpz_get_d_2exp(&ed, den);
		return log(dn) - log(dd) + (double)(en - ed) * log(2.0);
	}
}

/* Form a performance measure from reconstructed normalizing constants.
 * idx indexes the numerator constant; station>0 additionally applies the
 * queue-length weights m_k * L_kr, station<0 gives a throughput.  The
 * division by G is exact and happens once, here. */
static void ratio_exact(int idx, int station, mpz_t* num, mpz_t* den, mpq_t out)
{
	mpq_t b;
	mpq_init(b);
	mpz_set(mpq_numref(out), num[idx]); mpz_set(mpq_denref(out), den[idx]);
	mpq_canonicalize(out);
	mpz_set(mpq_numref(b), num[0]);     mpz_set(mpq_denref(b), den[0]);
	mpq_canonicalize(b);
	mpq_div(out, out, b);
	if (station > 0)
	{
		int r = (idx - 1 - qnm->R) % qnm->R + 1;
		if (mpz_cmp_ui(qnm->L[station-1][r-1], 0) == 0)
			mpq_set_ui(out, 0, 1);   /* a zero demand contributes no queue length */
		else
		{
			mpq_t w; mpq_init(w);
			mpz_set(mpq_numref(w), qnm->L[station-1][r-1]);
			mpz_set_ui(mpq_denref(w), 1);
			mpq_mul(out, out, w);
			mpq_set_si(w, qnm->mi[station-1], 1);
			mpq_mul(out, out, w);
			mpq_clear(w);
		}
	}
	mpq_clear(b);
}

static double ratio_double(int idx, int station, mpz_t* num, mpz_t* den)
{
	double d;
	mpq_t a; mpf_t f;
	mpq_init(a);
	ratio_exact(idx, station, num, den, a);
	mpf_init(f); mpf_set_q(f, a); d = mpf_get_d(f); mpf_clear(f);
	mpq_clear(a);
	return d;
}

/* ---- the CRT loop, over any per-prime image ---- */
static int S_NTHREADS = 1, S_PRIMEBITS = 61, S_NWITNESS = 0;
static long S_TOTPRIMES = 0, S_TOTUNLUCKY = 0;
static unsigned long S_MODBITS = 0;

typedef int (*image_fn)(qnmodel* q, zp_t p, zp_t* res, void* ctx);

/* The basis of each class depends only on (r,M), so it is built once per
 * model and shared read-only by every image. */
static combsrep** dn_build(const qnmodel* q)
{
	int r;
	combsrep** Dn = (combsrep**) calloc((size_t)q->R+1, sizeof(combsrep*));
	for (r=1;r<=q->R;r++)
	{
		Dn[r] = makeDn(r, q->M);
		if (Dn[r] == NULL) { fprintf(stderr,"Error: cannot build the class %d basis\n", r); exit(1); }
	}
	return Dn;
}

static void dn_free(combsrep** Dn, int R) { int r; for (r=1;r<=R;r++) freeDn(Dn[r]); free(Dn); }

static int img_direct(qnmodel* q, zp_t p, zp_t* res, void* ctx)
{
	return comom_image(q, (combsrep**) ctx, p, res);
}

/* Reconstruct the NRES(q) exported constants of q as exact rationals in
 * out[].  minprime is a lower bound on the primes (besides max N_r).
 * Returns 0 on success, 1 if the image failed modulo more than 40
 * consecutive primes (the system is singular over Q), 2 if the primes of
 * the requested size ran out (reported here). */
static int crt_solve(qnmodel* q, image_fn image, void* ctx, unsigned long minprime, mpq_t* out)
{
	int i, w, rc = 0, nres = NRES(q), nthreads = S_NTHREADS;
	unsigned long primefloor;
	mpz_t M, cand, pz, *acc, *num, *den, *pnum, *pden, *cnum, *cden;

	{
		int maxN = 0;
		for (i=0;i<q->R;i++) if (q->N[i] > maxN) maxN = q->N[i];
		primefloor = (unsigned long)(maxN > 2 ? maxN + 1 : 3);
		if (minprime > primefloor) primefloor = minprime;
	}
	mpz_init(cand);
	mpz_ui_pow_ui(cand, 2, (unsigned long)S_PRIMEBITS);
	if (mpz_cmp_ui(cand, primefloor) <= 0)
	{
		fprintf(stderr,"Error: -b %d gives no prime above %lu; raise -b\n", S_PRIMEBITS, primefloor);
		mpz_clear(cand);
		return 2;
	}

	acc  = (mpz_t*) malloc(nres*sizeof(mpz_t));
	cnum = (mpz_t*) malloc(nres*sizeof(mpz_t));
	cden = (mpz_t*) malloc(nres*sizeof(mpz_t));
	num  = (mpz_t*) malloc(nres*sizeof(mpz_t));
	den  = (mpz_t*) malloc(nres*sizeof(mpz_t));
	pnum = (mpz_t*) malloc(nres*sizeof(mpz_t));
	pden = (mpz_t*) malloc(nres*sizeof(mpz_t));
	mpz_init_set_ui(M,0); mpz_init(pz);
	for (i=0;i<nres;i++) { mpz_init(acc[i]); mpz_init(num[i]); mpz_init(den[i]); mpz_init(pnum[i]); mpz_init(pden[i]); mpz_init(cnum[i]); mpz_init(cden[i]); }

	/* Witness primes: solved like any other image but never folded into
	 * the CRT accumulator.  Drawn first from the same descending
	 * sequence, hence disjoint from the CRT primes. */
	NWITNESS = 0;
	for (w=0; w<S_NWITNESS; w++) WRES[w] = (zp_t*) malloc((size_t)nres*sizeof(zp_t));
	for (w=0; w<S_NWITNESS; w++)
	{
		int fails = 0;
		for (;;)
		{
			if (!next_prime_below(cand, primefloor, &WPRIME[w]))
			{
				fprintf(stderr,"Error: exhausted the primes of %d bits while drawing"
				               " witnesses; raise -b\n", S_PRIMEBITS);
				rc = 2; goto finish;
			}
			if (image(q, WPRIME[w], WRES[w], ctx) == 0) break;
			if (++fails > 40) { rc = 1; goto finish; }
		}
	}
	NWITNESS = S_NWITNESS;

	{
	int nprimes=0, unlucky=0, consec_fail=0, done=0;
	int next_try=4, gstable=0, confirmed=0;
	int batch = nthreads>1 ? nthreads : 1;
	zp_t* plist = (zp_t*) malloc(batch*sizeof(zp_t));
	zp_t* rbuf  = (zp_t*) malloc((size_t)batch*nres*sizeof(zp_t));
	int*  ok    = (int*)  malloc(batch*sizeof(int));

	while (!done && rc == 0)
	{
		int nb, bi;
		for (nb=0; nb<batch; nb++)
			if (!next_prime_below(cand, primefloor, &plist[nb]))
			{
				fprintf(stderr,
				        "Error: exhausted the primes of %d bits after %d of them"
				        " (%lu bits of modulus).\n"
				        "The reconstruction needs a larger modulus than this prime size"
				        " can supply; rerun with a larger -b\n"
				        "(the default -b 61 supplies on the order of 1e18 bits and is"
				        " never exhausted in practice).\n",
				        S_PRIMEBITS, nprimes, (unsigned long)mpz_sizeinbase(M,2));
				rc = 2; break;
			}
		if (rc) break;

		/* The first batch runs sequentially: nck() memoizes lazily into a
		 * shared table and is not thread-safe until it is warm, but is
		 * read-only thereafter. */
#ifdef _OPENMP
		#pragma omp parallel for num_threads(nthreads) schedule(dynamic,1) if(batch>1 && nprimes>0)
#endif
		for (bi=0; bi<batch; bi++)
			ok[bi] = (image(q, plist[bi], rbuf + (size_t)bi*nres, ctx) == 0);

		for (bi=0; bi<batch && !done; bi++)
		{
			if (!ok[bi]) { unlucky++; consec_fail++;
				if (consec_fail > 40) { rc = 1; break; }
				continue;
			}
			consec_fail = 0;
			nprimes++;
			for (i=0;i<nres;i++) zp_crt_step(acc[i], M, rbuf[(size_t)bi*nres+i], plist[bi]);
			mpz_set_ui(pz, (unsigned long) plist[bi]);
			if (mpz_sgn(M)==0) mpz_set(M, pz); else mpz_mul(M, M, pz);

			/* Reconstruction is quadratic in the size of the modulus, so
			 * attempting it after every prime would dominate the run.  Use
			 * a geometric schedule and probe G alone; only once G has
			 * settled is the full set reconstructed, and it must then
			 * reproduce itself at a later prime count to be accepted. */
			if (nprimes < next_try) continue;
			next_try = (int)(nprimes * 1.12) + 1;

			if (!reconstruct(num[0], den[0], acc[0], M, 0)) { gstable = 0; continue; }
			if (mpz_cmp(num[0],pnum[0])==0 && mpz_cmp(den[0],pden[0])==0) gstable++;
			else { gstable = 0; mpz_set(pnum[0],num[0]); mpz_set(pden[0],den[0]); }
			if (gstable < STABLE_ROUNDS) continue;

			{
				int allok = 1, same = 1;
				for (i=0;i<nres;i++)
				{
					if (!reconstruct(num[i], den[i], acc[i], M, i)) { allok = 0; break; }
					if (mpz_cmp(num[i],cnum[i])!=0 || mpz_cmp(den[i],cden[i])!=0) same = 0;
				}
				if (!allok) { gstable = 0; confirmed = 0; continue; }
				if (confirmed && same) { done = 1; break; }
				for (i=0;i<nres;i++) { mpz_set(cnum[i],num[i]); mpz_set(cden[i],den[i]); }
				confirmed = 1;
			}
		}
	}
	S_TOTPRIMES += nprimes; S_TOTUNLUCKY += unlucky;
	if (done)
	{
		S_MODBITS = (unsigned long) mpz_sizeinbase(M,2);
		for (i=0;i<nres;i++)
		{
			mpz_set(mpq_numref(out[i]), num[i]);
			mpz_set(mpq_denref(out[i]), den[i]);
			mpq_canonicalize(out[i]);
		}
	}
	free(plist); free(rbuf); free(ok);
	}

finish:
	for (w=0; w<S_NWITNESS; w++) free(WRES[w]);
	NWITNESS = 0;
	for (i=0;i<nres;i++) { mpz_clear(acc[i]); mpz_clear(num[i]); mpz_clear(den[i]); mpz_clear(pnum[i]); mpz_clear(pden[i]); mpz_clear(cnum[i]); mpz_clear(cden[i]); }
	free(acc); free(num); free(den); free(pnum); free(pden); free(cnum); free(cden);
	mpz_clear(M); mpz_clear(pz); mpz_clear(cand);
	return rc;
}

/* ---- singular models: see util/singular.h ---- */
static int sg_direct(qnmodel* q, mpq_t* out, void* ctx)
{
	int rc;
	combsrep** Dn = dn_build(q);
	(void)ctx;
	rc = crt_solve(q, img_direct, Dn, 0, out);
	dn_free(Dn, q->R);
	return rc;
}

/* Tier 4 per prime, as in mommod: each prime solves the perturbed models
 * at t=1,2,..., skips the t that fail, and evaluates the Lagrange
 * interpolant at t=0 in Z/p; the CRT then reconstructs g(0) exactly.  The
 * perturbed models share M and R with q, hence its class bases. */
typedef struct { int D, T; qnmodel** qt; combsrep** Dn; } interp_ctx;

static zp_t mulm(zp_t a, zp_t b, zp_t p) { return (zp_t)(((__uint128_t)a * b) % p); }
static zp_t powm(zp_t a, zp_t e, zp_t p) { zp_t r = 1; while (e) { if (e & 1) r = mulm(r, a, p); a = mulm(a, a, p); e >>= 1; } return r; }

static int img_interp(qnmodel* q, zp_t p, zp_t* res, void* vctx)
{
	interp_ctx* c = (interp_ctx*) vctx;
	int nres = NRES(q), D = c->D, got = 0, t, i, j, u;
	zp_t* buf = (zp_t*) malloc((size_t)(D+1)*nres*sizeof(zp_t));
	int* ts = (int*) malloc((D+1)*sizeof(int));
	for (t=1; t<=c->T && got<D+1; t++)
		if (comom_image(c->qt[t-1], c->Dn, p, buf + (size_t)got*nres) == 0) ts[got++] = t;
	if (got < D+1) { free(buf); free(ts); return -1; }
	for (u=0;u<nres;u++) res[u] = 0;
	for (i=0;i<=D;i++)
	{
		zp_t numw = 1, denw = 1;
		for (j=0;j<=D;j++) if (j != i)
		{
			numw = mulm(numw, p - (zp_t)ts[j], p);                         /* (0 - t_j) */
			denw = mulm(denw, ((zp_t)ts[i] + p - (zp_t)ts[j]) % p, p);     /* (t_i - t_j) */
		}
		zp_t wgt = mulm(numw, powm(denw, p-2, p), p);
		for (u=0;u<nres;u++) res[u] = (res[u] + mulm(wgt, buf[(size_t)i*nres+u], p)) % p;
	}
	free(buf); free(ts);
	return 0;
}

static int sg_interp_mod(qnmodel* q, mpq_t* out, void* vctx)
{
	interp_ctx c;
	int rc;
	(void)vctx;
	c.D = sg_interp_degree(q);
	c.T = sg_interp_points(q);
	c.qt = sg_perturb_build(q, c.T);
	c.Dn = dn_build(q);
	rc = crt_solve(q, img_interp, &c, (unsigned long)c.T + 1, out);
	dn_free(c.Dn, q->R);
	sg_perturb_free(c.qt, c.T);
	return rc;
}

static void usage(const char* prog)
{
	printf("USAGE: %s [-l|--log] [-e|--ex] [-g|--nc] [-t|--tput] [-q|--qlen] [-d|--debug]\n", prog);
	printf("       %*s [-v] [-j nthreads] [-b bits] [-B] [-W n] [--denbits n] model.qn\n",
	       (int)strlen(prog), "");
	printf("  -l, --log        : Print only log of normalizing constant as double\n");
	printf("  -e, --ex         : Print exact normalizing constant numerator and denominator\n");
	printf("  -g, --nc         : Print normalizing constant as double\n");
	printf("  -t, --tput       : Print only throughputs, one per row\n");
	printf("  -q, --qlen       : Print only queue lengths, one per row\n");
	printf("  -d, --debug      : Print performance metrics as exact rational numbers\n");
	printf("  -v               : Report prime count and timing on stderr\n");
	printf("  -j n             : Number of primes solved concurrently (OpenMP)\n");
	printf("  -b bits          : Bit size of the primes drawn (default 61)\n");
	printf("  -B, --bounded    : Bound the denominator during rational reconstruction.\n");
	printf("                     Recovers a B-bit constant from about B bits of modulus\n");
	printf("                     instead of 2B, roughly halving the number of primes.\n");
	printf("                     UNSOUND ON ITS OWN, so it is off by default and is only\n");
	printf("                     accepted when it also matches every witness prime.\n");
	printf("  -W, --witness n  : Number of primes solved but held out of the CRT, used\n");
	printf("                     solely to verify bounded reconstructions (default 1\n");
	printf("                     when -B is given, 0 otherwise).\n");
	printf("      --denbits n  : Denominator bound for -B, as a power of two.  The default\n");
	printf("                     is a third of the modulus accumulated so far, since CoMoM\n");
	printf("                     denominators grow with the population rather than staying\n");
	printf("                     a few dozen bits wide as MoM's do.\n");
}

int main(int argc, char**argv)
{
	int i, r, k, nres, nthreads = 1, primebits = 61;
	int nwitness = 0, denbits = 64;
	bool log_output=false, normconst_output=false, normconst_g_output=false;
	bool throughput_output=false, queue_output=false, debug_output=false;
	char* model_file=NULL;

	t0=CPUTIME;
	DEBUG=false; VERBOSE=true;

	if (argc < 2) { usage(argv[0]); return -1; }
	for (i=1;i<argc;i++)
	{
		if      (!strcmp(argv[i],"-l")||!strcmp(argv[i],"--log"))   log_output=true;
		else if (!strcmp(argv[i],"-e")||!strcmp(argv[i],"--ex"))    normconst_output=true;
		else if (!strcmp(argv[i],"-g")||!strcmp(argv[i],"--nc"))    normconst_g_output=true;
		else if (!strcmp(argv[i],"-t")||!strcmp(argv[i],"--tput"))  throughput_output=true;
		else if (!strcmp(argv[i],"-q")||!strcmp(argv[i],"--qlen"))  queue_output=true;
		else if (!strcmp(argv[i],"-d")||!strcmp(argv[i],"--debug")) debug_output=true;
		else if (!strcmp(argv[i],"-v"))                             VERBOSE_MOD=1;
		else if (!strcmp(argv[i],"-h")||!strcmp(argv[i],"--help"))  { usage(argv[0]); return 0; }
		else if (!strcmp(argv[i],"-j")) { if(i+1<argc) nthreads=atoi(argv[++i]); }
		else if (!strcmp(argv[i],"-b")) { if(i+1<argc) primebits=atoi(argv[++i]); }
		else if (!strcmp(argv[i],"-B")||!strcmp(argv[i],"--bounded")) { USE_BOUNDED=1; if(nwitness==0) nwitness=1; }
		else if (!strcmp(argv[i],"-W")||!strcmp(argv[i],"--witness")) { if(i+1<argc) nwitness=atoi(argv[++i]); }
		else if (!strcmp(argv[i],"--denbits")) { if(i+1<argc) { denbits=atoi(argv[++i]); DENBITS_SET=1; } }
		else if (argv[i][0]=='-') { /* ignore */ }
		else model_file=argv[i];
	}
	if (!model_file) { usage(argv[0]); return -1; }
	if (primebits < 8 || primebits > 62) { fprintf(stderr,"prime bit size must be in [8,62]\n"); return -1; }

	qnm = (qnmodel*)readmodel(model_file);
	qn_require_closed(qnm, "comommod", 0);
	nckinit(qnm->M+qnm->R+1,qnm->M);
	nres = NRES(qnm);

	mpz_init(DENBOUND);
	mpz_ui_pow_ui(DENBOUND, 2, (unsigned long) denbits);
	if (nwitness > 4) nwitness = 4;
	if (!USE_BOUNDED) nwitness = 0;
	S_NTHREADS = nthreads; S_PRIMEBITS = primebits; S_NWITNESS = nwitness;

	mpz_t *num = (mpz_t*) malloc(nres*sizeof(mpz_t));
	mpz_t *den = (mpz_t*) malloc(nres*sizeof(mpz_t));
	mpq_t *exact = sg_mpq_arr(nres);
	char tier[256];
	double tstart = CPUTIME;
	sg_ops ops = { sg_direct, sg_interp_mod, NULL, " per prime" };
	int rc = sg_solve_tiered(qnm, exact, tier, sizeof tier, &ops, 0);
	if (rc == 1)
	{
		fprintf(stderr,"Error: the CoMoM system is singular for this model, and every exact fallback\n"
		               "failed (other class orders, merged stations or classes, interpolation).\n"
		               "The normalizing constant may be zero or undefined; bin/ca solves any\n"
		               "well-posed model by convolution.\n");
		return 1;
	}
	if (rc != 0) return 1;
	for (i=0;i<nres;i++)
	{
		mpz_init_set(num[i], mpq_numref(exact[i]));
		mpz_init_set(den[i], mpq_denref(exact[i]));
	}
	sg_mpq_arr_free(exact, nres);

	double telapsed = CPUTIME - tstart;
	if (VERBOSE_MOD)
		fprintf(stderr,"comommod: %ld primes of %d bits (%ld unlucky discarded), modulus %lu bits, %.6f s, %d thread(s)%s; route: %s\n",
		        S_TOTPRIMES, primebits, S_TOTUNLUCKY, S_MODBITS, telapsed, nthreads,
		        USE_BOUNDED ? " [bounded reconstruction, witness-verified]" : "", tier);
	else if (strcmp(tier, "direct") != 0)
		fprintf(stderr,"comommod: singular CoMoM system, solved exactly via %s\n", tier);

	/* ---- report, in the same formats as bin/comom ---- */
	if (log_output)
		printf("%.15e\n", log_mpq(num[0],den[0]));
	else if (normconst_output)
	{
		gmp_printf("%Zd\n", num[0]);
		gmp_printf("%Zd\n", den[0]);
	}
	else if (normconst_g_output)
		printf("%.15e\n", to_double(num[0],den[0]));
	else if (throughput_output)
	{
		for (r=1;r<=qnm->R;r++) printf("%.15e\n", ratio_double(r,-1,num,den));
	}
	else if (queue_output)
	{
		for (k=1;k<=qnm->M;k++)
		{
			for (r=1;r<=qnm->R;r++)
			{
				printf("%.15e", ratio_double(1+qnm->R+(k-1)*qnm->R+(r-1),k,num,den));
				if (r<qnm->R) printf(" ");
			}
			printf("\n");
		}
	}
	else
	{
		mpq_t v, tot;
		mpq_init(v); mpq_init(tot);
		printmodel(qnm);
		printf("========== Performance Metrics ==========\n");
		if (debug_output)
		{
			mpq_t Gq; mpq_init(Gq);
			mpz_set(mpq_numref(Gq), num[0]); mpz_set(mpq_denref(Gq), den[0]);
			mpq_canonicalize(Gq);
			gmp_printf("G = %Qd\n", Gq);
			mpq_clear(Gq);
		}
		else
		{
			printf("G = %.15e\n", to_double(num[0],den[0]));
			printf("log(G) = %.15e\n", log_mpq(num[0],den[0]));
		}
		printf("\nX (throughputs):\n");
		for (r=1;r<=qnm->R;r++)
		{
			if (debug_output)
			{
				ratio_exact(r,-1,num,den,v);
				gmp_printf("X[%d] = %Qd\n", r, v);
			}
			else
				printf("X[%d] = %.15e\n", r, ratio_double(r,-1,num,den));
		}
		printf("\nQ (mean queue lengths):\n");
		for (k=1;k<=qnm->M;k++)
		{
			printf("Q[%d] =", k);
			if (debug_output)
			{
				mpq_set_ui(tot,0,1);
				for (r=1;r<=qnm->R;r++)
				{
					ratio_exact(1+qnm->R+(k-1)*qnm->R+(r-1),k,num,den,v);
					mpq_add(tot,tot,v);
					printf("\t"); gmp_printf("%Qd", v);
				}
				printf("\t(total: "); gmp_printf("%Qd", tot); printf(")\n");
			}
			else
			{
				double sum=0.0;
				for (r=1;r<=qnm->R;r++)
				{
					double d = ratio_double(1+qnm->R+(k-1)*qnm->R+(r-1),k,num,den);
					sum += d;
					printf("\t%.15e", d);
				}
				printf("\t(total: %.15e)\n", sum);
			}
		}
		printf("=========================================\n");
		printf("Primes: %ld of %d bits (%ld unlucky), modulus %lu bits, %d thread(s)\n",
		       S_TOTPRIMES, primebits, S_TOTUNLUCKY, S_MODBITS, nthreads);
		if (strcmp(tier, "direct") != 0) printf("Singular CoMoM system, solved exactly via %s\n", tier);
		t1 = CPUTIME;
		printf("Elapsed time (CoMoM-mod): %.6f s\n", t1-t0);
		mpq_clear(v); mpq_clear(tot);
	}

	return 0;
}
