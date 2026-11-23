/* The multi-modular CoMoM engine: CRT plus rational reconstruction over
 * comom_image() for many primes, and the tiers of util/singular.h with
 * tier 4 interpolated per prime.  Split out of main.c so that bin/comom
 * can link it for its own interpolation tier (see ccengine.h). */
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <gmp.h>
#include <zpla.h>
#include "comommod.h"
#include "singular.h"
#include "gbound.h"
#include "ccengine.h"

/* Number of consecutive prime counts that must reproduce the same
 * reconstructed rational for G before the full set is reconstructed. */
#define STABLE_ROUNDS 1
#ifdef _OPENMP
#include <omp.h>
#endif


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
/* A-priori modulus (util/gbound.h): on by default; --adaptive turns it off
 * and keeps only the adaptive schedule below. */
static int S_APRIORI = 1;
void comommod_apriori(int on) { S_APRIORI = on; }

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
	/* With bounds G <= G_up and den | D known before the first prime, a
	 * modulus M > 2 G_up D^2 makes the bounded reconstruction of every
	 * quantity exact (util/gbound.h): no probing, no stability rounds, no
	 * confirmation at a later prime count, and about half the primes. */
	mpz_t* DB = (mpz_t*) malloc(nres*sizeof(mpz_t));
	long tbits = S_APRIORI ? gb_modbits(q, DB) : -1;
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
				        " can supply; rerun with a larger --primebits\n"
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

			if (tbits >= 0) {
				int allok = 1;
				if ((long)mpz_sizeinbase(M,2) <= tbits) continue;
				for (i=0;i<nres && allok;i++)
					if (!zp_ratrecon_bounded(num[i], den[i], acc[i], M, DB[i]) || !witness_ok(num[i], den[i], i)) allok = 0;
				if (allok) { done = 1; break; }
				/* cannot happen for a valid bound; the adaptive loop takes over */
				fprintf(stderr, "comommod: a-priori bound of %ld bits did not reconstruct; continuing adaptively\n", tbits);
				tbits = -1;
				continue;
			}

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
	if (S_APRIORI) { for (i=0;i<nres;i++) mpz_clear(DB[i]); }
	free(DB);
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

/* ---- API (ccengine.h) ---- */
void comommod_config(int nthreads, int primebits, int nwitness, int use_bounded, int denbits)
{
	static int init = 0;
	if (!init) { mpz_init(DENBOUND); init = 1; }
	if (denbits > 0) { mpz_ui_pow_ui(DENBOUND, 2, (unsigned long) denbits); DENBITS_SET = 1; }
	else DENBITS_SET = 0;
	if (nwitness > 4) nwitness = 4;
	USE_BOUNDED = use_bounded;
	if (!use_bounded) nwitness = 0;
	S_NTHREADS = nthreads > 0 ? nthreads : 1;
	S_PRIMEBITS = primebits; S_NWITNESS = nwitness;
}

void comommod_stats(long* primes, long* unlucky, unsigned long* modbits)
{
	if (primes) *primes = S_TOTPRIMES;
	if (unlucky) *unlucky = S_TOTUNLUCKY;
	if (modbits) *modbits = S_MODBITS;
}

int comommod_solve_tiered(qnmodel* q, mpq_t* out, char* tier, size_t tsz)
{
	sg_ops ops = { sg_direct, sg_interp_mod, NULL, " per prime" };
	return sg_solve_tiered(q, out, tier, tsz, &ops, 0);
}

int comommod_interp(qnmodel* q, mpq_t* out, void* ctx)
{
	/* the per-prime recursion reads mi as an integer copy count */
	if (qn_mden(q) > 1) {
		fprintf(stderr, "comommod_interp: fractional multiplicities are not supported by the modular interpolation tier\n");
		return 1;
	}
	return sg_interp_mod(q, out, ctx);
}
