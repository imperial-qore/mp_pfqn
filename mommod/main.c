/*
 * mommod - multi-modular Method of Moments.
 *
 * The exact rational solver (bin/mom) carries normalizing constants whose
 * length grows as Theta(N log N) bits, so the cost of every arithmetic
 * operation in the population recursion grows with the population already
 * processed.  Here the whole recursion is instead run independently in
 * Z/p_i for a sequence of word-size primes, and the rational answers are
 * recovered at the end by incremental CRT plus rational reconstruction.
 *
 * Consequences:
 *   - no operand ever grows: every inner-loop operation is a single
 *     machine-word multiply-and-reduce, so the per-population cost is
 *     constant in N instead of growing with it;
 *   - the primes are completely independent, so the growth in total work
 *     is absorbed by parallelism rather than by wall clock;
 *   - memory per image is O(theta) words rather than O(theta * N log N)
 *     bits.
 *
 * Correctness caveats handled here:
 *   - p must exceed max_r N_r so that the division by n_r in the
 *     recursion is invertible;
 *   - a prime for which some diagonal block of C is singular ("unlucky
 *     prime") is discarded and another is drawn.  A genuinely singular C
 *     is singular modulo every prime and shows up as a run of consecutive
 *     failures; the model is then solved exactly by solve_tiered()
 *     (other class orders, merged stations or classes, interpolation);
 *   - the number of primes is not fixed in advance.  Primes are added
 *     until rational reconstruction of every exported quantity succeeds
 *     and repeats unchanged, which is the standard early-termination
 *     criterion and avoids needing an a priori height bound.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include <time.h>
#include <zpla.h>
#include "mommod.h"
#include "lvtab.h"

/* Number of consecutive primes that must reproduce the same reconstructed
 * rational before it is accepted.  One extra prime is the textbook
 * criterion; two makes a spurious agreement negligibly unlikely. */
#define STABLE_ROUNDS 1
#ifdef _OPENMP
#include <omp.h>
#endif

/* globals expected by profiling.h and the shared solver sources */
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

/* Lower replica levels: see util/lvtab.h. */
/* Advance levels 0..r-2 from population n-1_r to n.  glow points at the
 * old level r-1 (second part of g), which must not yet be overwritten. */
static void lower_step(qnmodel* qnm, lvtab* T, int r, zp_vec_t* low, zp_t* glow, zp_t ninv)
{
	int l, i, j, s, m = T->m;
	zp_t zr = zp_from_z(qnm->Z[r-1]);
	zp_t* Lr = (zp_t*)malloc(m*sizeof(zp_t));
	for (j=0;j<m;j++) Lr[j] = zp_from_z(qnm->L[j][r-1]);
	for (l=0;l<=r-2;l++) {
		zp_t* above = (l+1 == r-1) ? glow : low[l+1];
		for (i=0;i<T->card[l];i++) {
			for (s=0;s<r;s++) {
				zp_t acc = zp_mulv(zr, low[l][i*r+s]);
				for (j=0;j<m;j++) {
					if (Lr[j] == 0) continue;
					zp_t c = zp_mulv(Lr[j], zp_from_si(qnm->mi[j] + T->cb[l][i][j]));
					acc = zp_addv(acc, zp_mulv(c, above[T->up[l][i*m+j]*r+s]));
				}
				low[l][i*r+s] = zp_mulv(ninv, acc);
			}
		}
	}
	free(Lr);
}

/* ------------------------------------------------------------------ */
/* One modular image of the whole MoM recursion.                       */
/* Returns 0 on success, -1 if the prime is unlucky (or C is singular). */
/* ------------------------------------------------------------------ */
static int solve_mod(qnmodel* qnm, zp_t p, zp_t* res)
{
	int t,i,r,s,h,k;
	int *n;
	int cardG, cardGk;
	zp_vec_t g=NULL,Gk=NULL,G=NULL,b=NULL,b2=NULL,gr=NULL;
	LS* A=NULL;
	int m=qnm->M;
	int rc = 0;
	int R=qnm->R, l;
	lvtab* T=NULL;
	int full;
	zp_vec_t *low=NULL, *lowr=NULL; /* levels 0..r-2 and their snapshot at N-1_r */

	zp_setmod(p);
	full = mom_mixed_think(qnm); /* else the mdecrease descent is exact */
	T = full ? lvtab_new(m,R) : NULL;
	low  = (zp_vec_t*)calloc(R>1?R-1:1,sizeof(zp_vec_t));
	lowr = (zp_vec_t*)calloc(R>1?R-1:1,sizeof(zp_vec_t));

	n=(int*)int_vec(qnm->R,0);

	for(r=1;r<=qnm->R;r++)
	{
		cardGk=nck(m+r-1,r)*r;
		cardG=nck(m+r-2,r-1)*r;

		if(r==1)
		{
			A = setupls(qnm->L,n,qnm->Z,qnm->mi,m,r);
			if (A == NULL) { rc = -1; goto cleanup; }
			b  = zp_vec(m,0,1);
			b2 = zp_vec(m,0,1);
			g  = zp_vec(m+1,1,1);
			gr = zp_vec(m+1,0,1); /* G(n-1_1)=0 at n=0: stays 0 if N_1=0 */
			Gk = zp_vec(m,1,1);
			G  = zp_vec(1,1,1);
		}
		else
		{
			freels(A);
			A = setupls(qnm->L,n,qnm->Z,qnm->mi,m,r);
			if (A == NULL) { rc = -1; goto cleanup; }
			/* seed levels 0..r-2 for class r: entries for classes
			 * 1..r-2 carry over, the class r-1 entry is G(n-1_{r-1})
			 * taken from the snapshot at N-1_{r-1} */
			{
				int cardGk_old = nck(m+r-2,r-1)*(r-1);
				for (l=0;full && l<=r-2;l++) {
					zp_t* old  = (l == r-2) ? g  + cardGk_old : low[l];
					zp_t* snap = (l == r-2) ? gr + cardGk_old : lowr[l];
					zp_vec_t nw = zp_vec(T->card[l]*r,0,1);
					for (i=0;i<T->card[l];i++) {
						for (s=0;s<=r-2;s++) nw[i*r+s] = old[i*(r-1)+s];
						nw[i*r+r-1] = snap[i*(r-1)];
					}
					if (l < r-2) { zp_vecfree(low[l]); zp_vecfree(lowr[l]); }
					low[l]  = nw;
					lowr[l] = zp_vec(T->card[l]*r,0,1); /* stays 0 if N_r=0 */
				}
			}
			zp_vecfree(b); zp_vecfree(b2); zp_vecfree(G);
			b  = zp_vec(cardGk,0,1);
			b2 = zp_vec(cardGk,0,1);
			G  = zp_vec(cardG,0,1);
			for (i=0;i<nck(m+r-2,r-1);i++)
			{
				for (s=0;s<=r-2;s++) G[i*r+s] = g[i*(r-1)+s];
				G[i*r+r-1] = gr[i*(r-1)];
			}
			zp_mspvecmul(b,A->A12,G);
			for(t=0;t<cardGk;t++) b[t] = zp_negv(b[t]);
			Gk = blocksolve(A,b);
			if (Gk == NULL) { rc = -1; goto cleanup; }
			zp_vecfree(g); zp_vecfree(gr);
			g  = zp_vec(cardGk+cardG,0,1);
			gr = zp_vec(cardGk+cardG,0,1);
			for(t=0;t<cardGk;t++) g[t] = Gk[t];
			for(t=0;t<cardG;t++)  g[cardGk+t] = G[t];
			zp_vecfree(Gk); Gk=NULL;
		}

		for (n[r-1]=1;n[r-1]<=qnm->N[r-1];n[r-1]++)
		{
			if(n[r-1]==qnm->N[r-1]) {
				zp_vecdup(gr,g,cardGk+cardG);
				for (l=0;full && l<=r-2;l++) zp_vecdup(lowr[l],low[l],T->card[l]*r);
			}

			zp_mspvecmul(G,A->B2r,g);
			{
				zp_t ninv = zp_invv(zp_from_si(n[r-1]));
				for(t=0;t<cardG;t++) G[t] = zp_mulv(ninv, G[t]);
				if (full) lower_step(qnm,T,r,low,g+cardGk,ninv);
			}

			zp_mspvecmul(b,A->B1r,g);
			zp_mspvecmul(b2,A->A12,G);
			for(t=0;t<cardGk;t++) b[t] = zp_subv(b[t],b2[t]);

			Gk = blocksolve(A,b);
			if (Gk == NULL) { rc = -1; goto cleanup; }

			for(t=0;t<cardGk;t++) g[t] = Gk[t];
			for(t=0;t<cardG;t++)  g[cardGk+t] = G[t];
			zp_vecfree(Gk); Gk=NULL;
		}
		n[r-1]=qnm->N[r-1];
	}

	if (ZP_SINGULAR) { rc = -1; goto cleanup; }
	if (!full) mdecrease(qnm,G,Gk,g,gr,res);
	else
	/* read levels 0 and 1 at N into the layout exportres expects:
	 * G = [G(N), G(N-1_1..R)], Gk[k] = [G^{+k}(N), G^{+k}(N-1_1..R)] */
	{
		zp_t *L0, *S0, *L1, *S1;
		int cardGkR = nck(m+R-1,R)*R;
		L0 = (R >= 2) ? low[0]  : g + cardGkR;
		S0 = (R >= 2) ? lowr[0] : gr + cardGkR;
		L1 = (R >= 3) ? low[1]  : (R == 2 ? g + cardGkR  : g);
		S1 = (R >= 3) ? lowr[1] : (R == 2 ? gr + cardGkR : gr);
		zp_vecfree(G); zp_vecfree(Gk);
		G  = zp_vec(R+1,0,1);
		Gk = zp_vec(m*(R+1),0,1);
		for (s=0;s<R;s++) G[s] = L0[s];
		G[R] = S0[0];
		for (k=0;k<m;k++) {
			for (s=0;s<R;s++) Gk[k*(R+1)+s] = L1[k*R+s];
			Gk[k*(R+1)+R] = S1[k*R];
		}
		exportres(qnm,G,Gk,res);
	}
	if (ZP_SINGULAR) rc = -1;
	else { int u; for (u=0;u<NRES(qnm);u++) res[u] = zp_from_mont(res[u]); }

cleanup:
	(void)h; (void)k;
	freels(A);
	zp_vecfree(b); zp_vecfree(b2); zp_vecfree(g); zp_vecfree(gr);
	zp_vecfree(G); zp_vecfree(Gk);
	if (low)  for (l=0;l<R-1;l++) zp_vecfree(low[l]);
	if (lowr) for (l=0;l<R-1;l++) zp_vecfree(lowr[l]);
	free(low); free(lowr);
	lvtab_free(T);
	free(n);
	return rc;
}


/* Recover a rational from its residue.  The exported normalizing
 * constants have denominators of at most a few dozen bits, so the bounded
 * reconstruction succeeds from about half as many primes as plain Wang;
 * the unbounded routine remains the fallback for anything whose
 * denominator turns out to be larger than the bound. */
/* Reconstruction of a rational from its CRT residue.
 *
 * Plain Wang splits the modulus evenly between numerator and denominator,
 * so it needs 2B bits of modulus to recover a B-bit value even when the
 * denominator is only a few bits wide, which is the case for every
 * quantity exported here: the recursion divides only by class
 * populations (lower levels) and the pivots of the population-independent
 * diagonal blocks.  Bounding the denominator recovers the
 * value from about B bits and so from roughly half as many primes.
 *
 * A bounded reconstruction that succeeds is not necessarily correct.
 * While the modulus is still too small for the true value to satisfy the
 * bounds, the Euclid nevertheless returns some pair that does, and for a
 * quantity whose denominator exceeds the bound it can keep returning
 * different spurious pairs indefinitely.  The bounded route is therefore
 * sound only in combination with an independent check, which is what the
 * witness primes below provide: one or more primes are solved but held
 * out of the CRT accumulator, so their residues carry information the
 * reconstruction cannot have used.  A wrong candidate passes a witness
 * with probability about 1/p, i.e. 2^-61 per prime per quantity.
 *
 * This is disabled by default (-B enables it).  The default path is plain
 * Wang, whose size condition 2|num|den < M guarantees uniqueness with no
 * auxiliary check at all. */
static mpz_t DENBOUND;
static int   USE_BOUNDED = 0;
static int   NWITNESS    = 0;
static zp_t  WPRIME[4];
static zp_t* WRES[4];

/* Check num == c*den (mod p_w) for every held-out prime, where c is the
 * witness residue of the same quantity. */
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

static int reconstruct(mpz_t num, mpz_t den, const mpz_t c, const mpz_t M, int idx)
{
	if (USE_BOUNDED && NWITNESS > 0
	    && zp_ratrecon_bounded(num, den, c, M, DENBOUND)
	    && witness_ok(num, den, idx))
		return 1;
	return zp_ratrecon(num, den, c, M);
}

/* A bounded reconstruction can return a spurious pair while the modulus is
 * still too small for the true value to satisfy the bounds, so no single
 * reconstruction is trusted.  Every quantity must reconstruct to the same
 * rational at two separate prime counts before the run terminates. */


/* Draw the next prime below cand, descending.  Returns 0 once the
 * sequence would fall to or below the floor, i.e. once the primes of the
 * requested size are exhausted.
 *
 * Exhaustion is unreachable at the default -b 61: there are about
 * 2^60/ln(2^61) = 2.7e16 primes in [2^60,2^61], enough for a modulus of
 * roughly 1.7e18 bits.  It becomes reachable only when -b is set small.
 * Without this guard the search runs past 2 into negative candidates,
 * whose absolute values re-enter the sequence as primes already used; a
 * repeated prime makes M mod p zero, so zp_crt_step adds no information
 * while M keeps growing, the residue stops being a valid CRT residue for
 * the claimed modulus, and reconstruction never converges. */
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
static double log_mpq(mpz_t num, mpz_t den)
{
	/* log(num/den) without ever forming num/den as a double, which would
	 * overflow: mpz_get_d_2exp splits off the binary exponent. */
	signed long en, ed;
	double dn = mpz_get_d_2exp(&en, num);
	double dd = mpz_get_d_2exp(&ed, den);
	return log(dn) - log(dd) + (double)(en - ed) * log(2.0);
}

static double to_double(mpz_t num, mpz_t den)
{
	/* Same conversion path as bin/mom (mpq -> mpf -> double) so that the
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


/* Form a performance measure from reconstructed normalizing constants.
 * idx indexes the numerator constant; station>0 additionally applies the
 * queue-length weights m_k * rho_kr, station<0 gives a throughput.  The
 * division by G is exact and happens once, here. */
static double ratio_double(int idx, int station, mpz_t* num, mpz_t* den)
{
	double d;
	mpq_t a, b;
	mpq_init(a); mpq_init(b);
	mpz_set(mpq_numref(a), num[idx]); mpz_set(mpq_denref(a), den[idx]);
	mpq_canonicalize(a);
	mpz_set(mpq_numref(b), num[0]);   mpz_set(mpq_denref(b), den[0]);
	mpq_canonicalize(b);
	mpq_div(a, a, b);
	if (station > 0)
	{
		int r = (idx - 1 - qnm->R) % qnm->R + 1;
		if (mpz_cmp_ui(qnm->L[station-1][r-1], 0) == 0)
			mpq_set_ui(a, 0, 1);   /* a zero demand contributes no queue length */
		else
		{
			mpq_t w; mpq_init(w);
			mpz_set(mpq_numref(w), qnm->L[station-1][r-1]);
			mpz_set_ui(mpq_denref(w), 1);
			mpq_mul(a, a, w);
			mpq_set_si(w, qnm->mi[station-1], 1);
			mpq_mul(a, a, w);
			mpq_clear(w);
		}
	}
	{
		mpf_t f; mpf_init(f);
		mpf_set_q(f, a);
		d = mpf_get_d(f);
		mpf_clear(f);
	}
	mpq_clear(a); mpq_clear(b);
	return d;
}

/* ================================================================== */
/* Exact solution of models whose MoM system is singular.             */
/*                                                                    */
/* crt_solve() runs the prime/CRT/reconstruction loop for any model    */
/* and any per-prime "image" function.  solve_tiered() tries, cheapest */
/* first and every tier exact:                                         */
/*   1 the model as given;                                             */
/*   2 other class orders (the system can be singular for one class    */
/*     processed last and not for another);                            */
/*   3 structural reductions, each solved again by tiers 1-3:          */
/*     stations with identical demand rows folded into one station     */
/*     with summed multiplicity, and classes with proportional demands */
/*     and think times aggregated into one class;                      */
/*   4 the exact limit of a perturbation, by interpolation in Z/p.     */
/* ================================================================== */

static int  S_NTHREADS = 1, S_PRIMEBITS = 61, S_NWITNESS = 0;
static long S_TOTPRIMES = 0, S_TOTUNLUCKY = 0;
static unsigned long S_MODBITS = 0;

typedef int (*image_fn)(qnmodel* q, zp_t p, zp_t* res, void* ctx);

static int img_direct(qnmodel* q, zp_t p, zp_t* res, void* ctx)
{
	(void)ctx;
	return solve_mod(q, p, res);
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
	 * the CRT accumulator, so their residues remain independent evidence
	 * against which a bounded reconstruction can be checked.  Drawn first
	 * from the same descending sequence, hence disjoint from the CRT
	 * primes. */
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
		/* draw a batch of distinct primes, descending */
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

		/* the first batch runs sequentially: nck() memoizes lazily and is
		 * not thread-safe until its cache is warm */
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

			/* Reconstruction is O(d^2) in the size of the modulus, so
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

/* ---- derived models ---- */
static qnmodel* qn_new(int R, int M)
{
	int r, k;
	qnmodel* q = (qnmodel*) calloc(1, sizeof(qnmodel));
	q->R = R; q->M = M;
	q->N  = (int*) calloc(R, sizeof(int));
	q->mi = (int*) calloc(M, sizeof(int));
	q->Z  = (mpz_t*) malloc(R*sizeof(mpz_t));
	for (r=0;r<R;r++) mpz_init(q->Z[r]);
	q->L  = (mpz_t**) malloc(M*sizeof(mpz_t*));
	for (k=0;k<M;k++) { q->L[k] = (mpz_t*) malloc(R*sizeof(mpz_t)); for (r=0;r<R;r++) mpz_init(q->L[k][r]); }
	return q;
}

static void qn_del(qnmodel* q)
{
	int r, k;
	if (!q) return;
	for (k=0;k<q->M;k++) { for (r=0;r<q->R;r++) mpz_clear(q->L[k][r]); free(q->L[k]); }
	for (r=0;r<q->R;r++) mpz_clear(q->Z[r]);
	free(q->L); free(q->Z); free(q->N); free(q->mi); free(q);
}

static void qn_fix_nt(qnmodel* q) { int r; q->Nt = 0; for (r=0;r<q->R;r++) q->Nt += q->N[r]; }

static mpq_t* mpq_arr(int n) { int i; mpq_t* a = (mpq_t*) malloc(n*sizeof(mpq_t)); for (i=0;i<n;i++) mpq_init(a[i]); return a; }
static void mpq_arr_free(mpq_t* a, int n) { int i; if (!a) return; for (i=0;i<n;i++) mpq_clear(a[i]); free(a); }

/* index of G^{+k}(N-1_r) in the exported vector (k, r 0-based) */
#define IDX_GK(R,k,r) (1 + (R) + (k)*(R) + (r))

/* lexicographic next permutation; 0 when the last one was reached */
static int next_perm(int* a, int n)
{
	int i = n-2, j, t;
	while (i >= 0 && a[i] >= a[i+1]) i--;
	if (i < 0) return 0;
	j = n-1; while (a[j] <= a[i]) j--;
	t = a[i]; a[i] = a[j]; a[j] = t;
	for (i++, j=n-1; i<j; i++, j--) { t = a[i]; a[i] = a[j]; a[j] = t; }
	return 1;
}

#define MAX_CLASS_ORDERS 24

/* Tiers 1-2: the model as given, then other class orders. */
static int solve_order(qnmodel* q, mpq_t* out, int* reordered)
{
	int R = q->R, M = q->M, nres = NRES(q), i, j, k, tries = 1, rc;
	*reordered = 0;
	rc = crt_solve(q, img_direct, NULL, 0, out);
	if (rc != 1 || R < 2) return rc;
	int* perm = (int*) malloc(R*sizeof(int));
	for (i=0;i<R;i++) perm[i] = i;
	while (rc == 1 && tries < MAX_CLASS_ORDERS && next_perm(perm, R))
	{
		tries++;
		qnmodel* p = qn_new(R, M);
		for (j=0;j<R;j++) {               /* position j holds class perm[j] */
			p->N[j] = q->N[perm[j]];
			mpz_set(p->Z[j], q->Z[perm[j]]);
			for (k=0;k<M;k++) mpz_set(p->L[k][j], q->L[k][perm[j]]);
		}
		for (k=0;k<M;k++) p->mi[k] = q->mi[k];
		qn_fix_nt(p);
		mpq_t* o = mpq_arr(nres);
		rc = crt_solve(p, img_direct, NULL, 0, o);
		if (rc == 0) {
			mpq_set(out[0], o[0]);
			for (j=0;j<R;j++) {
				mpq_set(out[1+perm[j]], o[1+j]);
				for (k=0;k<M;k++) mpq_set(out[IDX_GK(R,k,perm[j])], o[IDX_GK(R,k,j)]);
			}
			*reordered = 1;
		}
		mpq_arr_free(o, nres);
		qn_del(p);
	}
	free(perm);
	return rc;
}

/* Fold stations with identical demand rows.  A station of multiplicity m
 * is m identical queues, so the folded station has the summed
 * multiplicity; G and G(N-1_r) are unchanged, and an extra replica of any
 * member station is an extra replica of the folded one, so
 * G^{+k}(N-1_r) = G'^{+g(k)}(N-1_r).  Returns NULL if nothing folds. */
static qnmodel* merge_stations(qnmodel* q, int* grp)
{
	int M = q->M, R = q->R, k, h, r, nm = 0;
	int* rep = (int*) malloc(M*sizeof(int));
	for (k=0;k<M;k++)
	{
		grp[k] = -1;
		for (h=0;h<nm && grp[k]<0;h++)
		{
			for (r=0;r<R;r++) if (mpz_cmp(q->L[k][r], q->L[rep[h]][r]) != 0) break;
			if (r == R) grp[k] = h;
		}
		if (grp[k] < 0) { rep[nm] = k; grp[k] = nm++; }
	}
	if (nm == M) { free(rep); return NULL; }
	qnmodel* m = qn_new(R, nm);
	for (r=0;r<R;r++) { m->N[r] = q->N[r]; mpz_set(m->Z[r], q->Z[r]); }
	for (h=0;h<nm;h++) for (r=0;r<R;r++) mpz_set(m->L[h][r], q->L[rep[h]][r]);
	for (k=0;k<M;k++) m->mi[grp[k]] += q->mi[k];
	qn_fix_nt(m);
	free(rep);
	return m;
}

/* c with L_kr = c L_ks at every station and Z_r = c Z_s; 0 if none. */
static int class_scale(qnmodel* q, int s, int r, mpq_t c)
{
	int k, have = 0;
	mpq_t t; mpq_init(t);
	for (k=0;k<=q->M;k++)
	{
		mpz_srcptr a = (k < q->M) ? q->L[k][s] : q->Z[s];
		mpz_srcptr b = (k < q->M) ? q->L[k][r] : q->Z[r];
		if (mpz_sgn(a) == 0 && mpz_sgn(b) == 0) continue;
		if (mpz_sgn(a) == 0 || mpz_sgn(b) == 0) { have = 0; break; }
		mpz_set(mpq_numref(t), b); mpz_set(mpq_denref(t), a); mpq_canonicalize(t);
		if (!have) { mpq_set(c, t); have = 1; }
		else if (!mpq_equal(c, t)) { have = 0; break; }
	}
	mpq_clear(t);
	return have;
}

/* Aggregate classes with proportional demands and think times.  With the
 * 1/n! convention the generating function of G,
 *   exp(sum_r Z_r x_r) prod_k (1 - sum_r L_kr x_r)^(-m_k),
 * depends on such a group only through x_s + sum_i c_i x_i, so
 *   G(N) = G'(n) F(N),  F(N) = prod_groups (sum N_i)!/prod N_i! prod c_i^{N_i},
 * and likewise G(N-1_i) = G'(n-1_g) F(N-1_i) and G^{+k}(N-1_i) =
 * G'^{+k}(n-1_g) F(N-1_i), where ' is the model with each group replaced
 * by its first class at the summed population.  Returns NULL if nothing
 * aggregates; otherwise grp[i] is the merged class of i and c[i] its
 * scale factor. */
static qnmodel* merge_classes(qnmodel* q, int* grp, mpq_t* c)
{
	int R = q->R, M = q->M, s, r, k, ng = 0;
	int* rep = (int*) malloc(R*sizeof(int));
	for (r=0;r<R;r++) grp[r] = -1;
	for (s=0;s<R;s++)
	{
		if (grp[s] >= 0) continue;
		rep[ng] = s; grp[s] = ng; mpq_set_ui(c[s], 1, 1);
		for (r=s+1;r<R;r++)
			if (grp[r] < 0 && class_scale(q, s, r, c[r])) grp[r] = ng;
		ng++;
	}
	if (ng == R) { free(rep); return NULL; }
	qnmodel* m = qn_new(ng, M);
	for (s=0;s<ng;s++) {
		mpz_set(m->Z[s], q->Z[rep[s]]);
		for (k=0;k<M;k++) mpz_set(m->L[k][s], q->L[k][rep[s]]);
	}
	for (r=0;r<R;r++) m->N[grp[r]] += q->N[r];
	for (k=0;k<M;k++) m->mi[k] = q->mi[k];
	qn_fix_nt(m);
	free(rep);
	return m;
}

/* F(Nv) of merge_classes, for the group structure grp/c over R classes */
static void class_factor(mpq_t F, int* Nv, int R, int ng, int* grp, mpq_t* c)
{
	int g, i;
	mpz_t f; mpz_init(f);
	mpq_t t; mpq_init(t);
	mpq_set_ui(F, 1, 1);
	for (g=0; g<ng; g++)
	{
		unsigned long n = 0;
		for (i=0;i<R;i++) if (grp[i] == g) n += Nv[i];
		mpz_fac_ui(f, n); mpq_set_z(t, f); mpq_mul(F, F, t);
		for (i=0;i<R;i++) if (grp[i] == g)
		{
			mpz_fac_ui(f, Nv[i]); mpq_set_z(t, f); mpq_div(F, F, t);
			mpz_pow_ui(mpq_numref(t), mpq_numref(c[i]), Nv[i]);
			mpz_pow_ui(mpq_denref(t), mpq_denref(c[i]), Nv[i]);
			mpq_canonicalize(t);
			mpq_mul(F, F, t);
		}
	}
	mpz_clear(f); mpq_clear(t);
}

/* Tiers 1-3.  tier receives a short description of the route taken. */
static int solve_reduce(qnmodel* q, mpq_t* out, char* tier, size_t tsz)
{
	int R = q->R, M = q->M, rc, reordered, i, k;
	rc = solve_order(q, out, &reordered);
	if (rc == 0) { if (reordered) snprintf(tier, tsz, "class order"); return 0; }
	if (rc != 1) return rc;

	/* identical stations */
	{
		int* grp = (int*) malloc(M*sizeof(int));
		qnmodel* ms = merge_stations(q, grp);
		if (ms)
		{
			int nres2 = NRES(ms);
			mpq_t* o = mpq_arr(nres2);
			char sub[96] = "";
			rc = solve_reduce(ms, o, sub, sizeof sub);
			if (rc == 0)
			{
				mpq_set(out[0], o[0]);
				for (i=0;i<R;i++) mpq_set(out[1+i], o[1+i]);
				for (k=0;k<M;k++) for (i=0;i<R;i++) mpq_set(out[IDX_GK(R,k,i)], o[IDX_GK(R,grp[k],i)]);
				snprintf(tier, tsz, "merged stations (M=%d->%d)%s%s", M, ms->M, sub[0] ? ", " : "", sub);
			}
			mpq_arr_free(o, nres2);
			qn_del(ms);
			if (rc != 1) { free(grp); return rc; }
		}
		free(grp);
	}

	/* proportional classes */
	{
		int* grp = (int*) malloc(R*sizeof(int));
		mpq_t* c = mpq_arr(R);
		qnmodel* mc = merge_classes(q, grp, c);
		if (mc)
		{
			int ng = mc->R, nres2 = NRES(mc);
			mpq_t* o = mpq_arr(nres2);
			char sub[96] = "";
			rc = solve_reduce(mc, o, sub, sizeof sub);
			if (rc == 0)
			{
				int* Nv = (int*) malloc(R*sizeof(int));
				mpq_t F; mpq_init(F);
				for (i=0;i<R;i++) Nv[i] = q->N[i];
				class_factor(F, Nv, R, ng, grp, c);
				mpq_mul(out[0], o[0], F);
				for (i=0;i<R;i++)
				{
					if (q->N[i] == 0) {
						mpq_set_ui(out[1+i], 0, 1);
						for (k=0;k<M;k++) mpq_set_ui(out[IDX_GK(R,k,i)], 0, 1);
						continue;
					}
					Nv[i]--;
					class_factor(F, Nv, R, ng, grp, c);
					Nv[i]++;
					mpq_mul(out[1+i], o[1+grp[i]], F);
					for (k=0;k<M;k++) mpq_mul(out[IDX_GK(R,k,i)], o[IDX_GK(ng,k,grp[i])], F);
				}
				mpq_clear(F); free(Nv);
				snprintf(tier, tsz, "merged classes (R=%d->%d)%s%s", R, ng, sub[0] ? ", " : "", sub);
			}
			mpq_arr_free(o, nres2);
			qn_del(mc);
		}
		mpq_arr_free(c, R); free(grp);
	}
	return rc;
}

/* ---- tier 4: exact limit of a perturbation ----
 * Every exported constant is a homogeneous polynomial in (L, Z) jointly,
 * of degree at most Ntot, so along (L, Z) + t (DL, DZ) it is a polynomial
 * g(t) of degree <= Ntot whose value at t=0 is the constant of the
 * original model.  Both L and Z must move: with Z=0 the system can stay
 * singular for every demand matrix.  For a generic direction the system
 * is singular at finitely many t only, so each prime solves the perturbed
 * model at t=1,2,..., skips the t that fail, and evaluates the Lagrange
 * interpolant through Ntot+1 points at t=0 in Z/p.  The CRT then
 * reconstructs g(0) exactly.  The direction is drawn from 1..1e6: small
 * integers create coincidences (equal demands) that are themselves MoM
 * degeneracies. */
typedef struct { int D, T; qnmodel** qt; } interp_ctx;

static zp_t mulm(zp_t a, zp_t b, zp_t p) { return (zp_t)(((__uint128_t)a * b) % p); }
static zp_t powm(zp_t a, zp_t e, zp_t p) { zp_t r = 1; while (e) { if (e & 1) r = mulm(r, a, p); a = mulm(a, a, p); e >>= 1; } return r; }

static int img_interp(qnmodel* q, zp_t p, zp_t* res, void* vctx)
{
	interp_ctx* c = (interp_ctx*) vctx;
	int nres = NRES(q), D = c->D, got = 0, t, i, j, u;
	zp_t* buf = (zp_t*) malloc((size_t)(D+1)*nres*sizeof(zp_t));
	int* ts = (int*) malloc((D+1)*sizeof(int));
	for (t=1; t<=c->T && got<D+1; t++)
		if (solve_mod(c->qt[t-1], p, buf + (size_t)got*nres) == 0) ts[got++] = t;
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

static int solve_interp(qnmodel* q, mpq_t* out)
{
	int R = q->R, M = q->M, r, k, t, rc;
	interp_ctx c;
	unsigned long long seed = 0x9E3779B97F4A7C15ULL;
	mpz_t *dz = (mpz_t*) malloc(R*sizeof(mpz_t)), **dl = (mpz_t**) malloc(M*sizeof(mpz_t*));
#define NEXTRAND() (seed = seed*6364136223846793005ULL + 1442695040888963407ULL, 1 + (unsigned long)((seed >> 33) % 1000000))
	for (k=0;k<M;k++) { dl[k] = (mpz_t*) malloc(R*sizeof(mpz_t)); for (r=0;r<R;r++) mpz_init_set_ui(dl[k][r], NEXTRAND()); }
	for (r=0;r<R;r++) mpz_init_set_ui(dz[r], NEXTRAND());
#undef NEXTRAND
	c.D = 0; for (r=0;r<R;r++) c.D += q->N[r];
	c.T = 4*(c.D+1) + 8;
	c.qt = (qnmodel**) malloc(c.T*sizeof(qnmodel*));
	for (t=1;t<=c.T;t++)
	{
		qnmodel* m = qn_new(R, M);
		for (r=0;r<R;r++) { m->N[r] = q->N[r]; mpz_set(m->Z[r], dz[r]); mpz_mul_ui(m->Z[r], m->Z[r], t); mpz_add(m->Z[r], m->Z[r], q->Z[r]); }
		for (k=0;k<M;k++) { m->mi[k] = q->mi[k];
			for (r=0;r<R;r++) { mpz_set(m->L[k][r], dl[k][r]); mpz_mul_ui(m->L[k][r], m->L[k][r], t); mpz_add(m->L[k][r], m->L[k][r], q->L[k][r]); } }
		qn_fix_nt(m);
		c.qt[t-1] = m;
	}
	rc = crt_solve(q, img_interp, &c, (unsigned long)c.T + 1, out);
	for (t=0;t<c.T;t++) qn_del(c.qt[t]);
	free(c.qt);
	for (k=0;k<M;k++) { for (r=0;r<R;r++) mpz_clear(dl[k][r]); free(dl[k]); }
	for (r=0;r<R;r++) mpz_clear(dz[r]);
	free(dl); free(dz);
	return rc;
}

/* All tiers.  Returns 0 and the exact exported constants of q in out[],
 * or 1 if every tier failed, or 2 on a hard error already reported. */
static int solve_tiered(qnmodel* q, mpq_t* out, char* tier, size_t tsz)
{
	int rc;
	snprintf(tier, tsz, "direct");
	rc = solve_reduce(q, out, tier, tsz);
	if (rc != 1) return rc;
	rc = solve_interp(q, out);
	if (rc == 0) {
		int r, D = 0;
		for (r=0;r<q->R;r++) D += q->N[r];
		snprintf(tier, tsz, "interpolation (%d+1 perturbed solves per prime)", D);
	}
	return rc;
}

static void usage(const char* prog)
{
	printf("USAGE: %s [-l|--log] [-e|--ex] [-g|--nc] [-t|--tput] [-q|--qlen]\n", prog);
	printf("       %*s [-v] [-j nthreads] [-b bits] [-B] [-W n] [--denbits n] model.qn\n",
	       (int)strlen(prog), "");
	printf("  -l, --log        : Print only log of normalizing constant as double\n");
	printf("  -e, --ex         : Print exact normalizing constant numerator and denominator\n");
	printf("  -g, --nc         : Print normalizing constant as double\n");
	printf("  -t, --tput       : Print only throughputs, one per row\n");
	printf("  -q, --qlen       : Print only queue lengths, one per row\n");
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
	printf("                     when -B is given, 0 otherwise). Each witness reduces\n");
	printf("                     the chance of accepting a wrong candidate to about 2^-61.\n");
	printf("      --denbits n  : Denominator bound for -B, as a power of two (default 64).\n");
}

int main(int argc, char**argv)
{
	int i, r, k, nres, nthreads = 1, primebits = 61;
	int nwitness = 0, denbits = 64;
	bool log_output=false, normconst_output=false, normconst_g_output=false;
	bool throughput_output=false, queue_output=false;
	char* model_file=NULL;

	t0=CPUTIME;
	DEBUG=false; VERBOSE=true;

	if (argc < 2) { usage(argv[0]); return -1; }
	for (i=1;i<argc;i++)
	{
		if      (!strcmp(argv[i],"-l")||!strcmp(argv[i],"--log"))  log_output=true;
		else if (!strcmp(argv[i],"-e")||!strcmp(argv[i],"--ex"))   normconst_output=true;
		else if (!strcmp(argv[i],"-g")||!strcmp(argv[i],"--nc"))   normconst_g_output=true;
		else if (!strcmp(argv[i],"-t")||!strcmp(argv[i],"--tput")) throughput_output=true;
		else if (!strcmp(argv[i],"-q")||!strcmp(argv[i],"--qlen")) queue_output=true;
		else if (!strcmp(argv[i],"-v"))                            VERBOSE_MOD=1;
		else if (!strcmp(argv[i],"-h")||!strcmp(argv[i],"--help")) { usage(argv[0]); return 0; }
		else if (!strcmp(argv[i],"-j")) { if(i+1<argc) nthreads=atoi(argv[++i]); }
		else if (!strcmp(argv[i],"-b")) { if(i+1<argc) primebits=atoi(argv[++i]); }
		else if (!strcmp(argv[i],"-B")||!strcmp(argv[i],"--bounded")) { USE_BOUNDED=1; if(nwitness==0) nwitness=1; }
		else if (!strcmp(argv[i],"-W")||!strcmp(argv[i],"--witness")) { if(i+1<argc) nwitness=atoi(argv[++i]); }
		else if (!strcmp(argv[i],"--denbits")) { if(i+1<argc) denbits=atoi(argv[++i]); }
		else if (argv[i][0]=='-') { /* ignore */ }
		else model_file=argv[i];
	}
	if (!model_file) { usage(argv[0]); return -1; }
	if (primebits < 8 || primebits > 62) { fprintf(stderr,"prime bit size must be in [8,62]\n"); return -1; }

	qnm = (qnmodel*)readmodel(model_file);
	qn_require_closed(qnm, "mommod", 0);
	nckinit(qnm->M+qnm->R-1,qnm->R);
	nres = NRES(qnm);

	/* Denominator bound for reconstruction.  The recursion divides only by
	 * class populations and the pivots of the diagonal blocks, so 2^256 is
	 * generous by a wide margin; anything beyond it falls back to
	 * unbounded reconstruction. */
	mpz_init(DENBOUND);
	mpz_ui_pow_ui(DENBOUND, 2, (unsigned long) denbits);
	if (nwitness > 4) nwitness = 4;
	if (!USE_BOUNDED) nwitness = 0;

	S_NTHREADS = nthreads; S_PRIMEBITS = primebits; S_NWITNESS = nwitness;

	mpz_t *num = (mpz_t*) malloc(nres*sizeof(mpz_t));
	mpz_t *den = (mpz_t*) malloc(nres*sizeof(mpz_t));
	mpq_t *exact = mpq_arr(nres);
	char tier[256];
	double tstart = CPUTIME;
	int rc = solve_tiered(qnm, exact, tier, sizeof tier);
	if (rc == 1)
	{
		fprintf(stderr,"Error: the MoM system is singular for this model, and every exact fallback\n"
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
	mpq_arr_free(exact, nres);

	double telapsed = CPUTIME - tstart;
	if (VERBOSE_MOD)
		fprintf(stderr,"mommod: %ld primes of %d bits (%ld unlucky discarded), modulus %lu bits, %.6f s, %d thread(s)%s; route: %s\n",
		        S_TOTPRIMES, primebits, S_TOTUNLUCKY, S_MODBITS, telapsed, nthreads,
		        USE_BOUNDED ? " [bounded reconstruction, witness-verified]" : "", tier);
	else if (strcmp(tier, "direct") != 0)
		fprintf(stderr,"mommod: singular MoM system, solved exactly via %s\n", tier);

	/* ---- report, in the same formats as bin/mom ---- */
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
		printmodel(qnm);
		printf("========== Performance Metrics ==========\n");
		printf("G = %.15e\n", to_double(num[0],den[0]));
		printf("log(G) = %.15e\n", log_mpq(num[0],den[0]));
		printf("\nX (throughputs):\n");
		for (r=1;r<=qnm->R;r++) printf("X[%d] = %.15e\n", r, ratio_double(r,-1,num,den));
		printf("\nQ (mean queue lengths):\n");
		for (k=1;k<=qnm->M;k++)
		{
			double tot=0.0;
			printf("Q[%d] =", k);
			for (r=1;r<=qnm->R;r++)
			{
				double v = ratio_double(1+qnm->R+(k-1)*qnm->R+(r-1),k,num,den);
				tot += v;
				printf("\t%.15e", v);
			}
			printf("\t(total: %.15e)\n", tot);
		}
		printf("=========================================\n");
		printf("Primes: %ld of %d bits (%ld unlucky), modulus %lu bits, %d thread(s)\n",
		       S_TOTPRIMES, primebits, S_TOTUNLUCKY, S_MODBITS, nthreads);
		if (strcmp(tier, "direct") != 0) printf("Singular MoM system, solved exactly via %s\n", tier);
		t1 = CPUTIME;
		printf("Elapsed time (MoM-mod): %.6f s\n", t1-t0);
	}

	return 0;
}
