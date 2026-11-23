/*
 * The multi-modular MoM engine: one image of the recursion per prime
 * (solve_mod), the CRT and rational-reconstruction loop over any per-prime
 * image (crt_solve), and the exact fallback for singular models
 * (solve_tiered, tiers 1-3 from util/singular.h, tier 4 by interpolation
 * per prime).  bin/mommod drives it from main.c; bin/mom links it, with
 * every symbol but the mmengine.h API made local, to run its own tier 4
 * at modular speed.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include <zpla.h>
#include "mommod.h"
#include "lvtab.h"
#include "singular.h"
#include "mmengine.h"

/* Number of consecutive primes that must reproduce the same reconstructed
 * rational before it is accepted.  One extra prime is the textbook
 * criterion; two makes a spurious agreement negligibly unlikely. */
#define STABLE_ROUNDS 1
#ifdef _OPENMP
#include <omp.h>
#endif

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
				        " can supply; rerun with a larger --primebits\n"
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

/* ---- singular models: see util/singular.h ---- */
static int sg_direct(qnmodel* q, mpq_t* out, void* ctx)
{
	(void)ctx;
	return crt_solve(q, img_direct, NULL, 0, out);
}

/* Tier 4 per prime: each prime solves the perturbed models at t=1,2,...,
 * skips the t that fail, and evaluates the Lagrange interpolant at t=0 in
 * Z/p; the CRT then reconstructs g(0) exactly.  This avoids the large
 * rationals of the perturbed models that interpolation over Q would
 * reconstruct one by one. */
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

static int sg_interp_mod(qnmodel* q, mpq_t* out, void* vctx)
{
	interp_ctx c;
	int rc;
	(void)vctx;
	c.D = sg_interp_degree(q);
	c.T = sg_interp_points(q);
	c.qt = sg_perturb_build(q, c.T);
	rc = crt_solve(q, img_interp, &c, (unsigned long)c.T + 1, out);
	sg_perturb_free(c.qt, c.T);
	return rc;
}

static int solve_tiered(qnmodel* q, mpq_t* out, char* tier, size_t tsz)
{
	sg_ops ops = { sg_direct, sg_interp_mod, NULL, " per prime" };
	return sg_solve_tiered(q, out, tier, tsz, &ops, 0);
}

/* ---- API (mmengine.h) ---- */
void mommod_config(int nthreads, int primebits, int nwitness, int use_bounded, int denbits)
{
	static int init = 0;
	if (!init) { mpz_init(DENBOUND); init = 1; }
	mpz_ui_pow_ui(DENBOUND, 2, (unsigned long) denbits);
	if (nwitness > 4) nwitness = 4;
	USE_BOUNDED = use_bounded;
	if (!use_bounded) nwitness = 0;
	S_NTHREADS = nthreads > 0 ? nthreads : 1;
	S_PRIMEBITS = primebits; S_NWITNESS = nwitness;
}

void mommod_stats(long* primes, long* unlucky, unsigned long* modbits)
{
	if (primes) *primes = S_TOTPRIMES;
	if (unlucky) *unlucky = S_TOTUNLUCKY;
	if (modbits) *modbits = S_MODBITS;
}

int mommod_solve_tiered(qnmodel* q, mpq_t* out, char* tier, size_t tsz)
{
	return solve_tiered(q, out, tier, tsz);
}

int mommod_interp(qnmodel* q, mpq_t* out, void* ctx)
{
	return sg_interp_mod(q, out, ctx);
}
