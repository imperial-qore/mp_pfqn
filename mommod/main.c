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
 *     is singular modulo every prime, and is reported after a run of
 *     consecutive failures;
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

	zp_setmod(p);

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
			gr = zp_vec(m+1,1,1);
			Gk = zp_vec(m,1,1);
			G  = zp_vec(1,1,1);
		}
		else
		{
			freels(A);
			A = setupls(qnm->L,n,qnm->Z,qnm->mi,m,r);
			if (A == NULL) { rc = -1; goto cleanup; }
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
			for(s=0;s<r-1;s++) if(mpz_cmp_ui(qnm->Z[s],0)==0)
				for(t=0;t<nck(m+r-2,r-1);t++) g[cardGk+t*r+1+s]=0;
		}

		for (n[r-1]=1;n[r-1]<=qnm->N[r-1];n[r-1]++)
		{
			if(n[r-1]==qnm->N[r-1]) zp_vecdup(gr,g,cardGk+cardG);

			zp_mspvecmul(G,A->B2r,g);
			{
				zp_t ninv = zp_invv(zp_from_si(n[r-1]));
				for(t=0;t<cardG;t++) G[t] = zp_mulv(ninv, G[t]);
			}

			zp_mspvecmul(b,A->B1r,g);
			zp_mspvecmul(b2,A->A12,G);
			for(t=0;t<cardGk;t++) b[t] = zp_subv(b[t],b2[t]);

			Gk = blocksolve(A,b);
			if (Gk == NULL) { rc = -1; goto cleanup; }

			for(t=0;t<cardGk;t++) g[t] = Gk[t];
			for(t=0;t<cardG;t++)  g[cardGk+t] = G[t];
			zp_vecfree(Gk); Gk=NULL;
			for(s=0;s<r-1;s++) if(mpz_cmp_ui(qnm->Z[s],0)==0)
				for(t=0;t<nck(m+r-2,r-1);t++) g[cardGk+t*r+1+s]=0;
		}
		n[r-1]=qnm->N[r-1];
	}

	if (ZP_SINGULAR) { rc = -1; goto cleanup; }
	mdecrease(qnm,G,Gk,g,gr,res);
	if (ZP_SINGULAR) rc = -1;
	else { int u; for (u=0;u<NRES(qnm);u++) res[u] = zp_from_mont(res[u]); }

cleanup:
	(void)h; (void)k;
	freels(A);
	zp_vecfree(b); zp_vecfree(b2); zp_vecfree(g); zp_vecfree(gr);
	zp_vecfree(G); zp_vecfree(Gk);
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
 * quantity exported here: mdecrease divides only by populations, station
 * multiplicities and think times.  Bounding the denominator recovers the
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
	unsigned long primefloor = 3;
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
	nckinit(qnm->M+qnm->R-1,qnm->R);
	nres = NRES(qnm);

	/* Denominator bound for reconstruction.  mdecrease divides only by
	 * populations, station multiplicities and think times, so 2^256 is
	 * generous by a wide margin; anything beyond it falls back to
	 * unbounded reconstruction. */
	mpz_init(DENBOUND);
	mpz_ui_pow_ui(DENBOUND, 2, (unsigned long) denbits);
	if (nwitness > 4) nwitness = 4;
	if (!USE_BOUNDED) nwitness = 0;

	/* Warm the shared nchoosek cache with a single sequential image before
	 * any concurrency: nck() memoizes lazily and is not otherwise
	 * thread-safe, but is read-only once populated. */

	mpz_t M, *acc, *num, *den, *pnum, *pden, pz;
	mpz_t *cnum, *cden;
	acc  = (mpz_t*) malloc(nres*sizeof(mpz_t));
	cnum = (mpz_t*) malloc(nres*sizeof(mpz_t));
	cden = (mpz_t*) malloc(nres*sizeof(mpz_t));
	num  = (mpz_t*) malloc(nres*sizeof(mpz_t));
	den  = (mpz_t*) malloc(nres*sizeof(mpz_t));
	pnum = (mpz_t*) malloc(nres*sizeof(mpz_t));
	pden = (mpz_t*) malloc(nres*sizeof(mpz_t));
	/* per-value stability tracking is folded into the probe on G */
	mpz_init_set_ui(M,0); mpz_init(pz);
	for (i=0;i<nres;i++) { mpz_init(acc[i]); mpz_init(num[i]); mpz_init(den[i]); mpz_init(pnum[i]); mpz_init(pden[i]); mpz_init(cnum[i]); mpz_init(cden[i]); }

	/* Primes are drawn downwards from 2^primebits so that they are
	 * distinct and comfortably larger than any N_r. */
	mpz_t cand; mpz_init(cand);
	mpz_ui_pow_ui(cand, 2, (unsigned long)primebits);
	/* Lowest usable modulus.  A prime must exceed every class population,
	 * otherwise the division by n_r in the recursion is not invertible and
	 * the image is discarded as unlucky; and it must be odd, because the
	 * Montgomery reduction inverts p modulo 2^64. */
	{
		int maxN = 0;
		for (i=0;i<qnm->R;i++) if (qnm->N[i] > maxN) maxN = qnm->N[i];
		primefloor = (unsigned long)(maxN > 2 ? maxN + 1 : 3);
		if (mpz_cmp_ui(cand, primefloor) <= 0)
		{
			fprintf(stderr,"Error: -b %d gives no prime above the largest class population %d;"
			               " raise -b\n", primebits, maxN);
			return 1;
		}
	}

	/* Witness primes: solved like any other image but never folded into
	 * the CRT accumulator, so their residues remain independent evidence
	 * against which a bounded reconstruction can be checked.  Drawn first
	 * from the same descending sequence, hence disjoint from the CRT
	 * primes. */
	{
		int w;
		for (w=0; w<nwitness; w++)
		{
			WRES[w] = (zp_t*) malloc((size_t)nres*sizeof(zp_t));
			for (;;)
			{
				if (!next_prime_below(cand, primefloor, &WPRIME[w]))
				{
					fprintf(stderr,"Error: exhausted the primes of %d bits while drawing"
					               " witnesses; raise -b\n", primebits);
					return 1;
				}
				if (solve_mod(qnm, WPRIME[w], WRES[w]) == 0) break;
			}
		}
		NWITNESS = nwitness;
	}

	int nprimes=0, unlucky=0, consec_fail=0, done=0;
	int next_try=4, gstable=0, confirmed=0;
	int batch = nthreads>1 ? nthreads : 1;
	zp_t* plist = (zp_t*) malloc(batch*sizeof(zp_t));
	zp_t* rbuf  = (zp_t*) malloc((size_t)batch*nres*sizeof(zp_t));
	int*  ok    = (int*)  malloc(batch*sizeof(int));

	double tstart = CPUTIME;

	while (!done)
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
				        primebits, nprimes, (unsigned long)mpz_sizeinbase(M,2));
				return 1;
			}

#ifdef _OPENMP
		#pragma omp parallel for num_threads(nthreads) schedule(dynamic,1) if(batch>1 && nprimes>0)
#endif
		for (bi=0; bi<batch; bi++)
			ok[bi] = (solve_mod(qnm, plist[bi], rbuf + (size_t)bi*nres) == 0);

		for (bi=0; bi<batch && !done; bi++)
		{
			if (!ok[bi]) { unlucky++; consec_fail++;
				if (consec_fail > 40) {
					fprintf(stderr,"Error: the linear system is singular modulo %d consecutive primes;\n"
					               "the coefficient matrix is singular over the rationals.\n", consec_fail);
					return 1;
				}
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

	double telapsed = CPUTIME - tstart;
	if (VERBOSE_MOD)
		fprintf(stderr,"mommod: %d primes of %d bits (%d unlucky discarded), modulus %lu bits, %.6f s, %d thread(s)%s\n",
		        nprimes, primebits, unlucky, (unsigned long)mpz_sizeinbase(M,2), telapsed, nthreads,
		        USE_BOUNDED ? " [bounded reconstruction, witness-verified]" : "");

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
		printf("Primes: %d of %d bits (%d unlucky), modulus %lu bits, %d thread(s)\n",
		       nprimes, primebits, unlucky, (unsigned long)mpz_sizeinbase(M,2), nthreads);
		t1 = CPUTIME;
		printf("Elapsed time (MoM-mod): %.6f s\n", t1-t0);
	}

	return 0;
}
