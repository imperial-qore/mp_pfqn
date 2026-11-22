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
#include "singular.h"
#include "mmengine.h"



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
	printf("       %*s [-v] [-j nthreads] [--primebits bits] [-B] [-W n] [--denbits n] model.qn\n",
	       (int)strlen(prog), "");
	printf("  -l, --log        : Print only log of normalizing constant as double\n");
	printf("  -e, --ex         : Print exact normalizing constant numerator and denominator\n");
	printf("  -g, --nc         : Print normalizing constant as double\n");
	printf("  -t, --tput       : Print only throughputs, one per row\n");
	printf("  -q, --qlen       : Print only queue lengths, one per row\n");
	printf("  -v               : Report prime count and timing on stderr\n");
	printf("  -j n             : Number of primes solved concurrently (OpenMP)\n");
	printf("  --primebits bits : Bit size of the primes drawn (default 61)\n");
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
	int nwitness = 0, denbits = 64, use_bounded = 0;
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
		else if (!strcmp(argv[i],"--primebits")) { if(i+1<argc) primebits=atoi(argv[++i]); }
		else if (!strcmp(argv[i],"-b")) {
			fprintf(stderr, "mommod: -b is now --primebits (-b means bounds in mom, momf and gmom)\n");
			return 1;
		}
		else if (!strcmp(argv[i],"-B")||!strcmp(argv[i],"--bounded")) { use_bounded=1; if(nwitness==0) nwitness=1; }
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
	mommod_config(nthreads, primebits, nwitness, use_bounded, denbits);

	mpz_t *num = (mpz_t*) malloc(nres*sizeof(mpz_t));
	mpz_t *den = (mpz_t*) malloc(nres*sizeof(mpz_t));
	mpq_t *exact = sg_mpq_arr(nres);
	char tier[256];
	double tstart = CPUTIME;
	int rc = mommod_solve_tiered(qnm, exact, tier, sizeof tier);
	long totprimes, totunlucky; unsigned long modbits;
	mommod_stats(&totprimes, &totunlucky, &modbits);
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
	sg_mpq_arr_free(exact, nres);

	double telapsed = CPUTIME - tstart;
	if (VERBOSE_MOD)
		fprintf(stderr,"mommod: %ld primes of %d bits (%ld unlucky discarded), modulus %lu bits, %.6f s, %d thread(s)%s; route: %s\n",
		        totprimes, primebits, totunlucky, modbits, telapsed, nthreads,
		        use_bounded ? " [bounded reconstruction, witness-verified]" : "", tier);
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
		       totprimes, primebits, totunlucky, modbits, nthreads);
		if (strcmp(tier, "direct") != 0) printf("Singular MoM system, solved exactly via %s\n", tier);
		t1 = CPUTIME;
		printf("Elapsed time (MoM-mod): %.6f s\n", t1-t0);
	}

	return 0;
}
