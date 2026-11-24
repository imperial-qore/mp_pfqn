#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <gmp.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "util.h"
#include "gld.h"
#include "replicas.h"

#define CPUTIME (getrusage(RUSAGE_SELF,&ruse), ruse.ru_utime.tv_sec + ruse.ru_stime.tv_sec + 1e-6 * (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec))

int main(int argc, char **argv)
{
	struct rusage ruse;
	double t0, t1;
	t0 = CPUTIME;

	bool log_output = false;
	bool normconst_output = false;
	bool normconst_g_output = false;
	bool tput_output = false, qlen_output = false;
	char *model_file = NULL;

	if (argc < 2) {
		printf("USAGE: %s [-l] [-e] [-g] [-t] [-q] [-h] model.qn\n", argv[0]);
		printf("  -l : log of normalizing constant\n");
		printf("  -e : exact normalizing constant (num/den)\n");
		printf("  -g : normalizing constant as double\n");
		printf("  -t : throughputs X_r = G(N-1_r)/G(N), one per row\n");
		printf("  -q : mean queue lengths, one row per station\n");
		return -1;
	}

	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-l") == 0)
			log_output = true;
		else if (strcmp(argv[i], "-e") == 0)
			normconst_output = true;
		else if (strcmp(argv[i], "-g") == 0)
			normconst_g_output = true;
		else if (strcmp(argv[i], "-t") == 0)
			tput_output = true;
		else if (strcmp(argv[i], "-q") == 0)
			qlen_output = true;
		else if (strcmp(argv[i], "-h") == 0) {
			printf("USAGE: %s [-l] [-e] [-g] [-t] [-q] [-h] model.qn\n", argv[0]);
			return 0;
		} else if (strcmp(argv[i], "-p") == 0 || strcmp(argv[i], "-s") == 0) {
			/* perturbation options of mom/comom: accepted and ignored */
			if (i + 1 < argc) i++;
		} else if (argv[i][0] == '-') {
			fprintf(stderr, "%s: unsupported option '%s' (use -l, -e, -g, -t or -q)\n", argv[0], argv[i]);
			return 1;
		} else {
			model_file = argv[i];
		}
	}

	if (model_file == NULL) {
		printf("USAGE: %s [-l] [-e] [-g] model.qn\n", argv[0]);
		return -1;
	}

	qnmodel *qn = readmodel(model_file);
	qn_require_closed(qn, "gld", 1);
	if (!log_output && !normconst_output && !normconst_g_output && !tput_output && !qlen_output)
		printmodel(qn);

	/* mi copies of a station as one load-dependent station (util/replicas.h) */
	if (qn_replica_mu(qn))
		return 1;

	if (tput_output || qlen_output) {
		/* exact measures from the convolution tables:
		 *   X_r    = G(N - 1_r) / G(N)
		 *   Q_{m,r} = sum_{k <= N} k_r Y_m(k) G_{-m}(N - k) / G(N),
		 * G_{-m} the table without station m; a folded station's Q is the
		 * total over its copies, as in every other solver */
		int R = qn->R, M = qn->M, r, m;
		long int P = gld_boxsize(qn->N, R), i, *stride = malloc(R * sizeof(long int));
		stride[0] = 1;
		for (r = 1; r < R; r++) stride[r] = stride[r-1] * (qn->N[r-1] + 1);
		mpq_t *T = malloc(P * sizeof(mpq_t)), *Tm = malloc(P * sizeof(mpq_t)), *Y = malloc(P * sizeof(mpq_t));
		for (i = 0; i < P; i++) { mpq_init(T[i]); mpq_init(Tm[i]); mpq_init(Y[i]); }
		gld_table(T, qn->L, qn->N, qn->mu, qn->Z, M, R, qn->Nt, -1);
		mpq_t v, acc;
		mpq_init(v); mpq_init(acc);
		if (mpq_sgn(T[P-1]) == 0) { fprintf(stderr, "gld: G(N) = 0\n"); return 1; }
		if (tput_output) {
			for (r = 0; r < R; r++) {
				if (qn->N[r] > 0) mpq_div(v, T[P-1-stride[r]], T[P-1]); else mpq_set_ui(v, 0, 1);
				printf("%.15e\n", mpq_get_d(v));
			}
		} else {
			int* k = calloc(R, sizeof(int));
			for (m = 0; m < M; m++) {
				gld_table(Tm, qn->L, qn->N, qn->mu, qn->Z, M, R, qn->Nt, m);
				gld_station_factor(Y, qn->L[m], qn->mu[m], qn->N, R, qn->Nt);
				for (r = 0; r < R; r++) {
					mpq_set_ui(acc, 0, 1);
					for (int s = 0; s < R; s++) k[s] = 0;
					for (i = 0; i < P; i++) {
						if (k[r] > 0 && mpq_sgn(Y[i])) {
							mpq_mul(v, Y[i], Tm[P-1-i]);   /* index of N - k is P-1-i */
							mpz_mul_ui(mpq_numref(v), mpq_numref(v), k[r]);
							mpq_canonicalize(v);
							mpq_add(acc, acc, v);
						}
						for (int s = 0; s < R; s++) { if (k[s] < qn->N[s]) { k[s]++; break; } k[s] = 0; }
					}
					mpq_div(acc, acc, T[P-1]);
					printf("%.15e", mpq_get_d(acc));
					if (r < R - 1) printf(" ");
				}
				printf("\n");
			}
			free(k);
		}
		mpq_clear(v); mpq_clear(acc);
		for (i = 0; i < P; i++) { mpq_clear(T[i]); mpq_clear(Tm[i]); mpq_clear(Y[i]); }
		free(T); free(Tm); free(Y); free(stride);
		freemodel(qn);
		return 0;
	}

	mpq_t G;
	mpq_init(G);

	gld_multi(G, qn->L, qn->N, qn->mu, qn->Z, qn->M, qn->R, qn->Nt);

	/* Compute logG (mpq_log does not overflow when G > DBL_MAX) */
	double logG = mpq_log(G);

	if (log_output) {
		printf("%.15e\n", logG);
	} else if (normconst_output) {
		mpq_canonicalize(G);
		mpz_t num, den;
		mpz_init(num);
		mpz_init(den);
		mpq_get_num(num, G);
		mpq_get_den(den, G);
		gmp_printf("%Zd\n", num);
		gmp_printf("%Zd\n", den);
		mpz_clear(num);
		mpz_clear(den);
	} else if (normconst_g_output) {
		mpf_t fval;
		mpf_init(fval);
		mpf_set_q(fval, G);
		printf("%.15e\n", mpf_get_d(fval));
		mpf_clear(fval);
	} else {
		printf("========== GLD Normalizing Constant ==========\n");
		mpf_t fval;
		mpf_init(fval);
		mpf_set_q(fval, G);
		printf("G = %.15e\n", mpf_get_d(fval));
		printf("log(G) = %.15e\n", logG);
		mpf_clear(fval);
		t1 = CPUTIME;
		printf("Elapsed time (GLD): %.6f s\n", t1 - t0);
	}

	mpq_clear(G);
	freemodel(qn);
	return 0;
}
