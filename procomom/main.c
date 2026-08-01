#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include <gmp.h>
#include <math.h>
#include <time.h>
#include <gmpla.h>
#include "procomom.h"

/* Global variables */
qnmodel* qnm;
double t0, t1;
struct rusage ruse;

/* --progress: per-step timing on stderr, so stdout stays machine-readable.
 * procomom can run for hours on a perturbed degenerate model and the ordinary
 * -q/-P modes are silent, which makes a slow run indistinguishable from a
 * hung one.  This reports where the time actually goes: matrix assembly, the
 * A^T A products, the LU factorisation, and the n-sweep of back-substitutions.
 * PROGRESS_STATION is the rotation currently being solved. */
static int PROGRESS = 0;
static int PROGRESS_STATION = 0, PROGRESS_NSTATIONS = 0;

/* Largest entry, in bits, over the whole pk table.  On a perturbed model the
 * rationals inflate step by step and that growth is what makes the LU and the
 * back-substitutions expensive; a single n-slice is a poor sample because the
 * high-n slices are still zero early in the sweep. */
static double peak_bits(mpq_vec_t* pk, int nmax, int len)
{
	size_t best = 0; int n, i;
	for (n = 0; n <= nmax; n++)
		for (i = 0; i < len; i++) {
			size_t b = mpz_sizeinbase(mpq_numref(pk[n][i]), 2)
			         + mpz_sizeinbase(mpq_denref(pk[n][i]), 2);
			if (b > best) best = b;
		}
	return (double) best;
}

static void printcompact(int* n, int R, double elapsed_time)
{
	int s;
	static int first_call = 1;
	if (!first_call) {
		fprintf(stdout, "\r\033[K");
	} else {
		first_call = 0;
	}
	fprintf(stdout, "n=(");
	for (s = 0; s < R-1; s++)
		fprintf(stdout, "%d,", n[s]);
	fprintf(stdout, "%d) - %.6f s", n[R-1], elapsed_time);
	fflush(stdout);
}

/* Free a PMatrices struct and all its dense matrices */
static void free_pmatrices(PMatrices* pm)
{
	int i, j;
	for (i = 0; i <= pm->numRows; i++) {
		for (j = 0; j <= pm->basisSize; j++) {
			mpq_clear(pm->A[i][j]);
			mpq_clear(pm->B[i][j]);
			mpq_clear(pm->DC[i][j]);
			mpq_clear(pm->DD[i][j]);
		}
		free(pm->A[i]);
		free(pm->B[i]);
		free(pm->DC[i]);
		free(pm->DD[i]);
	}
	free(pm->A);
	free(pm->B);
	free(pm->DC);
	free(pm->DD);
	free(pm);
}

/* Free a dense square matrix allocated by mpq_matzeros */
static void free_mpq_mat(mpq_mat_t mat, int rows, int cols)
{
	int i, j;
	for (i = 0; i <= rows; i++) {
		for (j = 0; j <= cols; j++)
			mpq_clear(mat[i][j]);
		free(mat[i]);
	}
	free(mat);
}

/* Components per shift for the current rotation, and the basis size that
 * follows.  The basis is the base component plus stations 1..M-1: with a
 * single-server reference station (the one rotated to position M) its own
 * replica term has weight mi_M - 1 = 0 everywhere, so giving it a column
 * would leave that column unconstrained and A^T A singular.  When the
 * reference IS replicated the weight is nonzero and the column is needed, so
 * the stride widens to M+1. */
static int set_stride(combsrep* Dn, qnmodel* qn)
{
	Dn->stride = (qn->mi[qn->M - 1] > 1) ? qn->M + 1 : qn->M;
	return Dn->card * Dn->stride;
}

/* Swap L rows and mi values for stations a and b (0-based) */
static void swap_stations(qnmodel* qnm, int a, int b)
{
	if (a == b) return;
	mpz_t* tmp_L = qnm->L[a];
	qnm->L[a] = qnm->L[b];
	qnm->L[b] = tmp_L;
	int tmp_mi = qnm->mi[a];
	qnm->mi[a] = qnm->mi[b];
	qnm->mi[b] = tmp_mi;
}

/*
 * Run the ProCoMoM algorithm and extract P(n_M = j) for j=0..sumN.
 * Station M (last in L matrix) is the reference station.
 *
 * result_dist[j] is set to the unnormalized probability P(ref has j customers)
 * Returns 0 on success, -1 on singular matrix.
 */
static int solve_procomom(qnmodel* qnm, combsrep* Dn, int sumN,
                          mpq_t* result_dist, int basisSize, int show_progress)
{
	int M = qnm->M;
	int R = qnm->R;
	int i, r, nr, n;

	/* Allocate pk[0..sumN] and pklast[0..sumN] */
	mpq_vec_t* pk = (mpq_vec_t*)calloc(sumN + 1, sizeof(mpq_vec_t));
	mpq_vec_t* pklast = (mpq_vec_t*)calloc(sumN + 1, sizeof(mpq_vec_t));
	for (i = 0; i <= sumN; i++) {
		pk[i] = mpq_vec(basisSize, 0, 1);
		pklast[i] = mpq_vec(basisSize, 0, 1);
	}

	/* Current population vector */
	int* Ncur = (int*)int_vec(R, 0);

	/* Initialize: pexact for empty network, n=0 */
	pexact(pk[0], Dn, M);

	/* Temporary vectors for solving */
	mpq_vec_t rhs = mpq_vec(basisSize, 0, 1);
	mpq_vec_t tmp1 = mpq_vec(basisSize, 0, 1);
	mpq_vec_t tmp2 = mpq_vec(basisSize, 0, 1);
	mpq_t n_q;
	mpq_init(n_q);
	int ret = 0;

	/* Main loop: for each class r, for each Nr */
	for (r = 1; r <= R; r++) {
		for (nr = 1; nr <= qnm->N[r-1]; nr++) {
			Ncur[r-1] = nr;
			double step_start = CPUTIME;
			double step_elapsed = step_start - t0;
			if (show_progress)
				printcompact(Ncur, R, step_elapsed);

			/* Swap pk and pklast */
			mpq_vec_t* temp_ptr = pk;
			pk = pklast;
			pklast = temp_ptr;

			/* Zero out pk */
			for (n = 0; n <= sumN; n++)
				for (i = 0; i < basisSize; i++)
					mpq_set_si(pk[n][i], 0, 1);

			/* Generate matrices A, B, DC, DD */
			double t_gen0 = CPUTIME;
			PMatrices* pm = genpmatrix(Dn, qnm, Ncur, r);
			int numRows_saved = pm->numRows;
			double t_gen = CPUTIME - t_gen0;

			/* Compute normal equation matrices: basisSize x basisSize */
			mpq_mat_t AtA = mpq_matzeros(basisSize, basisSize);
			mpq_mat_t AtB = mpq_matzeros(basisSize, basisSize);
			mpq_mat_t AtDC = mpq_matzeros(basisSize, basisSize);
			mpq_mat_t AtDD = mpq_matzeros(basisSize, basisSize);

			double t_mm0 = CPUTIME;
			mpq_mattransmul(AtA, pm->A, pm->A, pm->numRows, basisSize);
			mpq_mattransmul(AtB, pm->A, pm->B, pm->numRows, basisSize);
			mpq_mattransmul(AtDC, pm->A, pm->DC, pm->numRows, basisSize);
			mpq_mattransmul(AtDD, pm->A, pm->DD, pm->numRows, basisSize);
			double t_mm = CPUTIME - t_mm0;

			/* Free the rectangular matrices */
			free_pmatrices(pm);

			/* LU decompose AtA */
			double t_lu0 = CPUTIME;
			int* lu_indices = mpq_ludcmp(AtA, basisSize);
			double t_lu = CPUTIME - t_lu0;
			if (lu_indices == NULL) {
				fprintf(stderr, "\nSingular AtA at class %d, Nr=%d (numRows=%d, basisSize=%d)\n",
				        r, nr, numRows_saved, basisSize);
				free_mpq_mat(AtA, basisSize, basisSize);
				free_mpq_mat(AtB, basisSize, basisSize);
				free_mpq_mat(AtDC, basisSize, basisSize);
				free_mpq_mat(AtDD, basisSize, basisSize);
				ret = -1;
				goto cleanup;
			}

			int sumNcur = sum(Ncur, R);
			double t_solve0 = CPUTIME;

			/* n=0: AtA * pk[0] = AtB * pklast[0] */
			mpq_matvecmul(rhs, AtB, pklast[0], basisSize);
			mpq_vecdup(pk[0], rhs, basisSize);
			if (mpq_lubksb(AtA, pk[0], basisSize, lu_indices) < 0) {
				ret = -1;
				free_mpq_mat(AtA, basisSize, basisSize);
				free_mpq_mat(AtB, basisSize, basisSize);
				free_mpq_mat(AtDC, basisSize, basisSize);
				free_mpq_mat(AtDD, basisSize, basisSize);
				free(lu_indices);
				goto cleanup;
			}

			/* n=1..sum(Ncur) */
			for (n = 1; n <= sumNcur; n++) {
				mpq_matvecmul(rhs, AtB, pklast[n], basisSize);
				mpq_matvecmul(tmp1, AtDC, pk[n-1], basisSize);
				mpq_matvecmul(tmp2, AtDD, pklast[n-1], basisSize);

				mpq_set_si(n_q, n, 1);
				for (i = 0; i < basisSize; i++) {
					mpq_t term;
					mpq_init(term);
					if (mpq_sgn(tmp1[i]) != 0) {
						mpq_mul(term, n_q, tmp1[i]);
						mpq_add(rhs[i], rhs[i], term);
					}
					if (mpq_sgn(tmp2[i]) != 0) {
						mpq_mul(term, n_q, tmp2[i]);
						mpq_add(rhs[i], rhs[i], term);
					}
					mpq_clear(term);
				}

				mpq_vecdup(pk[n], rhs, basisSize);
				if (mpq_lubksb(AtA, pk[n], basisSize, lu_indices) < 0) {
					ret = -1;
					free_mpq_mat(AtA, basisSize, basisSize);
					free_mpq_mat(AtB, basisSize, basisSize);
					free_mpq_mat(AtDC, basisSize, basisSize);
					free_mpq_mat(AtDD, basisSize, basisSize);
					free(lu_indices);
					goto cleanup;
				}
			}

			free_mpq_mat(AtA, basisSize, basisSize);
			free_mpq_mat(AtB, basisSize, basisSize);
			free_mpq_mat(AtDC, basisSize, basisSize);
			free_mpq_mat(AtDD, basisSize, basisSize);
			free(lu_indices);

			if (show_progress) {
				step_elapsed = CPUTIME - step_start;
				printf(" [%.6f s]", step_elapsed);
				fflush(stdout);
			}
			if (PROGRESS) {
				double t_solve = CPUTIME - t_solve0;
				fprintf(stderr,
				  "[procomom] station %d/%d  class %d/%d  Nr %d/%d  n<=%d  basis %d"
				  "  gen %.2fs  AtA %.2fs  LU %.2fs  solve %.2fs  step %.2fs"
				  "  total %.1fs  peak %.0f bits\n",
				  PROGRESS_STATION, PROGRESS_NSTATIONS, r, R, nr, qnm->N[r-1],
				  sumNcur, basisSize, t_gen, t_mm, t_lu, t_solve,
				  CPUTIME - step_start, CPUTIME - t0,
				  peak_bits(pk, sumNcur, basisSize));
				fflush(stderr);
			}
		}
		if (show_progress)
			printf(" (Class %d done)\n", r);
	}

	/* Extract result: pk[j][0] gives P(station M has j customers) */
	for (n = 0; n <= sumN; n++)
		mpq_set(result_dist[n], pk[n][0]);

cleanup:
	mpq_clear(n_q);
	for (i = 0; i < basisSize; i++) {
		mpq_clear(rhs[i]);
		mpq_clear(tmp1[i]);
		mpq_clear(tmp2[i]);
	}
	free(rhs);
	free(tmp1);
	free(tmp2);
	for (i = 0; i <= sumN; i++) {
		for (int j = 0; j < basisSize; j++) {
			mpq_clear(pk[i][j]);
			mpq_clear(pklast[i][j]);
		}
		free(pk[i]);
		free(pklast[i]);
	}
	free(pk);
	free(pklast);
	free(Ncur);
	return ret;
}

int main(int argc, char** argv)
{
	t0 = CPUTIME;
	int i;

	/* Parse command line arguments */
	bool queue_output = false;
	bool prob_output = false;
	int perturbation_digit = 0;
	int perturbation_seed = 23000;
	char* model_file = NULL;
	bool auto_perturbation = false;
	bool show_progress = true;

	if (argc < 2) {
		printf("USAGE: %s [-q|--qlen] [-P|--prob] [-p digit] [-s seed] [-h|--help] model.qn\n", argv[0]);
		printf("  -q, --qlen  : Print mean queue lengths (one per station, matching CoMoM format)\n");
		printf("  -P, --prob  : Print marginal probability distributions per station\n");
		printf("  -p digit    : Apply perturbation at the specified digit\n");
		printf("  -s seed     : Set perturbation seed (default: 23000)\n");
		printf("  --progress  : Report per-step timing on stderr (diagnostic)\n");
		printf("  -h, --help  : Print this help message\n");
		return -1;
	}

	for (i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-q") == 0 || strcmp(argv[i], "--qlen") == 0) {
			queue_output = true;
		} else if (strcmp(argv[i], "-P") == 0 || strcmp(argv[i], "--prob") == 0) {
			prob_output = true;
		} else if (strcmp(argv[i], "--progress") == 0) {
			PROGRESS = 1;
		} else if (strcmp(argv[i], "-p") == 0) {
			if (i + 1 < argc) {
				perturbation_digit = atoi(argv[++i]);
				if (perturbation_digit < 1) {
					printf("Error: Perturbation digit must be at least 1\n");
					return -1;
				}
			} else {
				printf("Error: -p option requires a digit argument\n");
				return -1;
			}
		} else if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--seed") == 0) {
			if (i + 1 < argc) {
				perturbation_seed = atoi(argv[++i]);
			} else {
				printf("Error: -s option requires a seed argument\n");
				return -1;
			}
		} else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
			printf("USAGE: %s [-q|--qlen] [-P|--prob] [-p digit] [-s seed] [-h|--help] model.qn\n", argv[0]);
			printf("  -q, --qlen  : Print mean queue lengths (one per station, matching CoMoM format)\n");
			printf("  -P, --prob  : Print marginal probability distributions per station\n");
			printf("  -p digit    : Apply perturbation at the specified digit\n");
			printf("  -s seed     : Set perturbation seed (default: 23000)\n");
			printf("  -h, --help  : Print this help message\n");
			return 0;
		} else if (argv[i][0] == '-') {
			/* Skip unknown options */
		} else {
			model_file = argv[i];
		}
	}

	if (model_file == NULL) {
		printf("Error: No model file specified\n");
		return -1;
	}

	qnm = (qnmodel*)readmodel(model_file);
	if (queue_output || prob_output)
		show_progress = false;

solve_attempt:;
	/* Apply perturbation if requested (same logic as comom) */
	mpz_t scale_factor;
	mpz_init(scale_factor);
	mpz_set_ui(scale_factor, 1);

	mpz_t** original_L = NULL;
	mpz_t* original_Z = NULL;

	if (perturbation_digit > 0) {
		size_t r_size = (size_t)qnm->R;
		size_t m_size = (size_t)qnm->M;
		original_L = (mpz_t**)calloc(m_size, sizeof(mpz_t*));
		original_Z = (mpz_t*)calloc(r_size, sizeof(mpz_t));
		for (int ii = 0; ii < qnm->M; ii++) {
			original_L[ii] = (mpz_t*)calloc(r_size, sizeof(mpz_t));
			for (int j = 0; j < qnm->R; j++) {
				mpz_init(original_L[ii][j]);
				mpz_set(original_L[ii][j], qnm->L[ii][j]);
			}
		}
		for (int j = 0; j < qnm->R; j++) {
			mpz_init(original_Z[j]);
			mpz_set(original_Z[j], qnm->Z[j]);
		}

		mpz_set_ui(scale_factor, 10);
		mpz_pow_ui(scale_factor, scale_factor, perturbation_digit);

		int* perm = (int*)calloc(qnm->M + 1, sizeof(int));
		for (int j = 0; j < qnm->R; j++) {
			for (int ii = 0; ii < qnm->M + 1; ii++)
				perm[ii] = ii + 1;
			unsigned int rand_state = perturbation_seed + j * 1000;
			for (int ii = qnm->M; ii > 0; ii--) {
				rand_state = rand_state * 1103515245 + 12345;
				int kk = (rand_state / 65536) % (ii + 1);
				int temp = perm[ii];
				perm[ii] = perm[kk];
				perm[kk] = temp;
			}
			for (int ii = 0; ii < qnm->M; ii++) {
				mpz_mul(qnm->L[ii][j], qnm->L[ii][j], scale_factor);
				mpz_add_ui(qnm->L[ii][j], qnm->L[ii][j], perm[ii]);
			}
			mpz_mul(qnm->Z[j], qnm->Z[j], scale_factor);
			mpz_add_ui(qnm->Z[j], qnm->Z[j], perm[qnm->M]);
		}
		free(perm);

		if (show_progress) {
			if (auto_perturbation)
				printf("\nNote: Auto perturbation at digit %d.\n", perturbation_digit);
			else
				printf("\nWarning: Perturbation at digit %d.\n\n", perturbation_digit);
		}
	}

	if (show_progress) {
		if (perturbation_digit > 0)
			printmodel_with_perturbation(qnm, perturbation_digit, scale_factor, perturbation_seed, original_L, original_Z);
		else
			printmodel(qnm);
	}

	int M = qnm->M;
	int R = qnm->R;
	nckinit(M + R + 1, M);

	/* Fixed Dn for all classes: multichoose(R, M) with column R-1 zeroed */
	combsrep* Dn = (combsrep*)calloc(1, sizeof(combsrep));
	Dn->n = R;
	Dn->k = M;
	Dn->card = nck(M + R - 1, M);
	Dn->combs = (int**)multichoose(R, M);
	for (int d = 0; d < Dn->card; d++)
		Dn->combs[d][R-1] = 0;
	Dn->combs = (int**)sortbynnzpos(Dn->combs, Dn->card, R);
	/* Width of a shift block.  The reference station is rotated to position M,
	 * so the stride is recomputed per rotation in the loops below. */
	Dn->stride = M;

	int basisSize = Dn->card * Dn->stride;
	int sumN = sum(qnm->N, R);

	/* Allocate result distribution array */
	mpq_t* dist = (mpq_t*)calloc(sumN + 1, sizeof(mpq_t));
	for (i = 0; i <= sumN; i++)
		mpq_init(dist[i]);

	if (queue_output) {
		/*
		 * Compute Q for ALL stations by rotating L so each station
		 * takes a turn as the reference (station M).
		 */
		mpf_t fval;
		mpf_init(fval);

		/* Buffer the results: on a singular solve this function jumps back to
		 * solve_attempt and recomputes every station, so anything already
		 * printed would be emitted twice.  Print only once all M stations
		 * have succeeded. */
		double* qbuf = (double*) calloc(M, sizeof(double));

		for (int k = 1; k <= M; k++) {
			/* Swap station k-1 with station M-1 to make station k the reference */
			swap_stations(qnm, k-1, M-1);
			basisSize = set_stride(Dn, qnm);
			PROGRESS_STATION = k; PROGRESS_NSTATIONS = M;
			if (PROGRESS) fprintf(stderr, "[procomom] === station %d/%d, basis %d ===\n", k, M, basisSize);

			/* Run algorithm */
			int ret = solve_procomom(qnm, Dn, sumN, dist, basisSize, 0);
			if (ret < 0) {
				/* Swap back before handling error */
				swap_stations(qnm, k-1, M-1);
				goto singular_matrix_detected;
			}

			/* Swap back */
			swap_stations(qnm, k-1, M-1);

			/* Normalize and compute Q */
			mpq_t total_prob, Q_k, term;
			mpq_init(total_prob);
			mpq_init(Q_k);
			mpq_init(term);
			mpq_set_si(total_prob, 0, 1);
			mpq_set_si(Q_k, 0, 1);

			for (int j = 0; j <= sumN; j++)
				mpq_add(total_prob, total_prob, dist[j]);

			for (int j = 0; j <= sumN; j++) {
				if (mpq_sgn(dist[j]) != 0) {
					mpq_set_si(term, j, 1);
					mpq_mul(term, term, dist[j]);
					mpq_add(Q_k, Q_k, term);
				}
			}
			if (mpq_sgn(total_prob) != 0)
				mpq_div(Q_k, Q_k, total_prob);

			/* mi_k > 1 means mi_k replicated single-server queues, and the
			 * distribution above is the marginal at ONE copy, so the station
			 * total is mi_k times its mean.  Same convention as promom -q,
			 * and what makes this comparable with ca -q summed over classes. */
			{ mpq_t mq; mpq_init(mq); mpq_set_si(mq, qnm->mi[k-1], 1);
			  mpq_mul(Q_k, Q_k, mq); mpq_clear(mq); }
			mpf_set_q(fval, Q_k);
			qbuf[k-1] = mpf_get_d(fval);

			mpq_clear(total_prob);
			mpq_clear(Q_k);
			mpq_clear(term);
		}

		for (int k = 1; k <= M; k++) printf("%.15e\n", qbuf[k-1]);
		free(qbuf);
		mpf_clear(fval);
	} else if (prob_output) {
		/*
		 * Print raw probabilities per station (machine-readable)
		 * Run once per station with rotation.
		 */
		mpf_t fval;
		mpf_init(fval);

		/* buffered for the same reason as the -q branch above */
		double* pbuf = (double*) calloc((size_t)M * (sumN+1), sizeof(double));

		for (int k = 1; k <= M; k++) {
			swap_stations(qnm, k-1, M-1);
			basisSize = set_stride(Dn, qnm);
			PROGRESS_STATION = k; PROGRESS_NSTATIONS = M;
			if (PROGRESS) fprintf(stderr, "[procomom] === station %d/%d, basis %d ===\n", k, M, basisSize);

			int ret = solve_procomom(qnm, Dn, sumN, dist, basisSize, 0);
			if (ret < 0) {
				swap_stations(qnm, k-1, M-1);
				goto singular_matrix_detected;
			}
			swap_stations(qnm, k-1, M-1);

			/* Normalize */
			mpq_t total_prob;
			mpq_init(total_prob);
			mpq_set_si(total_prob, 0, 1);
			for (int j = 0; j <= sumN; j++)
				mpq_add(total_prob, total_prob, dist[j]);

			for (int j = 0; j <= sumN; j++) {
				mpq_t p;
				mpq_init(p);
				mpq_set(p, dist[j]);
				if (mpq_sgn(total_prob) != 0)
					mpq_div(p, p, total_prob);
				mpf_set_q(fval, p);
				pbuf[(size_t)(k-1)*(sumN+1) + j] = mpf_get_d(fval);
				mpq_clear(p);
			}
			mpq_clear(total_prob);
		}

		for (int k = 1; k <= M; k++) {
			for (int j = 0; j <= sumN; j++) {
				printf("%.15e", pbuf[(size_t)(k-1)*(sumN+1) + j]);
				if (j < sumN) printf(" ");
			}
			printf("\n");
		}
		free(pbuf);
		mpf_clear(fval);
	} else {
		/*
		 * Default: run for all stations and print distributions + Q
		 */
		printf("========== Marginal Queue-Length Probabilities ==========\n");
		mpf_t fval;
		mpf_init(fval);

		for (int k = 1; k <= M; k++) {
			swap_stations(qnm, k-1, M-1);
			basisSize = set_stride(Dn, qnm);
			PROGRESS_STATION = k; PROGRESS_NSTATIONS = M;
			if (PROGRESS) fprintf(stderr, "[procomom] === station %d/%d, basis %d ===\n", k, M, basisSize);

			int ret = solve_procomom(qnm, Dn, sumN, dist, basisSize,
			                          (k == 1) ? show_progress : 0);
			if (ret < 0) {
				swap_stations(qnm, k-1, M-1);
				goto singular_matrix_detected;
			}
			swap_stations(qnm, k-1, M-1);

			printf("\nStation %d:\n", k);

			/* Normalize */
			mpq_t total_prob, p_norm, Q_k, term;
			mpq_init(total_prob);
			mpq_init(p_norm);
			mpq_init(Q_k);
			mpq_init(term);
			mpq_set_si(total_prob, 0, 1);
			mpq_set_si(Q_k, 0, 1);

			for (int j = 0; j <= sumN; j++)
				mpq_add(total_prob, total_prob, dist[j]);

			for (int j = 0; j <= sumN; j++) {
				mpq_set(p_norm, dist[j]);
				if (mpq_sgn(total_prob) != 0)
					mpq_div(p_norm, p_norm, total_prob);

				mpf_set_q(fval, p_norm);
				double pval = mpf_get_d(fval);
				if (pval > 1e-20 || j <= 5)
					printf("  P(n_%d = %d) = %.15e\n", k, j, pval);

				if (mpq_sgn(p_norm) != 0) {
					mpq_set_si(term, j, 1);
					mpq_mul(term, term, p_norm);
					mpq_add(Q_k, Q_k, term);
				}
			}

			{ mpq_t mq; mpq_init(mq); mpq_set_si(mq, qnm->mi[k-1], 1);
			  mpq_mul(Q_k, Q_k, mq); mpq_clear(mq); }
			mpf_set_q(fval, Q_k);
			printf("  Q[%d] = %.15e\n", k, mpf_get_d(fval));

			mpq_clear(total_prob);
			mpq_clear(p_norm);
			mpq_clear(Q_k);
			mpq_clear(term);
		}

		mpf_clear(fval);
		printf("\n=========================================================\n");
		t1 = CPUTIME;
		printf("Elapsed time (ProCoMoM): %g s\n", t1 - t0);
	}

	/* Cleanup */
	for (i = 0; i <= sumN; i++)
		mpq_clear(dist[i]);
	free(dist);
	mpz_clear(scale_factor);

	return 0;

singular_matrix_detected:
	if (!auto_perturbation && perturbation_digit == 0) {
		fprintf(stderr, "\nWarning: Singular matrix. Applying perturbation at digit 20.\n");
		perturbation_digit = 20;
		auto_perturbation = true;

		if (original_L == NULL) {
			size_t r_size = (size_t)qnm->R;
			size_t m_size = (size_t)qnm->M;
			original_L = (mpz_t**)calloc(m_size, sizeof(mpz_t*));
			original_Z = (mpz_t*)calloc(r_size, sizeof(mpz_t));
			for (int ii = 0; ii < qnm->M; ii++) {
				original_L[ii] = (mpz_t*)calloc(r_size, sizeof(mpz_t));
				for (int j = 0; j < qnm->R; j++) {
					mpz_init(original_L[ii][j]);
					mpz_set(original_L[ii][j], qnm->L[ii][j]);
				}
			}
			for (int j = 0; j < qnm->R; j++) {
				mpz_init(original_Z[j]);
				mpz_set(original_Z[j], qnm->Z[j]);
			}
		}

		qnm = (qnmodel*)readmodel(model_file);

		for (i = 0; i <= sumN; i++)
			mpq_clear(dist[i]);
		free(dist);
		mpz_clear(scale_factor);

		goto solve_attempt;
	} else {
		fprintf(stderr, "\nError: Singular matrix with perturbation at digit %d.\n", perturbation_digit);
		for (i = 0; i <= sumN; i++)
			mpq_clear(dist[i]);
		free(dist);
		mpz_clear(scale_factor);
		return 1;
	}
}
