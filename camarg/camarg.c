#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <gmp.h>
#include <gmpla.h>
#include "util.h"

/* camarg - exact marginal queue-length distributions by convolution.
 *
 * bin/promom and bin/procomom compute the same quantity from a moment basis,
 * and both inherit the singularities of that basis: promom declines a
 * degenerate model, procomom perturbs it.  Convolution has no coefficient
 * matrix, so it never goes singular.  camarg is therefore
 *
 *   - an exact oracle for promom/procomom on EVERY closed model, including
 *     the degenerate ones where neither can be trusted, and
 *   - the exact Tier 3 that safe_promom lacks (safe_promom/README.md).
 *
 * It is the analogue of bin/ca for distributions rather than means: the same
 * exponential-in-R lattice cost, used as a reference, not as the fast path.
 *
 * Quantity.  For a single-server station k, conditioning on its contents:
 *
 *   q_k(N, n) = sum_{|j|=n} (n! / prod_r j_r!) prod_r L_{k,r}^{j_r}
 *                           G_{-k}(N - j)
 *
 * with G_{-k} the constant of the network with ONE copy of station k
 * removed, and the multinomial counting the orderings of j jobs in that
 * station's queue.  P(n_k = n) = q_k(N,n) / G(N), and sum_n q_k(N,n) = G(N)
 * is the identity used as an internal check.
 *
 * mi_k > 1 means mi_k replicated single-server queues (as bin/ca expands
 * them), so the distribution is the marginal at ONE copy and the station
 * total is mi_k times its mean -- the same convention as promom and procomom.
 *
 * Method.  Per reference station, one Buzen convolution over the population
 * lattice omitting one copy of it, then one pass over the lattice
 * accumulating the sum above.  Cost O(M * Meff * prod_r (N_r+1)).
 */

static qnmodel* qnm;

/* g[pop] for the network built from the given virtual stations, plus the
 * delay station when any Z is nonzero.  Virtual stations are given as
 * indices into qnm->L rows; a station with mi copies appears mi times. */
static mpq_t* convolve(int* virt, int nvirt, int* planesizes, int lattice)
{
	int R = qnm->R, *N = qnm->N;
	int i, r, t;
	mpq_t* g = (mpq_t*) calloc(lattice, sizeof(mpq_t));
	for (t = 0; t < lattice; t++) mpq_init(g[t]);

	int hasZ = 0;
	for (r = 0; r < R; r++) if (mpz_sgn(qnm->Z[r]) != 0) hasZ = 1;

	int* n = initpop(R);
	if (hasZ) {
		/* delay station: prod_r Z_r^{n_r} / n_r! */
		do {
			t = popindex(n, R, planesizes);
			mpq_t term; mpq_init(term); mpq_set_si(term, 1, 1);
			for (r = 0; r < R; r++) {
				mpq_t zr; mpq_init(zr); mpq_set_z(zr, qnm->Z[r]);
				int e;
				for (e = 0; e < n[r]; e++) mpq_mul(term, term, zr);
				/* 1/n_r! exactly: factorial() returns long int and overflows
				 * past 20!, and these populations go well beyond that */
				mpz_t fz; mpz_init(fz); mpz_fac_ui(fz, (unsigned long) n[r]);
				mpq_t f; mpq_init(f); mpq_set_z(f, fz); mpq_inv(f, f);
				mpq_mul(term, term, f);
				mpz_clear(fz); mpq_clear(zr); mpq_clear(f);
			}
			mpq_set(g[t], term);
			mpq_clear(term);
		} while (nextpop(n, N, R) != -1);
	} else {
		/* empty network: G(0) = 1, everything else 0 */
		for (r = 0; r < R; r++) n[r] = 0;
		mpq_set_si(g[popindex(n, R, planesizes)], 1, 1);
	}

	/* Buzen: g_m(pop) = g_{m-1}(pop) + sum_r L_mr g_m(pop - 1_r).
	 * The lattice is swept in increasing population, which nextpop does,
	 * so the g_m(pop - 1_r) on the right is already the updated value. */
	for (i = 0; i < nvirt; i++) {
		int m = virt[i];
		for (r = 0; r < R; r++) n[r] = 0;
		do {
			t = popindex(n, R, planesizes);
			mpq_t acc; mpq_init(acc); mpq_set(acc, g[t]);
			for (r = 0; r < R; r++) {
				if (n[r] == 0) continue;
				n[r]--;
				int tp = popindex(n, R, planesizes);
				n[r]++;
				mpq_t term; mpq_init(term);
				mpq_set_z(term, qnm->L[m][r]);
				mpq_mul(term, term, g[tp]);
				mpq_add(acc, acc, term);
				mpq_clear(term);
			}
			mpq_set(g[t], acc);
			mpq_clear(acc);
		} while (nextpop(n, N, R) != -1);
	}
	free(n);
	return g;
}

static void free_g(mpq_t* g, int lattice)
{
	int t;
	for (t = 0; t < lattice; t++) mpq_clear(g[t]);
	free(g);
}

int main(int argc, char** argv)
{
	int i, r, t;
	bool prob_output = false, queue_output = false;
	char* model_file = NULL;

	for (i = 1; i < argc; i++) {
		if      (!strcmp(argv[i], "-P") || !strcmp(argv[i], "--prob")) prob_output = true;
		else if (!strcmp(argv[i], "-q") || !strcmp(argv[i], "--qlen")) queue_output = true;
		else if (argv[i][0] != '-') model_file = argv[i];
	}
	if (!model_file) {
		printf("USAGE: %s [-P|-q] model.qn\n", argv[0]);
		printf("  camarg: exact marginal queue-length distributions by convolution.\n");
		printf("  Never singular, so it is the reference for bin/promom and\n");
		printf("  bin/procomom on degenerate models, where those two decline or\n");
		printf("  perturb.  Cost is the convolution lattice, exponential in R.\n");
		printf("  -P, --prob : distribution per station, one row per station\n");
		printf("  -q, --qlen : mean queue length per station (station total)\n");
		return -1;
	}

	qnm = (qnmodel*) readmodel(model_file);
	int M = qnm->M, R = qnm->R, *N = qnm->N;
	{ int c; for (c = 0; c < R; c++) if (N[c] < 0) {
		fprintf(stderr, "camarg: open/mixed models are not supported\n"); return 1; } }
	if (qnm->hasOpen) { fprintf(stderr, "camarg: open/mixed models are not supported\n"); return 1; }
	if (qnm->isLD)    { fprintf(stderr, "camarg: load-dependent stations are not supported\n"); return 1; }

	int* planesizes = getplanesizes(N, R);
	int lattice = planesizes[0];
	int sumN = sum(N, R);

	/* virtual station list: mi_k copies of each station k */
	int Meff = 0;
	for (i = 0; i < M; i++) Meff += qnm->mi[i];
	int* allvirt = (int*) calloc(Meff, sizeof(int));
	{ int v = 0, c; for (i = 0; i < M; i++) for (c = 0; c < qnm->mi[i]; c++) allvirt[v++] = i; }

	mpq_t* q = (mpq_t*) calloc(sumN+1, sizeof(mpq_t));
	for (t = 0; t <= sumN; t++) mpq_init(q[t]);

	int k;
	for (k = 0; k < M; k++) {
		/* every virtual station except ONE copy of station k */
		int* virt = (int*) calloc(Meff, sizeof(int));
		int nvirt = 0, dropped = 0;
		for (i = 0; i < Meff; i++) {
			if (!dropped && allvirt[i] == k) { dropped = 1; continue; }
			virt[nvirt++] = allvirt[i];
		}
		mpq_t* gm = convolve(virt, nvirt, planesizes, lattice);

		for (t = 0; t <= sumN; t++) mpq_set_si(q[t], 0, 1);

		/* q[n] += (n!/prod j_r!) prod L_kr^{j_r} G_-k(N-j) over the lattice */
		int* j = initpop(R);
		do {
			int n = sum(j, R);
			int* comp = (int*) calloc(R, sizeof(int));
			for (r = 0; r < R; r++) comp[r] = N[r] - j[r];
			mpq_t term; mpq_init(term);
			mpq_set(term, gm[popindex(comp, R, planesizes)]);
			free(comp);
			if (mpq_sgn(term) != 0) {
				for (r = 0; r < R; r++) {
					mpq_t lr; mpq_init(lr); mpq_set_z(lr, qnm->L[k][r]);
					int e;
					for (e = 0; e < j[r]; e++) mpq_mul(term, term, lr);
					mpq_clear(lr);
				}
				/* multinomial n!/prod j_r!, exact: n reaches sum(N), far past
				 * the 20! where a long-int factorial silently overflows */
				mpz_t num, den, fz;
				mpz_init(num); mpz_init(den); mpz_init(fz);
				mpz_fac_ui(num, (unsigned long) n);
				mpz_set_ui(den, 1);
				for (r = 0; r < R; r++) {
					mpz_fac_ui(fz, (unsigned long) j[r]);
					mpz_mul(den, den, fz);
				}
				mpz_divexact(num, num, den);
				mpq_t mult; mpq_init(mult); mpq_set_z(mult, num);
				mpq_mul(term, term, mult);
				mpq_clear(mult); mpz_clear(num); mpz_clear(den); mpz_clear(fz);
				mpq_add(q[n], q[n], term);
			}
			mpq_clear(term);
		} while (nextpop(j, N, R) != -1);
		free(j);
		free_g(gm, lattice);
		free(virt);

		/* normalise: sum_n q[n] = G(N) */
		mpq_t total; mpq_init(total); mpq_set_si(total, 0, 1);
		for (t = 0; t <= sumN; t++) mpq_add(total, total, q[t]);
		if (mpq_sgn(total) == 0) {
			fprintf(stderr, "camarg: zero normalising constant at station %d\n", k+1);
			return 1;
		}

		mpq_t Qk, p, tmp;
		mpq_init(Qk); mpq_init(p); mpq_init(tmp);
		mpq_set_si(Qk, 0, 1);
		for (t = 0; t <= sumN; t++) {
			mpq_set(p, q[t]);
			mpq_div(p, p, total);
			if (prob_output) {
				printf("%.15e%s", mpq_get_d(p), t < sumN ? " " : "");
			} else if (!queue_output) {
				double pv = mpq_get_d(p);
				if (t == 0) printf("\nStation %d:\n", k+1);
				if (pv > 1e-20 || t <= 5) printf("  P(n_%d = %d) = %.15e\n", k+1, t, pv);
			}
			mpq_set_si(tmp, t, 1);
			mpq_mul(tmp, tmp, p);
			mpq_add(Qk, Qk, tmp);
		}
		if (prob_output) printf("\n");
		/* station total = mi_k x the per-copy mean */
		mpq_set_si(tmp, qnm->mi[k], 1);
		mpq_mul(Qk, Qk, tmp);
		if (queue_output)      printf("%.15e\n", mpq_get_d(Qk));
		else if (!prob_output) printf("  Q[%d] = %.15e\n", k+1, mpq_get_d(Qk));
		mpq_clear(Qk); mpq_clear(p); mpq_clear(tmp); mpq_clear(total);
	}

	for (t = 0; t <= sumN; t++) mpq_clear(q[t]);
	free(q); free(allvirt); free(planesizes);
	return 0;
}
