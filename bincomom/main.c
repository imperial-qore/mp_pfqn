/*
 * bincomom - CoMoM with the population steps of each class multiplied
 * together by binary splitting.
 *
 * This is bin/comom (its main.c compiled here as comom_cli_main) with one
 * hook set: once the system of class r is set up, COMOM_ADVANCE may move
 * the basis from n_r = 0 to n_r = N_r - (M+1) in one go.  The last M+1
 * steps stay with comom, because they feed the next class's initial basis
 * and the final marginals, so every report, flag and the singular-model
 * fallback are comom's own.
 *
 * Within class r, comom's step (comom_step) is
 *
 *   G(n)  = B2 g(n-1) / n,     Gk(n) = A11^{-1} (B1 g(n-1) - A12 G(n))
 *
 * (A11 solved by BTF or LU, with the impossible-combination handling of
 * comom), every operation linear in g and every pivot chosen from the
 * coefficients alone.  Hence
 *
 *   n g(n) = (n P + Q) g(n-1)
 *
 * with P and Q fixed within the class.  They are read off comom_step
 * itself: column j of P + Q is 1*step(e_j, 1), of 2P + Q is
 * 2*step(e_j, 2).  Reading them through the step, instead of rebuilding
 * them from A11, A12, B1 and B2, keeps every special path of comom in
 * the matrices.  A check against one real step guards the linearity.
 * The steps are then multiplied by util/binsplit.h.
 *
 * A class is multiplied out only when N_r - (M+1) >= cross * d^2: the step
 * matrices are dense d x d (A11^{-1} fills them in) while comom steps on
 * sparse blocks, so below that the sequential steps are cheaper.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <gmp.h>
#include <gmpla.h>
#include "comom.h"
#include "binsplit.h"

int comom_cli_main(int argc, char** argv);

static double CROSS = 2.0;
static int    VERB = 0;

static void vfree(mpq_vec_t v, long n) { long i; for (i=0;i<=n;i++) mpq_clear(v[i]); free(v); }

static int bc_advance(LS* ls, combsrep* Dn, int r, long target, mpq_vec_t g, long cardG, long cardGk)
{
	int d = (int)(cardG + cardGk), i, j, rc = 1;
	if ((double) target < CROSS * (double) d * (double) d) {
		if (VERB) fprintf(stderr, "bincomom: class %d  d=%d  N_r=%ld  sequential (comom steps)\n",
		                  r, d, target + qnm->M + 1);
		return 0;
	}
	double t_set = CPUTIME;
	mpq_vec_t e   = (mpq_vec_t) mpq_vec(d, 0, 1);
	mpq_vec_t G   = (mpq_vec_t) mpq_vec(cardG, 0, 1);
	mpq_vec_t b1  = (mpq_vec_t) mpq_vec(cardGk, 0, 1);
	mpq_vec_t b1b = (mpq_vec_t) mpq_vec(cardGk, 0, 1);
	mpq_t* P = malloc((size_t)d*d*sizeof(mpq_t));
	mpq_t* Q = malloc((size_t)d*d*sizeof(mpq_t));
	mpq_t two, t;
	mpq_init(two); mpq_init(t); mpq_set_ui(two, 2, 1);
	for (i=0;i<d*d;i++) { mpq_init(P[i]); mpq_init(Q[i]); }

	/* P[:,j] = v2 - v1, Q[:,j] = 2 v1 - v2 with v1 = step(e_j,1), v2 = 2 step(e_j,2) */
	for (j=0;j<d && rc>0;j++) {
		mpq_set_ui(e[j], 1, 1);
		if (comom_step(ls, Dn, r, 1, e, G, b1, b1b, cardG, cardGk, false) < 0) { rc = -1; break; }
		for (i=0;i<cardGk;i++) mpq_set(Q[i*d+j], b1[i]);
		for (i=0;i<cardG;i++)  mpq_set(Q[(cardGk+i)*d+j], G[i]);
		if (comom_step(ls, Dn, r, 2, e, G, b1, b1b, cardG, cardGk, false) < 0) { rc = -1; break; }
		for (i=0;i<cardGk;i++) mpq_mul(P[i*d+j], b1[i], two);
		for (i=0;i<cardG;i++)  mpq_mul(P[(cardGk+i)*d+j], G[i], two);
		for (i=0;i<d;i++) {
			/* Q holds v1, P holds v2 */
			mpq_mul(t, Q[i*d+j], two);
			mpq_sub(t, t, P[i*d+j]);          /* 2 v1 - v2 */
			mpq_sub(P[i*d+j], P[i*d+j], Q[i*d+j]); /* v2 - v1 */
			mpq_set(Q[i*d+j], t);
		}
		mpq_set_ui(e[j], 0, 1);
	}

	/* linearity guard: (P + Q) g must equal comom's own first step */
	if (rc > 0) {
		if (comom_step(ls, Dn, r, 1, g, G, b1, b1b, cardG, cardGk, false) < 0) rc = -1;
		for (i=0;i<d && rc>0;i++) {
			mpq_set_ui(e[i], 0, 1);
			for (j=0;j<d;j++) {
				if (mpq_sgn(g[j]) == 0) continue;
				mpq_add(t, P[i*d+j], Q[i*d+j]);
				mpq_mul(t, t, g[j]);
				mpq_add(e[i], e[i], t);
			}
			if (!mpq_equal(e[i], i < cardGk ? b1[i] : G[i-cardGk])) {
				fprintf(stderr, "bincomom: class %d step is not linear as expected; using comom's steps\n", r);
				rc = 0;
			}
		}
	}

	if (rc > 0) {
		zmat Pz, Qz;
		mpz_t D, den;
		mpz_init(D); mpz_init(den);
		bs_integerize((const mpq_t*) P, (const mpq_t*) Q, d, &Pz, &Qz, D);
		double t_run = CPUTIME;
		mpz_t* v = malloc((size_t)d*sizeof(mpz_t));
		for (i=0;i<d;i++) mpz_init(v[i]);
		bs_vec_to_z((const mpq_t*) g, d, v, den);
		bs_steps_apply(&Pz, &Qz, 1, target, v, 1);
		bs_scale_den(den, D, 1, target);
		bs_vec_from_z(g, d, (const mpz_t*) v, den);
		if (VERB)
			fprintf(stderr, "bincomom: class %d  d=%d  N_r=%ld  binary splitting to n_r=%ld  D %lu bits  setup %.3fs  steps %.3fs\n",
			        r, d, target + qnm->M + 1, target, (unsigned long) mpz_sizeinbase(D, 2),
			        t_run - t_set, CPUTIME - t_run);
		for (i=0;i<d;i++) mpz_clear(v[i]);
		free(v);
		zm_free(Pz); zm_free(Qz);
		mpz_clear(D); mpz_clear(den);
	}

	for (i=0;i<d*d;i++) { mpq_clear(P[i]); mpq_clear(Q[i]); }
	free(P); free(Q);
	mpq_clear(two); mpq_clear(t);
	vfree(e, d); vfree(G, cardG); vfree(b1, cardGk); vfree(b1b, cardGk);
	return rc;
}

static void usage(const char* p)
{
	printf("USAGE: %s [comom options] [-v] [--cross c] [--leaf k] [--seq] model.qn\n", p);
	printf("  bincomom: CoMoM whose population steps are multiplied together by binary\n");
	printf("  splitting (exact, integer arithmetic); output identical to bin/comom.\n");
	printf("  -l -e -g -t -q -d -X -p -s : as bin/comom\n");
	printf("  -j n       : threads for the matrix products and for comom's tier 4\n");
	printf("  -v         : per-class basis size, route and timing on stderr\n");
	printf("  --cross c  : binary splitting for a class when N_r-M-1 >= c * d^2 (default 2)\n");
	printf("  --leaf k   : steps multiplied directly at a leaf (default 8)\n");
	printf("  --seq      : never split (comom's sequential steps throughout)\n");
}

int main(int argc, char** argv)
{
	char** av = malloc((size_t)(argc + 1) * sizeof(char*));
	int ac = 0, i;
	av[ac++] = argv[0];
	for (i=1;i<argc;i++) {
		if      (!strcmp(argv[i],"-v")) VERB = 1;
		else if (!strcmp(argv[i],"--cross") && i+1<argc) CROSS = atof(argv[++i]);
		else if (!strcmp(argv[i],"--leaf") && i+1<argc) { BS_LEAF = atoi(argv[++i]); if (BS_LEAF < 1) BS_LEAF = 1; }
		else if (!strcmp(argv[i],"--seq")) CROSS = 1e300;
		else if (!strcmp(argv[i],"-h") || !strcmp(argv[i],"--help")) { usage(argv[0]); free(av); return 0; }
		else {
			if (!strcmp(argv[i],"-j") && i+1<argc) BS_NTHREADS = atoi(argv[i+1]);
			av[ac++] = argv[i];
		}
	}
	av[ac] = NULL;
	COMOM_ADVANCE = bc_advance;
	int rc = comom_cli_main(ac, av);
	free(av);
	return rc;
}
