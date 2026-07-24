#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include <gmpla.h>
#include "gmom.h"

struct rusage ruse;
double t0, t1;

extern void gmva(mpq_t G, mpz_t** dem, int nq, int* pop, mpz_t* Z, int r);
extern void gmom_measures(qnmodel* qn, mpq_vec_t vlM, mpq_vec_t grM, int oe, int og, int ol, int ot, int oq, int* cperm, mpz_t scale_factor);

/* ---- small dense helpers over rectangular mpq matrices ---------------- */

/* rop[nrows] = A[nrows x ncols] * v[ncols] */
static void matvec_rect(mpq_vec_t rop, mpq_mat_t A, int nrows, int ncols, mpq_vec_t v)
{
	int i, j; mpq_t t; mpq_init(t);
	for (i = 0; i < nrows; i++) {
		mpq_set_si(rop[i], 0, 1);
		for (j = 0; j < ncols; j++)
			if (mpq_sgn(A[i][j]) && mpq_sgn(v[j])) { mpq_mul(t, A[i][j], v[j]); mpq_add(rop[i], rop[i], t); }
	}
	mpq_clear(t);
}

/* rop[ncols] = A^T[ncols x nrows] * v[nrows] */
static void tvec_rect(mpq_vec_t rop, mpq_mat_t A, int nrows, int ncols, mpq_vec_t v)
{
	int i, j; mpq_t t; mpq_init(t);
	for (j = 0; j < ncols; j++) {
		mpq_set_si(rop[j], 0, 1);
		for (i = 0; i < nrows; i++)
			if (mpq_sgn(A[i][j]) && mpq_sgn(v[i])) { mpq_mul(t, A[i][j], v[i]); mpq_add(rop[j], rop[j], t); }
	}
	mpq_clear(t);
}

/* Solve the overdetermined A x = rhs (A: nrows x ncols) by normal
 * equations (A^T A) x = A^T rhs, exactly.  Returns a fresh x[ncols]. */
static mpq_vec_t solve_normal(mpq_mat_t A, int nrows, int ncols, mpq_vec_t rhs)
{
	mpq_mat_t AtA = mpq_matzeros(ncols, ncols);
	mpq_vec_t Atb = mpq_vec(ncols, 0, 1);
	mpq_mattransmul(AtA, A, A, nrows, ncols);
	tvec_rect(Atb, A, nrows, ncols, rhs);
	int* indx = mpq_ludcmp(AtA, ncols);
	if (!indx) { fprintf(stderr, "gmom: SINGULAR A^T A  nrows=%d ncols=%d\n", nrows, ncols);
		if (getenv("GMOM_DBG")) { fprintf(stderr,"A =\n"); mpq_matprint(A, nrows, ncols); }
		exit(2); }
	mpq_lubksb(AtA, Atb, ncols, indx);
	free(indx);
	int i, j;
	for (i = 0; i < ncols+1; i++) { for (j = 0; j < ncols+1; j++) mpq_clear(AtA[i][j]); free(AtA[i]); }
	free(AtA);
	return Atb; /* length ncols; caller frees */
}

/* G of the prefix sub-model on queues 0..m-1 with (1+combo[j]) copies of
 * queue j (Ladd), at population pop, think times Z, r classes. */
static void base_G(mpq_t G, mpz_t** L, int m, int* combo, int r, int* pop, mpz_t* Z)
{
	int j, c, nq = 0;
	for (j = 0; j < m; j++) nq += 1 + combo[j];
	mpz_t** dem = (mpz_t**) malloc(nq * sizeof(mpz_t*));
	c = 0;
	for (j = 0; j < m; j++) { int t; for (t = 0; t <= combo[j]; t++) dem[c++] = L[j]; }
	gmva(G, dem, nq, pop, Z, r);
	free(dem);
}

static void oner(int* pop, int* nvec, int r, int s)
{
	int c; for (c = 0; c < r; c++) pop[c] = nvec[c];
	if (s >= 1) pop[s-1] -= 1;
}

/* free/alloc an mpq vector slot */
static mpq_vec_t vnew(int n) { return mpq_vec(n, 0, 1); }
static void vfree(mpq_vec_t v, int n) { int i; if (!v) return; for (i = 0; i < n+1; i++) mpq_clear(v[i]); free(v); }
static void fp_copy_vec(mpq_vec_t d, mpq_vec_t s, int n) { int i; for (i = 0; i < n; i++) mpq_set(d[i], s[i]); }


/* Fill V[m] both slots directly from exact convolution at the current
 * nvec, for any level m.  Used to initialise a new class exactly (the
 * cheap shift-based init is only valid when intermediate classes have no
 * think time). */
static void fill_level_all(mpq_vec_t* vl, mpq_vec_t* vl1, int* lenl, int* lenl1,
                           mpz_t** L, int m, int r, int* nvec, mpz_t* Z)
{
	int i, s;
	int* pop = (int*) malloc(r * sizeof(int));
	int szIk = nck(m + (r-1) - 1, m-1);
	int szI  = nck(m + (r-2) - 1, m-1);
	int** Ik = sortbynnzpos(multichoose(m, r-1), szIk, m);
	int** I  = sortbynnzpos(multichoose(m, r-2), szI, m);
	vfree(vl[m], lenl[m]); lenl[m] = szIk*r; vl[m] = vnew(szIk*r);
	for (i = 0; i < szIk; i++) for (s = 0; s < r; s++) { oner(pop, nvec, r, s); base_G(vl[m][i*r+s], L, m, Ik[i], r, pop, Z); }
	vfree(vl1[m], lenl1[m]); lenl1[m] = szI*r; vl1[m] = vnew(szI*r);
	for (i = 0; i < szI; i++) for (s = 0; s < r; s++) { oner(pop, nvec, r, s); base_G(vl1[m][i*r+s], L, m, I[i], r, pop, Z); }
	free(Ik); free(I); free(pop);
}
/* Fill V[m=1] base levels (both slots) at the current nvec for class r. */
static void fill_level1(mpq_vec_t* vl, mpq_vec_t* vl1, int* lenl, int* lenl1,
                        mpz_t** L, int r, int* nvec, mpz_t* Z)
{
	int s;
	int* combo0 = (int*) calloc(1, sizeof(int)); /* multichoose(1,k) = [k] */
	int* pop = (int*) malloc(r * sizeof(int));

	/* slot l_1: multichoose(1, r-2) = [r-2] (empty combo if r-2<=0 -> [max(0,r-2)]) */
	int k1 = (r-2 < 0) ? 0 : r-2;
	vfree(vl1[1], lenl1[1]); lenl1[1] = r; vl1[1] = vnew(r);
	combo0[0] = k1;
	for (s = 0; s < r; s++) { oner(pop, nvec, r, s); base_G(vl1[1][s], L, 1, combo0, r, pop, Z); }

	/* slot l: multichoose(1, r-1) = [r-1] */
	vfree(vl[1], lenl[1]); lenl[1] = r; vl[1] = vnew(r);
	combo0[0] = r-1;
	for (s = 0; s < r; s++) { oner(pop, nvec, r, s); base_G(vl[1][s], L, 1, combo0, r, pop, Z); }

	free(combo0); free(pop);
}

void generate_permutation(int* perm, int n, unsigned int seed) {
	// Initialize with 1 to n
	for(int i = 0; i < n; i++) {
		perm[i] = i + 1;
	}
	
	// Use a simple linear congruential generator for reproducible randomness
	unsigned int rand_state = seed;
	
	// Fisher-Yates shuffle
	for(int i = n - 1; i > 0; i--) {
		// Generate random number using LCG
		rand_state = rand_state * 1103515245 + 12345;
		int j = (rand_state / 65536) % (i + 1);
		
		// Swap elements i and j
		int temp = perm[i];
		perm[i] = perm[j];
		perm[j] = temp;
	}
}

void apply_perturbation_to_model(qnmodel* qnm, int perturbation_digit, int perturbation_seed, mpz_t scale_factor) {
	// Calculate scale factor: 10^d where d is the perturbation digit
	mpz_set_ui(scale_factor, 10);
	mpz_pow_ui(scale_factor, scale_factor, perturbation_digit);
	
	// For each class, generate a random permutation of 1 to M+1
	int* perm = (int*)calloc(qnm->M + 1, sizeof(int));
	
	for(int j = 0; j < qnm->R; j++) {
		// Generate permutation for class j using seed based on class index
		generate_permutation(perm, qnm->M + 1, perturbation_seed + j * 1000);
		
		// Apply perturbation to L[i][j] for all stations i
		for(int i = 0; i < qnm->M; i++) {
			mpz_mul(qnm->L[i][j], qnm->L[i][j], scale_factor);
			mpz_add_ui(qnm->L[i][j], qnm->L[i][j], perm[i]);
		}
		
		// Apply perturbation to Z[j] using the last element of permutation
		mpz_mul(qnm->Z[j], qnm->Z[j], scale_factor);
		mpz_add_ui(qnm->Z[j], qnm->Z[j], perm[qnm->M]);
	}
	
	free(perm);
}

int main(int argc, char** argv)
{
	int i, s, m, r, nr;
	bool out_e = false, out_g = false, out_l = false, out_t = false, out_q = false, validate = false;
	int pdigit = 0, pseed = 23000;
	char* model_file = NULL;

	t0 = CPUTIME;
	for (i = 1; i < argc; i++) {
		if      (!strcmp(argv[i], "-e") || !strcmp(argv[i], "--ex"))  out_e = true;
		else if (!strcmp(argv[i], "-g") || !strcmp(argv[i], "--nc"))  out_g = true;
		else if (!strcmp(argv[i], "-l") || !strcmp(argv[i], "--log")) out_l = true;
		else if (!strcmp(argv[i], "-t") || !strcmp(argv[i], "--tput")) out_t = true;
		else if (!strcmp(argv[i], "-q") || !strcmp(argv[i], "--qlen")) out_q = true;
		else if (!strcmp(argv[i], "--validate")) validate = true;
		else if (!strcmp(argv[i], "-p")) { if (i+1 < argc) pdigit = atoi(argv[++i]); }
		else if (!strcmp(argv[i], "-s") || !strcmp(argv[i], "--seed")) { if (i+1 < argc) pseed = atoi(argv[++i]); }
		else if (argv[i][0] != '-') model_file = argv[i];
	}
	if (!model_file) {
		printf("USAGE: %s [-e|-g|-l|-t|-q] [--validate] [-p digit] [-s seed] model.qn\n", argv[0]);
		printf("  gmom: generalized (divide-and-conquer, b=1) Method of Moments.\n");
		printf("  Outputs the top normalizing-constant basis entry V{M,l}(1),\n");
		printf("  matching the reference mbmom1.  --validate checks the whole\n");
		printf("  basis against exact convolution of the augmented models.\n");
		return -1;
	}

	qnmodel* qn = readmodel(model_file);
	int M = qn->M, R = qn->R;
	if (R < 2) { fprintf(stderr, "gmom requires at least two classes\n"); return 1; }
	{ int c; for (c = 0; c < R; c++) if (qn->N[c] < 0) {
		fprintf(stderr, "gmom: open/mixed models are not supported (class %d is open); use a closed model\n", c+1); return 1; } }
	if (qn->hasOpen) { fprintf(stderr, "gmom: open/mixed models are not supported\n"); return 1; }
	if (qn->isLD)    { fprintf(stderr, "gmom: load-dependent stations are not supported\n"); return 1; }
	{ int q; for (q = 0; q < M; q++) if (qn->mi[q] > 1) {
		fprintf(stderr, "gmom: multiserver/replicated stations (mi>1 at queue %d) are not supported; the b=1 basis assumes distinct single-server queues\n", q+1); return 1; } }
	mpz_t** L = qn->L; mpz_t* Z = qn->Z; int* N = qn->N;
	nckinit(M + R - 1, R);

	mpz_t scale_factor; mpz_init_set_ui(scale_factor, 1);
	if (pdigit > 0) apply_perturbation_to_model(qn, pdigit, pseed, scale_factor);

	/* Init-time singularity screen (loading-only, one factorisation per
	 * level).  If the recursion would hit a degenerate matrix, move the
	 * offending class to the recursion position (its loadings then leave
	 * the coefficient matrix); if no single-class reorder removes it, the
	 * model is genuinely degenerate for gmom's basis. */
	int* cperm = (int*) malloc(R * sizeof(int));
	for (i = 0; i < R; i++) cperm[i] = i;
	if (pfqn_recursion_singular(L, M, R)) {
		int fix = pfqn_nonsingular_recclass(L, M, R);
		if (fix < 0) {
			if (out_e) {
				fprintf(stderr, "gmom: model is singular under every class order (genuine loading degeneracy);\n  no exact G available. Use bin/ca or bin/safe_comom, or -g/-t/-q with -p for a perturbed approximation.\n");
				return 1;
			}
			/* auto-perturb at digit 20 (as mom/comom do) for an approximate answer */
			if (mpz_cmp_ui(scale_factor, 1) == 0) {
				fprintf(stderr, "gmom: singular system; automatically applying perturbation at digit 20 (approximate result).\n");
				apply_perturbation_to_model(qn, 20, pseed, scale_factor);
				L = qn->L; Z = qn->Z; N = qn->N;
			}
		} else {
		/* cperm[j] = original class placed at position j (fix goes last) */
		int c = 0, j;
		for (j = 0; j < R; j++) if (j != fix) cperm[c++] = j;
		cperm[R-1] = fix;
		/* apply the permutation to L columns, N and Z (work in place) */
		mpz_t** L2 = (mpz_t**) malloc(M * sizeof(mpz_t*));
		int k;
		for (k = 0; k < M; k++) { L2[k] = (mpz_t*) malloc(R * sizeof(mpz_t)); for (j = 0; j < R; j++) { mpz_init_set(L2[k][j], L[k][cperm[j]]); } }
		int* N2 = (int*) malloc(R * sizeof(int)); mpz_t* Z2 = (mpz_t*) malloc(R * sizeof(mpz_t));
		for (j = 0; j < R; j++) { N2[j] = N[cperm[j]]; mpz_init_set(Z2[j], Z[cperm[j]]); }
		L = L2; N = N2; Z = Z2;
		qn->L = L2; qn->N = N2; qn->Z = Z2;   /* measures reads from qn */
		}
	}

	int* nvec = (int*) calloc(R, sizeof(int));

	/* V slots for m=1..M */
	mpq_vec_t* vl  = (mpq_vec_t*) calloc(M+1, sizeof(mpq_vec_t));
	mpq_vec_t* vl1 = (mpq_vec_t*) calloc(M+1, sizeof(mpq_vec_t));
	int* lenl  = (int*) calloc(M+1, sizeof(int));
	int* lenl1 = (int*) calloc(M+1, sizeof(int));
	mpq_vec_t grM = NULL; int lengrM = 0;   /* top basis at N-1_R, for mdecrease */

	/* ================= class 1 initialisation (r=2) ================= */
	nvec[0] = N[0];
	{
		r = 2;
		int* pop = (int*) malloc(r * sizeof(int));
		for (m = 1; m <= M; m++) {
			/* slot l_1: multichoose(m, r-2=0) -> 1 combo of zeros */
			int szI0 = 1;
			int** I0 = sortbynnzpos(multichoose(m, 0), 1, m);
			vfree(vl1[m], lenl1[m]); lenl1[m] = szI0*r; vl1[m] = vnew(szI0*r);
			for (i = 0; i < szI0; i++) for (s = 0; s < r; s++) {
				oner(pop, nvec, r, s); base_G(vl1[m][i*r+s], L, m, I0[i], r, pop, Z);
			}
			free(I0);
			/* slot l: multichoose(m, r-1=1) */
			int szIk = nck(m + (r-1) - 1, m-1);
			int** Ik = sortbynnzpos(multichoose(m, 1), szIk, m);
			vfree(vl[m], lenl[m]); lenl[m] = szIk*r; vl[m] = vnew(szIk*r);
			for (i = 0; i < szIk; i++) for (s = 0; s < r; s++) {
				oner(pop, nvec, r, s); base_G(vl[m][i*r+s], L, m, Ik[i], r, pop, Z);
			}
			free(Ik);
		}
		free(pop);
	}

	/* ================= class loop r = 2..R ========================== */
	for (r = 2; r <= R; r++) {
		/* build matrices for m=2..M with the current nvec (nvec[r-1]=0 here) */
		LS1** LS = (LS1**) calloc(M+1, sizeof(LS1*));
		for (m = 2; m <= M; m++) LS[m] = setup1(L, m, r, nvec, Z);

		if (r > 2) {
			/* Initialise the new class exactly from convolution at the
			 * current nvec (class r at 0).  The reference mbmom1 uses a
			 * cheap component shift here, but that is only valid when the
			 * intermediate classes have no think time (the reference
			 * errors on Z<>0); the exact fill is correct for all Z and
			 * costs only R-1 convolutions, off the population sweep. */
			for (m = 1; m <= M; m++) fill_level_all(vl, vl1, lenl, lenl1, L, m, r, nvec, Z);
		}

		/* population sweep of class r */
		for (nr = 1; nr <= N[r-1]; nr++) {
			nvec[r-1] = nr;
			/* capture the top basis at N-1_R just before the final step */
			if (r == R && nr == N[R-1]) {
				vfree(grM, lengrM); lengrM = lenl[M]; grM = vnew(lenl[M]);
				fp_copy_vec(grM, vl[M], lenl[M]);
			}
			fill_level1(vl, vl1, lenl, lenl1, L, r, nvec, Z);
			for (m = 2; m <= M; m++) {
				int colA = LS[m]->colA, colC = LS[m]->colC, colD = LS[m]->colD, nrows = LS[m]->nrows;
				/* u = (E*vl[m] + F*vl1[m]) / nr  ; keep old vl[m] */
				mpq_vec_t ue = vnew(colD), uf = vnew(colD), u = vnew(colD);
				matvec_rect(ue, LS[m]->E, colD, colA, vl[m]);
				matvec_rect(uf, LS[m]->F, colD, colD, vl1[m]);
				mpq_t inv; mpq_init(inv); mpq_set_si(inv, 1, nr);
				for (i = 0; i < colD; i++) { mpq_add(u[i], ue[i], uf[i]); mpq_mul(u[i], u[i], inv); }
				mpq_clear(inv);
				vfree(vl1[m], lenl1[m]); vl1[m] = u; lenl1[m] = colD;
				vfree(ue, colD); vfree(uf, colD);
				/* rhs = B*vl[m] + C*vl[m-1] + D*vl1[m] */
				mpq_vec_t rb = vnew(nrows), rc = vnew(nrows), rd = vnew(nrows), rhs = vnew(nrows);
				matvec_rect(rb, LS[m]->B, nrows, colA, vl[m]);
				matvec_rect(rc, LS[m]->C, nrows, colC, vl[m-1]);
				matvec_rect(rd, LS[m]->D, nrows, colD, vl1[m]);
				for (i = 0; i < nrows; i++) { mpq_add(rhs[i], rb[i], rc[i]); mpq_add(rhs[i], rhs[i], rd[i]); }
				mpq_vec_t x = solve_normal(LS[m]->A, nrows, colA, rhs);
				vfree(vl[m], lenl[m]); vl[m] = x; lenl[m] = colA;
				vfree(rb, nrows); vfree(rc, nrows); vfree(rd, nrows); vfree(rhs, nrows);
			}
		}
		for (m = 2; m <= M; m++) ls1_free(LS[m], r);
		free(LS);
	}

	/* ===================== output / validation ===================== */
	if (validate) {
		int* pop = (int*) malloc(R * sizeof(int));
		mpq_t ref; mpq_init(ref);
		int total = 0, totbad = 0;
		for (m = 1; m <= M; m++) {
			int szIk = nck(m + (R-1) - 1, m-1);
			int** Ik = sortbynnzpos(multichoose(m, R-1), szIk, m);
			int badm = 0;
			for (i = 0; i < szIk; i++) for (s = 0; s < R; s++) {
				oner(pop, N, R, s);
				base_G(ref, L, m, Ik[i], R, pop, Z);
				total++;
				if (mpq_cmp(ref, vl[m][i*R+s]) != 0) {
					badm++; totbad++;
					if (badm <= 2) gmp_fprintf(stderr, "  m=%d i=%d s=%d gmom=%Qd conv=%Qd\n", m, i, s, vl[m][i*R+s], ref);
				}
			}
			free(Ik);
		}
		mpq_clear(ref); free(pop);
		if (totbad) { printf("VALIDATE: %d/%d WRONG\n", totbad, total); return 1; }
		printf("VALIDATE: all %d basis entries match exact convolution\n", total);
		return 0;
	}

	/* plain G(N), X, Q via the mdecrease replica-descent */
	if (!grM) { fprintf(stderr, "gmom: internal error, grM not captured\n"); return 2; }
	gmom_measures(qn, vl[M], grM, out_e, out_g, out_l, out_t, out_q, cperm, scale_factor);

	return 0;
}
