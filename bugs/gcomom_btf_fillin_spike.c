/* gcomom spike (PLAN_extensions.md Step 4 gate).
 *
 * Question: gmom solves each level's overdetermined system by exact normal
 * equations (A^T A) x = A^T b.  CoMoM's per-class matrix A11 is upper block
 * triangular (comom/btf_decompose.c) and that is where its complexity
 * advantage comes from.  Does the block-triangular zero pattern survive the
 * normal equations?
 *
 * This measures, on comom's own A11 for a given model and class:
 *   nnz(A11), nnz(A11^T A11)
 *   nnz below the block diagonal in A11 (must be 0 by BTF)
 *   nnz below the block diagonal in A11^T A11 (the fill-in)
 */
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <gmpla.h>
#include "comom.h"

double t0, t1;
int USE_LINBOX = 0;

int main(int argc, char** argv)
{
	if (argc < 2) { printf("usage: %s model.qn\n", argv[0]); return 1; }
	qnmodel* qn = readmodel(argv[1]);
	int M = qn->M, R = qn->R, r;
	nckinit(M + R + 1, M);

	for (r = 2; r <= R; r++) {
		combsrep* Dn = (combsrep*) calloc(1, sizeof(combsrep));
		Dn->n = r; Dn->k = M;
		Dn->card = nck(M + r - 1, M);
		Dn->combs = (int**) multichoose(r, M);
		{ int d; for (d = 0; d < Dn->card; d++) Dn->combs[d][r-1] = 0; }
		Dn->combs = (int**) sortbynnzpos(Dn->combs, Dn->card, r);

		long cardG  = nck(M + r - 1, M);
		long cardGk = cardG * M;
		int* N = (int*) int_vec(R, 0);
		{ int c; for (c = 0; c < R; c++) N[c] = qn->N[c]; }

		LS* ls = setupls(Dn, qn, N, r, 0.0, 0);
		btf_info* btf = btf_decompose(ls->A11, Dn, r, M, cardGk);
		if (!btf) { printf("class %d: btf_decompose failed\n", r); continue; }

		/* block id per index */
		int* blk = (int*) calloc(cardGk, sizeof(int));
		int b, i, j;
		for (b = 0; b < btf->num_blocks; b++)
			for (i = btf->block_starts[b]; i < btf->block_starts[b+1]; i++)
				blk[i] = b;

		long nnzA = 0, belowA = 0;
		for (i = 0; i < cardGk; i++)
			for (j = 0; j < cardGk; j++)
				if (mpq_sgn(ls->A11[i][j])) { nnzA++; if (blk[i] > blk[j]) belowA++; }

		mpq_mat_t AtA = mpq_matzeros(cardGk, cardGk);
		mpq_mattransmul(AtA, ls->A11, ls->A11, cardGk, cardGk);

		long nnzG = 0, belowG = 0;
		for (i = 0; i < cardGk; i++)
			for (j = 0; j < cardGk; j++)
				if (mpq_sgn(AtA[i][j])) { nnzG++; if (blk[i] > blk[j]) belowG++; }

		long belowcells = 0;
		for (i = 0; i < cardGk; i++)
			for (j = 0; j < cardGk; j++)
				if (blk[i] > blk[j]) belowcells++;

		printf("class r=%d  size=%ld  blocks=%d\n", r, cardGk, btf->num_blocks);
		printf("  A11    : nnz=%6ld (%.1f%% dense)   below-diagonal-block nnz=%ld\n",
		       nnzA, 100.0*nnzA/(double)(cardGk*cardGk), belowA);
		printf("  A11^T*A: nnz=%6ld (%.1f%% dense)   below-diagonal-block nnz=%ld of %ld cells (%.1f%% filled)\n",
		       nnzG, 100.0*nnzG/(double)(cardGk*cardGk), belowG, belowcells,
		       belowcells ? 100.0*belowG/(double)belowcells : 0.0);
		free(blk);
	}
	return 0;
}
