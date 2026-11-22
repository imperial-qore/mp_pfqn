#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <gmp.h>
#include <zpla.h>
#include "comommod.h"

/* One modular image of the whole CoMoM recursion.
 *
 * A transliteration of the class loop in comom/main.c: same basis, same
 * order of classes and populations, same linear systems, with mpq_t
 * replaced by a residue in a machine word.  What is dropped is the
 * progress reporting, the perturbation machinery (an image is exact by
 * construction and a singular system is reported rather than perturbed)
 * and the dead overdetermined path guarded by n_nf_combos, which
 * comom/main.c never enables.
 *
 * Returns 0 on success and -1 when the prime is unlucky, meaning some
 * matrix that is non-singular over Q became singular modulo p.  The
 * caller distinguishes an unlucky prime from a genuine singularity by
 * trying further primes: a matrix singular over Q is singular modulo
 * every prime.
 */
int comom_image(qnmodel* qnm, combsrep** Dn_all, zp_t p, zp_t* res)
{
	int r,s,d,i,k,nr;
	int M = qnm->M, R = qnm->R;
	int rc = 0;
	int* n = NULL;
	combsrep *Dn = NULL, *Dn_old = NULL;
	zp_vec_t g=NULL, b1=NULL, b1b=NULL, G=NULL, g_prev=NULL;
	zp_vec_t* g_m = NULL;
	zp_t* marginal_G = NULL;
	zp_t** marginal_Gk = NULL;
	LS* linsys = NULL;
	long int max_cardG, max_size;
	long int cardG = 0, cardGk = 0;
	long int finalCardGk;

	zp_setmod(p);

	/* g_m holds the last M+1 bases of the previous class, sized once for
	 * the largest basis so that no class needs to reallocate. */
	max_cardG = nck(M+R-1,M);
	max_size  = max_cardG + max_cardG*M;

	g_m = (zp_vec_t*) calloc((size_t)M+1, sizeof(zp_vec_t));
	for (i=0;i<=M;i++) g_m[i] = zp_vec(max_size,0,1);

	marginal_G  = (zp_t*)  calloc((size_t)R, sizeof(zp_t));
	marginal_Gk = (zp_t**) calloc((size_t)M, sizeof(zp_t*));
	for (k=0;k<M;k++) marginal_Gk[k] = (zp_t*) calloc((size_t)R, sizeof(zp_t));

	for (r=1;r<=R;r++)
	{
		Dn = Dn_all[r];

		free(n);
		n = (int*) int_vec(R,0);
		for (s=1;s<=r-1;s++) n[s-1] = qnm->N[s-1];

		cardG  = nck(M+r-1,M);
		cardGk = cardG*M;

		/* comom leaks these once per class; an image runs once per prime,
		 * so they are released here instead. */
		zp_vecfree(g); zp_vecfree(b1); zp_vecfree(b1b); zp_vecfree(G);
		g   = zp_vec(cardG+cardGk,0,1);
		b1  = zp_vec(cardGk,0,1);
		b1b = zp_vec(cardGk,0,1);
		G   = zp_vec(cardG,0,1);

		/* initialize g */
		if (r==1)
		{
			int* comb = (int*) int_vec(r,0);
			for (i=0;i<=M;i++) g[hash(Dn,comb,i)-1] = zp_ratio_si(1,1);
			free(comb);
		}
		else
		{
			int* temp_comb = (int*) calloc((size_t)r-1, sizeof(int));
			for (d=1;d<=Dn->card;d++)
				for (i=0;i<=M;i++)
				{
					int col = hash(Dn,Dn->combs[d-1],i)-1;
					int col_old, j;
					s = Dn->combs[d-1][r-2];
					for (j=0;j<r-1;j++) temp_comb[j] = Dn->combs[d-1][j];
					temp_comb[r-2] = 0;
					col_old = hash(Dn_old,temp_comb,i)-1;
					g[col] = g_m[s][col_old];
				}
			free(temp_comb);
		}

		/* When N_r <= M the inner loop does not produce enough saved
		 * bases, so the initial one stands in for population n_r = 0. */
		if (r < R && qnm->N[r-1] <= M)
		{
			zp_vecdup(g_m[qnm->N[r-1]], g, cardG+cardGk);
			Dn_old = Dn;
		}

		for (nr=1;nr<=qnm->N[r-1];nr++)
		{
			n[r-1] = nr;

			if (nr==1)
			{
				freels(linsys);
				linsys = setupls(Dn,qnm,n,r);
				linsys->btf = btf_decompose(linsys->A11,Dn,r,M,(int)cardGk);
				if (!linsys->btf || btf_factorize(linsys->btf) < 0)
				{
					/* A block that is singular mod p need not be singular
					 * over Q, so fall back to a direct LU of A11 before
					 * giving up on this prime. */
					if (linsys->btf) { btf_free(linsys->btf); linsys->btf = NULL; }
					linsys->A11_lu = (zp_mat_t) zp_matzeros(cardGk,cardGk);
					for (i=0;i<cardGk;i++)
						memcpy(linsys->A11_lu[i], linsys->A11[i], (size_t)cardGk*sizeof(zp_t));
					linsys->A11_lu_indices = zp_ludcmp(linsys->A11_lu,(int)cardGk);
					if (linsys->A11_lu_indices == NULL) { rc = -1; goto cleanup; }
				}
			}

			/* G = B2*g/n_r.  n_r < p is guaranteed by the prime floor, so
			 * the inverse always exists. */
			zp_mspvecmul(G, linsys->B2, g);
			{
				zp_t ninv = zp_invv(zp_from_si(nr));
				for (i=0;i<cardG;i++) G[i] = zp_mulv(ninv, G[i]);
			}

			/* zero G at impossible positions before forming the RHS */
			if (r >= 2)
				for (d=0; d<Dn->card; d++)
				{
					int is_imp = 0;
					for (s=0; s<=r-2; s++)
						if (Dn->combs[d][s] > qnm->N[s]) { is_imp = 1; break; }
					if (is_imp) G[d] = 0;
				}

			/* b1 = B1*g - A12*G */
			zp_mspvecmul(b1,  linsys->B1,  g);
			zp_mspvecmul(b1b, linsys->A12, G);
			for (i=0;i<cardGk;i++) b1[i] = zp_subv(b1[i], b1b[i]);

			/* solve A11 Gk = b1 */
			if (linsys->btf)
			{
				if (btf_solve(linsys->btf, b1, (int)cardGk) < 0) { rc = -1; goto cleanup; }
			}
			else
			{
				if (zp_lubksb(linsys->A11_lu, b1, (int)cardGk, linsys->A11_lu_indices) < 0)
					{ rc = -1; goto cleanup; }
			}

			/* Fix impossible-combo contamination: re-solve on the
			 * sub-system spanned by the reachable combinations only. */
			if (r >= 2)
			{
				int has_impossible = 0;
				int n_possible = 0;
				int* combo_possible = (int*) calloc((size_t)Dn->card, sizeof(int));

				for (d=1; d<=Dn->card; d++)
				{
					int is_imp = 0;
					for (s=0; s<=r-2; s++)
						if (Dn->combs[d-1][s] > qnm->N[s]) { is_imp = 1; break; }
					combo_possible[d-1] = !is_imp;
					if (!is_imp) n_possible++; else has_impossible = 1;
				}

				if (has_impossible)
				{
					int n_poss_Gk = n_possible * M;
					int* col_map = (int*) calloc((size_t)n_poss_Gk, sizeof(int));
					int* row_included = (int*) calloc((size_t)cardGk, sizeof(int));
					int n_sub_rows = 0;
					int ci = 0, rank_ok = 1;
					zp_vec_t b1_rhs, a12g, b_sub, x_sub;
					zp_mat_t A_sub;

					for (d=0; d<Dn->card; d++)
					{
						if (!combo_possible[d]) continue;
						for (k=1;k<=M;k++) col_map[ci++] = hash(Dn,Dn->combs[d],k)-1;
					}

					for (d=0; d<Dn->card; d++)
						if (!combo_possible[d]) G[d] = 0;

					/* RHS recomputed with G zeroed at the impossible slots */
					b1_rhs = zp_vec(cardGk,0,1);
					a12g   = zp_vec(cardGk,0,1);
					zp_mspvecmul(b1_rhs, linsys->B1,  g);
					zp_mspvecmul(a12g,   linsys->A12, G);
					for (i=0;i<cardGk;i++) b1_rhs[i] = zp_subv(b1_rhs[i], a12g[i]);

					{
						int row = 0;
						for (d=1; d<=Dn->card; d++)
						{
							int dp;
							if (int_vecsubsum(Dn->combs[d-1], 0, r-1) >= M) continue;
							dp = combo_possible[d-1];
							for (k=0;k<M;k++)
							{
								if (dp) { row_included[row] = 1; n_sub_rows++; }
								row++;
							}
							for (s=0;s<r-1;s++)
							{
								if (dp && Dn->combs[d-1][s] + 1 <= qnm->N[s])
									{ row_included[row] = 1; n_sub_rows++; }
								row++;
							}
						}
					}

					A_sub = (zp_mat_t) zp_matzeros(n_sub_rows, n_poss_Gk);
					b_sub = zp_vec(n_sub_rows,0,1);
					{
						int sr = 0, orig_row, j;
						for (orig_row=0; orig_row<cardGk; orig_row++)
						{
							if (!row_included[orig_row]) continue;
							for (j=0;j<n_poss_Gk;j++)
								A_sub[sr][j] = linsys->A11[orig_row][col_map[j]];
							b_sub[sr] = b1_rhs[orig_row];
							sr++;
						}
					}

					/* Gaussian elimination, pivoting on the first non-zero.
					 * A pivot that is missing modulo p but present over Q
					 * would silently change the rank and hence the
					 * particular solution, so the image is discarded
					 * instead: over Q the rank is at least the rank mod p,
					 * and a rank that is genuinely deficient fails for
					 * every prime and is reported by the caller. */
					{
						int col, row, j;
						for (col=0; col<n_poss_Gk; col++)
						{
							int piv = -1;
							for (row=col; row<n_sub_rows; row++)
								if (A_sub[row][col] != 0) { piv = row; break; }
							if (piv < 0) { rank_ok = 0; break; }
							if (piv != col)
							{
								zp_t* t = A_sub[col]; A_sub[col] = A_sub[piv]; A_sub[piv] = t;
								{ zp_t w = b_sub[col]; b_sub[col] = b_sub[piv]; b_sub[piv] = w; }
							}
							for (row=col+1; row<n_sub_rows; row++)
							{
								zp_t fac;
								if (A_sub[row][col] == 0) continue;
								fac = zp_divv(A_sub[row][col], A_sub[col][col]);
								for (j=col;j<n_poss_Gk;j++)
									A_sub[row][j] = zp_subv(A_sub[row][j], zp_mulv(fac,A_sub[col][j]));
								b_sub[row] = zp_subv(b_sub[row], zp_mulv(fac,b_sub[col]));
							}
						}
					}

					x_sub = zp_vec(n_poss_Gk,0,1);
					if (rank_ok)
					{
						int col, j;
						for (col=n_poss_Gk-1; col>=0; col--)
						{
							zp_t acc = b_sub[col];
							for (j=col+1;j<n_poss_Gk;j++)
								if (A_sub[col][j] != 0)
									acc = zp_subv(acc, zp_mulv(A_sub[col][j], x_sub[j]));
							x_sub[col] = zp_divv(acc, A_sub[col][col]);
						}
						for (j=0;j<n_poss_Gk;j++) b1[col_map[j]] = x_sub[j];
						for (d=0; d<Dn->card; d++)
							if (!combo_possible[d])
								for (k=1;k<=M;k++) b1[hash(Dn,Dn->combs[d],k)-1] = 0;
					}

					zp_matfree(A_sub, n_sub_rows);
					zp_vecfree(b_sub); zp_vecfree(x_sub);
					zp_vecfree(b1_rhs); zp_vecfree(a12g);
					free(col_map);
					free(row_included);

					if (!rank_ok) { free(combo_possible); rc = -1; goto cleanup; }
				}
				free(combo_possible);
			}

			/* G_1: the basis one population short of N, used for class R */
			if (r == R && nr == qnm->N[r-1])
			{
				if (g_prev == NULL) g_prev = zp_vec(cardG+cardGk,0,1);
				zp_vecdup(g_prev, g, cardG+cardGk);
			}

			for (i=0;i<cardGk;i++) g[i] = b1[i];
			for (i=0;i<cardG;i++)  g[cardGk+i] = G[i];

			/* keep the last M+1 bases for the next class */
			if (r < R && qnm->N[r-1]-nr <= M)
			{
				zp_vecdup(g_m[qnm->N[r-1]-nr], g, cardG+cardGk);
				Dn_old = Dn;
			}

			/* marginal constants for the performance measures */
			if (r == R && nr == qnm->N[r-1]-1)
			{
				marginal_G[R-1] = g[cardGk];
				for (k=0;k<M;k++) marginal_Gk[k][R-1] = g[k];
			}
			else if (r == R && nr == qnm->N[r-1])
			{
				int* comb = (int*) calloc((size_t)r, sizeof(int));
				for (s=1;s<=r;s++) comb[s-1] = Dn->combs[0][s-1];
				for (s=1;s<=R-1;s++)
				{
					int index;
					comb[s-1]++;
					index = hash(Dn,comb,0)-1;
					marginal_G[s-1] = g[index];
					for (k=1;k<=M;k++)
					{
						index = hash(Dn,comb,k)-1;
						marginal_Gk[k-1][s-1] = g[index];
					}
					comb[s-1]--;
				}
				free(comb);
			}
		}

		/* N_R = 0: the loop above never ran.  g is still the class-R basis
		 * at n_R = 0, i.e. G(N - comb), which is what the extraction above
		 * reads at nr = N_R. */
		if (r == R && qnm->N[r-1] == 0)
		{
			int* comb = (int*) calloc((size_t)r, sizeof(int));
			for (s=1;s<=r;s++) comb[s-1] = Dn->combs[0][s-1];
			for (s=1;s<=R-1;s++)
			{
				int index;
				comb[s-1]++;
				index = hash(Dn,comb,0)-1;
				marginal_G[s-1] = g[index];
				for (k=1;k<=M;k++)
				{
					index = hash(Dn,comb,k)-1;
					marginal_Gk[k-1][s-1] = g[index];
				}
				comb[s-1]--;
			}
			free(comb);
		}
		/* a class with no jobs: G(N-e_s) is an unreachable population */
		if (r == R)
			for (s=1;s<=R-1;s++)
				if (qnm->N[s-1] == 0)
				{
					marginal_G[s-1] = 0;
					for (k=0;k<M;k++) marginal_Gk[k][s-1] = 0;
				}
	}

	if (ZP_SINGULAR) { rc = -1; goto cleanup; }

	finalCardGk = nck(M+R-1,M)*M;
	res[0] = g[finalCardGk];
	for (r=1;r<=R;r++)
		res[r] = (r<R) ? marginal_G[r-1]
		               : (g_prev ? g_prev[finalCardGk] : marginal_G[r-1]);
	for (k=1;k<=M;k++)
		for (r=1;r<=R;r++)
			res[1+R+(k-1)*R+(r-1)] = (r<R) ? marginal_Gk[k-1][r-1]
			                               : (g_prev ? g_prev[k-1] : marginal_Gk[k-1][r-1]);
	for (i=0;i<NRES(qnm);i++) res[i] = zp_from_mont(res[i]);

cleanup:
	freels(linsys);
	zp_vecfree(g); zp_vecfree(b1); zp_vecfree(b1b); zp_vecfree(G);
	zp_vecfree(g_prev);
	if (g_m) { for (i=0;i<=M;i++) zp_vecfree(g_m[i]); free(g_m); }
	if (marginal_Gk) { for (k=0;k<M;k++) free(marginal_Gk[k]); free(marginal_Gk); }
	free(marginal_G);
	free(n);
	return rc;
}
