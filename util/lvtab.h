#ifndef LVTAB_H
#define LVTAB_H

/* Replica-level tables for the full MoM basis (Casale, Perf. Eval. 68
 * (2011), Def. 1 and App. C).
 *
 * While class r is processed the MoM recursion holds, for every multiset
 * d of extra queue replicas with |d| = l (l = 0..r), the r values
 *   G^{+d}(n), G^{+d}(n-1_1), ..., G^{+d}(n-1_{r-1}).
 * Levels r and r-1 are the two levels the linear system works on; levels
 * 0..r-2 are carried alongside and advanced by the class-r population
 * constraint
 *   n_r G^{+d}(n) = Z_r G^{+d}(n-1_r) + sum_j L_jr (m_j+d_j) G^{+d+1_j}(n-1_r),
 * which holds for any think time.  Levels 0 and 1 at N then give G(N),
 * G(N-1_s) and G^{+k}(N-1_s) directly.  The earlier downward pass
 * (mdecrease) instead used, for Z_s = 0, an identity that holds only when
 * no class has a positive think time, so it was wrong for models mixing
 * zero and positive think times.
 *
 * This header holds the arithmetic-independent part: the multisets of
 * each level, in the sortbynnzpos(multichoose()) order used by setupls,
 * and the map from (multiset, station) to the index of multiset+1_j in
 * the next level. */

#include <stdlib.h>
#include "util.h"

long int int_matmatchrow(int** mat, long int nrows, long int ncols, int* row);

typedef struct {
	int R, m;
	long *card;   /* card[l] = number of multisets of size l, l=0..R */
	int ***cb;    /* cb[l][i] = i-th multiset of size l */
	int **up;     /* up[l][i*m+j] = index of cb[l][i]+1_j in level l+1 */
} lvtab;

static inline lvtab* lvtab_new(int m, int R)
{
	int l, i, j;
	lvtab* T = (lvtab*)malloc(sizeof(lvtab));
	T->R = R; T->m = m;
	T->card = (long*)malloc((R+1)*sizeof(long));
	T->cb = (int***)malloc((R+1)*sizeof(int**));
	T->up = (int**)malloc((R+1)*sizeof(int*));
	for (l=0;l<=R;l++) {
		T->card[l] = nck(m+l-1,l);
		T->cb[l] = sortbynnzpos(multichoose(m,l),T->card[l],m);
	}
	for (l=0;l<R;l++) {
		T->up[l] = (int*)malloc(T->card[l]*m*sizeof(int));
		for (i=0;i<T->card[l];i++)
			for (j=0;j<m;j++) {
				T->cb[l][i][j]++;
				T->up[l][i*m+j] = int_matmatchrow(T->cb[l+1],T->card[l+1],m,T->cb[l][i]);
				T->cb[l][i][j]--;
			}
	}
	T->up[R] = NULL;
	return T;
}

static inline void lvtab_free(lvtab* T)
{
	int l, i;
	if (!T) return;
	for (l=0;l<=T->R;l++) {
		for (i=0;i<T->card[l];i++) free(T->cb[l][i]);
		free(T->cb[l]);
		free(T->up[l]);
	}
	free(T->cb); free(T->up); free(T->card); free(T);
}

#endif
