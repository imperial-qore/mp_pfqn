#ifndef PROMOM_H
#define PROMOM_H

#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <gmp.h>
#include <gmpla.h>
#include "util.h"

#define CPUTIME (getrusage(RUSAGE_SELF,&ruse), ruse.ru_utime.tv_sec + ruse.ru_stime.tv_sec + 1e-6 * (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec))

/* promom - marginal queue-length probabilities on the MoM basis.
 *
 * procomom does this on the CoMoM basis, whose elements are indexed by a
 * class-population shift and a station.  The MoM basis is indexed instead by
 * a station-REPLICA combination and a single class decrement, so the port is
 * a different index map and different CE/PC rows, not a different algorithm.
 *
 * Quantity.  Fix a reference station (rotated to be station M).  For the
 * model augmented with the replica vector d (station j has mi_j + d_j
 * servers), at population P, write
 *
 *   q_d(P, n) = sum_{|k|=n} (n!/prod_r k_r!) prod_r L_{M,r}^{k_r} g_{d,-M}(P-k)
 *
 * the unnormalized probability that the reference copy of station M holds
 * exactly n jobs; g_{d,-M} is the constant of the same model with that copy
 * removed.  Summing over n gives the ordinary constant, q_d(P,.) is what the
 * MoM basis carries here in place of g_d(P), and
 *
 *   P(n_M = n) = q(N, n) / G(N).
 *
 * Three identities close the recursion (derivations in README.md):
 *
 *   (n-recursion)  q_d(P, n) = sum_r L_{M,r} q_d(P - 1_r, n-1),   n >= 1
 *   (CE)           q_d(P, n) = q_{d-1_j}(P, n) + sum_r L_{j,r} q_d(P - 1_r, n)
 *   (PC)           P_s q_d(P,n) - n L_{M,s} q_d(P-1_s, n-1)
 *                    = Z_s q_d(P-1_s, n)
 *                    + sum_j m'_j L_{j,s} q_{d+1_j}(P-1_s, n),
 *                  m'_j = mi_j + d_j - [j = M]
 *
 * The CE is mom's convolution equation verbatim (a replica of station j does
 * not change the count at the reference copy).  The PC is mom's population
 * constraint plus one n-coupled term -- exactly the role of procomom's DC/DD
 * blocks -- and with the reference copy held out of the station sum, which
 * is what makes q(P,0) the constant of the model WITHOUT that copy rather
 * than G(P).  So the per-step system has the same shape as
 * procomom's,
 *
 *   A q(n) = B q_prev(n) + n [ DC q(n-1) + DD q_prev(n-1) ],
 *
 * only assembled over the MoM basis.  Setting n = 0 everywhere recovers
 * mom's own system, which is the check that the port is faithful.
 */

/* Combinations with repetition (independent from mom/comom). */
typedef struct
{
	int n;       /* vector length (M stations)                */
	int k;       /* total to distribute (replica level)       */
	int card;    /* nchoosek(n+k-1, k)                        */
	int **combs; /* card x n, sorted by nnz position          */
} qcombs;

/* The four dense rectangular blocks of one population step.
 * Columns are the MoM basis at the current point:
 *     [0, cardk*r)                level-r part,   index w*r + s
 *     [cardk*r, (cardk+cardi)*r)  level-(r-1) part, index (cardk+i)*r + s
 * with s = 0 meaning P and s >= 1 meaning P - 1_s. */
typedef struct
{
	mpq_mat_t A;    /* current point,  n     */
	mpq_mat_t B;    /* previous point, n     */
	mpq_mat_t DC;   /* current point,  n-1   */
	mpq_mat_t DD;   /* previous point, n-1   */
	int nrows;
	int ncols;
} QMatrices;

extern struct rusage ruse;
extern double t0, t1;

/* combination helpers */
qcombs* qcombs_new(int nvec, int level);
void    qcombs_free(qcombs* c);
int     qcombs_index(qcombs* c, int* comb);

/* per-step system over the MoM basis; Ncur is the population at the current
 * point, r the class being swept, refL the reference station's demand row */
QMatrices* genqmatrix(qcombs* Ik, qcombs* I, qnmodel* qnm, int* Ncur, int r);
void       free_qmatrices(QMatrices* qm);

/* replica descent: q at level R/R-1 down to level 0, for every n */
int qdecrease(qnmodel* qnm, mpq_vec_t* gk, mpq_vec_t* g, mpq_vec_t* grv,
              int sumN, mpq_t* dist);

#endif
