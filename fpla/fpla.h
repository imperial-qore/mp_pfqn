#ifndef FPLA_H
#define FPLA_H

/* fpla - fixed-precision floating-point linear algebra.
 *
 * Mirrors the gmpla mpq_* API one-for-one, but every quantity is an MPFR
 * float carried at a run-time settable working precision FP_WPREC, with
 * round-to-nearest.  At FP_WPREC=53 the arithmetic is bit-for-bit IEEE
 * double precision except for the exponent range, which MPFR leaves
 * effectively unbounded.  That separation is deliberate: it lets us
 * measure the loss caused by finite *precision* alone, without the
 * overflow of normalizing constants confounding the experiment.
 *
 * A second precision FP_RPREC is used only for residual evaluation in
 * Wilkinson/Moler iterative refinement.
 */

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <mpfr.h>
#include <math.h>

/* integer helpers, reused from gmpla (same objects are linked in) */
int** int_mat(int rows, int cols, int initval);
long int int_matmatchrow(int** mat, long int nrows, long int ncols,int* row);
void int_matprint(int** mat, int n, int m);
int* int_vec(int n, int value);
void int_vecprint(int* vec, int n);

extern mpfr_prec_t FP_WPREC; /* working precision, bits */
extern mpfr_prec_t FP_RPREC; /* residual precision, bits */

typedef mpfr_t  fp_t;
typedef mpfr_t* fp_vec_t;
typedef mpfr_t** fp_mat_t;

typedef struct struct_fp_msp
{
	fp_vec_t coeff;
	int* pos_row;
	int* pos_col;
	int rows;
	int cols;
	int nnz;
	int lastnnz;
} *fp_msp_t;

/* --- scalar operations: MPFR with round-to-nearest ------------------- */
#define fp_init(x)        mpfr_init2((x),FP_WPREC)
#define fp_clear(x)       mpfr_clear(x)
#define fp_set(r,a)       mpfr_set((r),(a),MPFR_RNDN)
#define fp_set_d(r,d)     mpfr_set_d((r),(d),MPFR_RNDN)
#define fp_set_z(r,z)     mpfr_set_z((r),(z),MPFR_RNDN)
#define fp_add(r,a,b)     mpfr_add((r),(a),(b),MPFR_RNDN)
#define fp_sub(r,a,b)     mpfr_sub((r),(a),(b),MPFR_RNDN)
#define fp_mul(r,a,b)     mpfr_mul((r),(a),(b),MPFR_RNDN)
#define fp_div(r,a,b)     mpfr_div((r),(a),(b),MPFR_RNDN)
#define fp_neg(r,a)       mpfr_neg((r),(a),MPFR_RNDN)
#define fp_abs(r,a)       mpfr_abs((r),(a),MPFR_RNDN)
#define fp_cmp(a,b)       mpfr_cmp((a),(b))
#define fp_sgn(a)         mpfr_sgn(a)
#define fp_get_d(a)       mpfr_get_d((a),MPFR_RNDN)

void fp_set_si(fp_t x, long num, long den);
void fp_set_ui(fp_t x, unsigned long num, unsigned long den);

/* --- containers ------------------------------------------------------ */
fp_vec_t fp_vec(int rows, int valuenum, int valueden);
fp_mat_t fp_matzeros(int rows, int cols);
void     fp_matfree(fp_mat_t m, int rows, int cols);
fp_msp_t fp_msp(int rows, int cols, int numnnz);
fp_t*    fp_mspget(fp_msp_t msp, int row, int col);
void     fp_matset_si(fp_mat_t mat, int row, int col, int coeffnum, int coeffden);
void     fp_mspset_si(fp_msp_t msp, int row, int col, int coeffnum, int coeffden);
void     fp_mspset_z(fp_msp_t msp, int row, int col, mpz_t coeffnum, int coeffden);
void     fp_mspvecmul(fp_vec_t rop, fp_msp_t op1, fp_vec_t op2);
void     fp_vecdup(fp_vec_t v1, fp_vec_t v2, int n);
void     fp_vecprint(fp_vec_t vec, int n);
void     fp_matprint(fp_mat_t mat, int n, int m);
void     fp_mspprint(fp_msp_t msp);

/* --- dense LU -------------------------------------------------------- */
int*     fp_ludcmp(fp_mat_t A, int N);
int      fp_lubksb(fp_mat_t A, fp_vec_t b, int N, int* indx);

/* --- Wilkinson/Moler iterative refinement ---------------------------- */

/* Duplicate a square matrix (used to keep a pristine copy of the block
 * before fp_ludcmp overwrites it with its factors). */
fp_mat_t fp_matdup(fp_mat_t A, int N);

/* r = b - A*x, accumulated at precision FP_RPREC then rounded back to
 * the working precision.  This is the only step of refinement that must
 * be carried above the working precision: it is where the cancellation
 * lives. */
void fp_residual(fp_vec_t r, fp_mat_t A, fp_vec_t x, fp_vec_t b, int N);

/* Solve A x = b to working-precision accuracy.  Aorig is the pristine
 * matrix, LU/indx its factorization at working precision.  On entry b
 * holds the right-hand side; on exit x holds the refined solution.
 * Performs at most maxit refinement sweeps, stopping early when the
 * correction is below tol in relative infinity norm.  Returns the number
 * of sweeps actually performed, or -1 on a singular factor.
 * *lastrel, if non-NULL, receives the final relative correction norm. */
double fp_cond_inf(fp_mat_t Aorig, fp_mat_t LU, int* indx, int N);

int fp_lurefine(fp_mat_t Aorig, fp_mat_t LU, int* indx, fp_vec_t b,
                fp_vec_t x, int N, int maxit, double tol, double* lastrel);

#endif
