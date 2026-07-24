#ifndef ZPLA_H
#define ZPLA_H

/* zpla - linear algebra over the prime field Z/p, in machine words.
 *
 * Mirrors the gmpla mpq_* and fpla fp_* APIs so that the MoM recursion
 * can be run verbatim in a modular image.  Every value is a residue in
 * [0,p) held in a single 64-bit word, so no operation ever grows an
 * operand: the digit growth that dominates the exact rational solver is
 * moved out of the inner loop entirely and into the number of primes,
 * which is embarrassingly parallel.
 *
 * The modulus is thread-local so that independent primes can be run
 * concurrently in one process.
 *
 * Primes must satisfy p > max(N_r) so that the division by n_r in the
 * population recursion is invertible, and p < 2^62 so that the
 * __uint128_t product does not overflow.  A prime for which some
 * diagonal block is singular mod p ("unlucky prime") is detected and
 * discarded by the caller.
 */

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

typedef unsigned long long zp_t;
typedef zp_t*  zp_vec_t;
typedef zp_t** zp_mat_t;

typedef struct struct_zp_msp
{
	zp_vec_t coeff;
	int* pos_row;
	int* pos_col;
	int rows;
	int cols;
	int nnz;
	int lastnnz;
} *zp_msp_t;

/* integer helpers, reused from gmpla (same objects are linked in) */
int** int_mat(int rows, int cols, int initval);
long int int_matmatchrow(int** mat, long int nrows, long int ncols,int* row);
void int_matprint(int** mat, int n, int m);
int* int_vec(int n, int value);
void int_vecprint(int* vec, int n);

/* Residues are held in Montgomery form (a*2^64 mod p) so that the modular
 * multiply in the inner loop is three 64x64 multiplies rather than a
 * 128-by-64 hardware division, which costs an order of magnitude more.
 * Addition and subtraction are unaffected by the representation.  Values
 * cross the boundary through zp_from_z/zp_from_si on the way in and
 * zp_from_mont on the way out. */
extern __thread zp_t ZP_P;        /* current modulus, odd */
extern __thread zp_t ZP_PINV;     /* -p^{-1} mod 2^64 */
extern __thread zp_t ZP_R2;       /* (2^64)^2 mod p */
extern __thread int  ZP_SINGULAR; /* set when an inverse of 0 was requested */

void zp_setmod(zp_t p);
zp_t zp_to_mont(zp_t a_normal);
zp_t zp_from_mont(zp_t a_mont);

/* --- scalar arithmetic ----------------------------------------------- */
zp_t zp_addv(zp_t a, zp_t b);
zp_t zp_subv(zp_t a, zp_t b);
zp_t zp_mulv(zp_t a, zp_t b);
zp_t zp_invv(zp_t a);
zp_t zp_divv(zp_t a, zp_t b);
zp_t zp_negv(zp_t a);
zp_t zp_from_si(long v);
zp_t zp_from_z(mpz_t z);
zp_t zp_ratio_si(long num, long den);

/* Mutating forms, named to match the fp_/mpq_ API so that the solver
 * sources are a mechanical transliteration. */
#define zp_init(x)        ((x) = 0)
#define zp_clear(x)       ((void)0)
#define zp_set(r,a)       ((r) = (a))
#define zp_set_d(r,d)     ((r) = zp_from_si((long)(d)))
#define zp_set_z(r,z)     ((r) = zp_from_z(z))
#define zp_set_si(r,n,d)  ((r) = zp_ratio_si((long)(n),(long)(d)))
#define zp_set_ui(r,n,d)  ((r) = zp_ratio_si((long)(n),(long)(d)))
#define zp_add(r,a,b)     ((r) = zp_addv((a),(b)))
#define zp_sub(r,a,b)     ((r) = zp_subv((a),(b)))
#define zp_mul(r,a,b)     ((r) = zp_mulv((a),(b)))
#define zp_div(r,a,b)     ((r) = zp_divv((a),(b)))
#define zp_neg(r,a)       ((r) = zp_negv(a))
#define zp_sgn(a)         ((a) != 0)

/* --- containers ------------------------------------------------------ */
zp_vec_t zp_vec(int rows, int valuenum, int valueden);
zp_mat_t zp_matzeros(int rows, int cols);
void     zp_matfree(zp_mat_t m, int rows);
void     zp_vecfree(zp_vec_t v);
zp_msp_t zp_msp(int rows, int cols, int numnnz);
void     zp_mspfree(zp_msp_t msp);
zp_t*    zp_mspget(zp_msp_t msp, int row, int col);
void     zp_matset_si(zp_mat_t mat, int row, int col, int coeffnum, int coeffden);
void     zp_mspset_si(zp_msp_t msp, int row, int col, int coeffnum, int coeffden);
void     zp_mspset_z(zp_msp_t msp, int row, int col, mpz_t coeffnum, int coeffden);
void     zp_mspvecmul(zp_vec_t rop, zp_msp_t op1, zp_vec_t op2);
void     zp_vecdup(zp_vec_t v1, zp_vec_t v2, int n);
void     zp_matprint(zp_mat_t mat, int n, int m);

/* --- dense LU over Z/p ------------------------------------------------
 * Crout elimination with the same row-swap bookkeeping as the rational
 * and floating-point versions.  Pivoting is by first nonzero rather than
 * by magnitude: over a finite field there is no magnitude, and no growth
 * to control, so any nonzero pivot is as good as any other.  A block
 * with no nonzero pivot is singular mod p; that is either a genuine
 * singularity or an unlucky prime, and the caller distinguishes the two
 * by trying further primes. */
int* zp_ludcmp(zp_mat_t A, int N);
int  zp_lubksb(zp_mat_t A, zp_vec_t b, int N, int* indx);

/* --- reconstruction --------------------------------------------------- */

/* Incremental CRT (Garner).  acc holds a residue modulo M, r is the
 * residue of the same quantity modulo the new prime p.  Updates acc in
 * place to a residue modulo M*p; the caller advances M once per prime,
 * after all residues of that prime have been folded in. */
void zp_crt_step(mpz_t acc, const mpz_t M, zp_t r, zp_t p);

/* Rational reconstruction (Wang).  Recovers num/den with
 * num = c*den (mod M) and 2*|num|*den < M.  Returns 1 on success. */
int zp_ratrecon(mpz_t num, mpz_t den, const mpz_t c, const mpz_t M);

/* Rational reconstruction with an explicit denominator bound.  Plain
 * Wang splits the modulus evenly between numerator and denominator, so it
 * needs 2B bits to recover a B-bit value even when the denominator is
 * tiny.  The normalizing constants exported here have denominators of at
 * most a few dozen bits (mdecrease divides only by populations, station
 * multiplicities and think times), so bounding the denominator recovers
 * them from about B bits, halving the number of primes required.
 * Uniqueness still holds: |num| <= M/(2*denbound) and den <= denbound
 * together give 2*|num|*den <= M.  Returns 1 on success. */
int zp_ratrecon_bounded(mpz_t num, mpz_t den, const mpz_t c, const mpz_t M,
                        const mpz_t denbound);

/* Symmetric lift: the representative of c in (-M/2, M/2]. */
void zp_symlift(mpz_t out, const mpz_t c, const mpz_t M);

#endif
