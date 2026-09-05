#include "zpla.h"
#include <string.h>

__thread zp_t ZP_P = 0;
__thread zp_t ZP_PINV = 0;
__thread zp_t ZP_R2 = 0;
__thread int  ZP_SINGULAR = 0;

/* --- plain (non-Montgomery) helpers, used at setup and by the CRT ---- */
static zp_t mulmod_plain(zp_t a, zp_t b, zp_t p)
{
	return (zp_t)(((__uint128_t)a * (__uint128_t)b) % (__uint128_t)p);
}

static zp_t inv_plain(zp_t a, zp_t p)
{
	long long t = 0, newt = 1;
	long long r = (long long) p, newr = (long long)(a % p);
	if (newr == 0) return 0;
	while (newr != 0)
	{
		long long q = r / newr;
		long long tmp = t - q * newt; t = newt; newt = tmp;
		tmp = r - q * newr;          r = newr; newr = tmp;
	}
	if (t < 0) t += (long long) p;
	return (zp_t) t;
}

/* Montgomery reduction of T < p*2^64. */
static inline zp_t redc(__uint128_t T)
{
	zp_t m = (zp_t)((zp_t)T * ZP_PINV);
	zp_t t = (zp_t)((T + (__uint128_t)m * (__uint128_t)ZP_P) >> 64);
	return (t >= ZP_P) ? t - ZP_P : t;
}

void zp_setmod(zp_t p)
{
	zp_t inv = p, R;
	int i;
	ZP_P = p;
	/* Newton iteration for p^{-1} mod 2^64, starting from p^{-1} mod 2^3. */
	for (i=0;i<6;i++) inv *= 2 - p*inv;
	ZP_PINV = (zp_t)(0 - inv);            /* -p^{-1} mod 2^64 */
	R = (zp_t)(((__uint128_t)1 << 64) % (__uint128_t)p);
	ZP_R2 = mulmod_plain(R, R, p);
	ZP_SINGULAR = 0;
}

zp_t zp_to_mont(zp_t a_normal) { return redc((__uint128_t)a_normal * (__uint128_t)ZP_R2); }
zp_t zp_from_mont(zp_t a_mont) { return redc((__uint128_t)a_mont); }

zp_t zp_addv(zp_t a, zp_t b)
{
	zp_t s = a + b;
	return (s >= ZP_P) ? s - ZP_P : s;
}

zp_t zp_subv(zp_t a, zp_t b)
{
	return (a >= b) ? a - b : a + ZP_P - b;
}

zp_t zp_mulv(zp_t a, zp_t b)
{
	return redc((__uint128_t)a * (__uint128_t)b);
}

zp_t zp_negv(zp_t a)
{
	return a ? ZP_P - a : 0;
}

zp_t zp_invv(zp_t a)
{
	zp_t an, ai;
	if (a == 0) { ZP_SINGULAR = 1; return 0; }
	an = zp_from_mont(a);
	ai = inv_plain(an, ZP_P);
	return zp_to_mont(ai);
}

zp_t zp_divv(zp_t a, zp_t b)
{
	return zp_mulv(a, zp_invv(b));
}

zp_t zp_from_si(long v)
{
	long long m = (long long)(v % (long long)ZP_P);
	if (m < 0) m += (long long) ZP_P;
	return zp_to_mont((zp_t) m);
}

zp_t zp_from_z(mpz_t z)
{
	/* mpz_fdiv_ui returns the non-negative remainder, so negative inputs
	 * map correctly without an extra correction. */
	return zp_to_mont((zp_t) mpz_fdiv_ui(z, (unsigned long) ZP_P));
}

zp_t zp_ratio_si(long num, long den)
{
	long long n = (long long)(num % (long long)ZP_P);
	long long d;
	if (n < 0) n += (long long) ZP_P;
	if (den == 1) return zp_to_mont((zp_t) n);
	d = (long long)(den % (long long)ZP_P);
	if (d < 0) d += (long long) ZP_P;
	if (d == 0) { ZP_SINGULAR = 1; return 0; }
	return zp_to_mont(mulmod_plain((zp_t)n, inv_plain((zp_t)d, ZP_P), ZP_P));
}

zp_vec_t zp_vec(int rows, int valuenum, int valueden)
{
	int i;
	zp_t* v = (zp_t*) calloc(rows+1, sizeof(zp_t));
	zp_t val = zp_ratio_si(valuenum, valueden);
	for (i=0;i<rows+1;i++) v[i] = val;
	return v;
}

void zp_vecfree(zp_vec_t v) { free(v); }

zp_mat_t zp_matzeros(int rows, int cols)
{
	int i;
	zp_mat_t m = (zp_t**) calloc(rows+1, sizeof(zp_t*));
	for (i=0;i<rows+1;i++) m[i] = (zp_t*) calloc(cols+1, sizeof(zp_t));
	return m;
}

void zp_matfree(zp_mat_t m, int rows)
{
	int i;
	if (!m) return;
	for (i=0;i<rows+1;i++) free(m[i]);
	free(m);
}

zp_msp_t zp_msp(int rows, int cols, int numnnz)
{
	zp_msp_t msp = (zp_msp_t) calloc(1, sizeof(struct struct_zp_msp));
	msp->coeff = zp_vec(numnnz,0,1);
	msp->pos_row = (int*) calloc(numnnz, sizeof(int));
	msp->pos_col = (int*) calloc(numnnz, sizeof(int));
	msp->nnz = numnnz;
	msp->lastnnz = -1;
	msp->rows = rows;
	msp->cols = cols;
	return msp;
}

void zp_mspfree(zp_msp_t msp)
{
	if (!msp) return;
	free(msp->coeff); free(msp->pos_row); free(msp->pos_col); free(msp);
}

zp_t* zp_mspget(zp_msp_t msp, int row, int col)
{
	int k;
	for (k=0;k<=msp->lastnnz;k++)
		if (msp->pos_row[k]==row && msp->pos_col[k]==col) return &msp->coeff[k];
	return NULL;
}

void zp_matset_si(zp_mat_t mat, int row, int col, int coeffnum, int coeffden)
{
	mat[row][col] = zp_ratio_si(coeffnum, coeffden);
}

void zp_mspset_si(zp_msp_t msp, int row, int col, int coeffnum, int coeffden)
{
	zp_t* elem = zp_mspget(msp,row,col);
	if (elem == NULL)
	{
		msp->lastnnz++;
		msp->pos_row[msp->lastnnz]=row;
		msp->pos_col[msp->lastnnz]=col;
		msp->coeff[msp->lastnnz] = zp_ratio_si(coeffnum, coeffden);
	}
	else
		*elem = zp_ratio_si(coeffnum, coeffden);
}

void zp_mspset_z(zp_msp_t msp, int row, int col, mpz_t coeffnum, int coeffden)
{
	zp_t v = zp_from_z(coeffnum);
	zp_t* elem;
	if (coeffden != 1) v = zp_mulv(v, zp_invv(zp_from_si(coeffden)));
	elem = zp_mspget(msp,row,col);
	if (elem == NULL)
	{
		msp->lastnnz++;
		msp->pos_row[msp->lastnnz]=row;
		msp->pos_col[msp->lastnnz]=col;
		msp->coeff[msp->lastnnz] = v;
	}
	else
		*elem = v;
}

void zp_mspvecmul(zp_vec_t rop, zp_msp_t op1, zp_vec_t op2)
{
	int t,i;
	for (i=0;i<op1->rows;i++) rop[i] = 0;
	for (t=0;t<op1->nnz;t++)
		rop[op1->pos_row[t]] = zp_addv(rop[op1->pos_row[t]],
		                               zp_mulv(op1->coeff[t], op2[op1->pos_col[t]]));
}

void zp_vecdup(zp_vec_t v1, zp_vec_t v2, int n)
{
	memcpy(v1, v2, (size_t)n * sizeof(zp_t));
}

void zp_matprint(zp_mat_t mat, int n, int m)
{
	int i,j;
	for (i=0;i<n;i++)
	{
		for (j=0;j<m;j++) printf("%llu ", mat[i][j]);
		printf("\n");
	}
	printf("\n");
}

int* zp_ludcmp(zp_mat_t A, int N)
{
	int i,j,k,imax;
	int* indx = (int*) calloc(N+1, sizeof(int));
	zp_t sum, w;

	for (j=1;j<=N;j++)
	{
		for (i=1;i<=j-1;i++)
		{
			sum = A[i-1][j-1];
			for (k=1;k<=i-1;k++)
				if (A[i-1][k-1] && A[k-1][j-1])
					sum = zp_subv(sum, zp_mulv(A[i-1][k-1], A[k-1][j-1]));
			A[i-1][j-1] = sum;
		}
		imax = 0;
		for (i=j;i<=N;i++)
		{
			sum = A[i-1][j-1];
			for (k=1;k<=j-1;k++)
				if (A[i-1][k-1] && A[k-1][j-1])
					sum = zp_subv(sum, zp_mulv(A[i-1][k-1], A[k-1][j-1]));
			A[i-1][j-1] = sum;
			if (imax == 0 && sum != 0) imax = i;
		}
		if (imax == 0) { free(indx); return NULL; } /* singular mod p */
		if (j != imax)
			for (k=1;k<=N;k++)
			{
				w = A[imax-1][k-1];
				A[imax-1][k-1] = A[j-1][k-1];
				A[j-1][k-1] = w;
			}
		indx[j-1] = imax;
		if (j != N)
		{
			if (A[j-1][j-1] == 0) { free(indx); return NULL; }
			w = zp_invv(A[j-1][j-1]);
			for (i=j+1;i<=N;i++)
				if (A[i-1][j-1]) A[i-1][j-1] = zp_mulv(A[i-1][j-1], w);
		}
	}
	if (A[N-1][N-1] == 0) { free(indx); return NULL; }
	return indx;
}

int zp_lubksb(zp_mat_t A, zp_vec_t b, int N, int* indx)
{
	int i,j;
	int l=0;
	zp_t sum;
	for (i=1;i<=N;i++)
	{
		int p = indx[i-1];
		sum = b[p-1];
		b[p-1] = b[i-1];
		if (l != 0)
		{
			for (j=l;j<=i-1;j++)
				if (A[i-1][j-1]) sum = zp_subv(sum, zp_mulv(A[i-1][j-1], b[j-1]));
		}
		else if (sum != 0) l = i;
		b[i-1] = sum;
	}
	for (i=N;i>=1;i--)
	{
		sum = b[i-1];
		for (j=i+1;j<=N;j++)
			if (A[i-1][j-1]) sum = zp_subv(sum, zp_mulv(A[i-1][j-1], b[j-1]));
		if (A[i-1][i-1] == 0) return -1;
		b[i-1] = zp_mulv(sum, zp_invv(A[i-1][i-1]));
	}
	return 0;
}

void zp_crt_step(mpz_t acc, const mpz_t M, zp_t r, zp_t p)
{
	mpz_t t, delta;
	zp_t amodp, mmodp, d;

	if (mpz_sgn(M) == 0) { mpz_set_ui(acc, (unsigned long) r); return; }

	mpz_init(t); mpz_init(delta);
	amodp = (zp_t) mpz_fdiv_ui(acc, (unsigned long) p);
	mmodp = (zp_t) mpz_fdiv_ui(M,   (unsigned long) p);
	/* plain arithmetic: this runs modulo the new prime, not the one the
	 * Montgomery context is currently set up for */
	d = mulmod_plain((r % p >= amodp) ? (r % p) - amodp : (r % p) + p - amodp,
	                 inv_plain(mmodp, p), p);
	mpz_set_ui(delta, (unsigned long) d);
	mpz_mul(t, M, delta);
	mpz_add(acc, acc, t);
	mpz_clear(t); mpz_clear(delta);
}

int zp_ratrecon(mpz_t num, mpz_t den, const mpz_t c, const mpz_t M)
{
	/* Half-extended Euclid on (M, c), stopping when the remainder drops
	 * below sqrt(M/2).  The pair (remainder, cofactor) is then the unique
	 * rational with 2*|num|*den < M congruent to c. */
	mpz_t r0,r1,t0,t1,q,tmp,bound,chk;
	int ok = 0;
	mpz_inits(r0,r1,t0,t1,q,tmp,bound,chk,NULL);

	mpz_fdiv_q_2exp(bound, M, 1);
	mpz_sqrt(bound, bound);

	mpz_set(r0, M); mpz_set(r1, c);
	mpz_set_ui(t0, 0); mpz_set_ui(t1, 1);

	while (mpz_cmp(r1, bound) > 0 && mpz_sgn(r1) != 0)
	{
		mpz_fdiv_q(q, r0, r1);
		mpz_mul(tmp, q, r1); mpz_sub(tmp, r0, tmp); mpz_set(r0, r1); mpz_set(r1, tmp);
		mpz_mul(tmp, q, t1); mpz_sub(tmp, t0, tmp); mpz_set(t0, t1); mpz_set(t1, tmp);
	}

	if (mpz_sgn(t1) != 0)
	{
		mpz_set(num, r1);
		mpz_set(den, t1);
		if (mpz_sgn(den) < 0) { mpz_neg(num, num); mpz_neg(den, den); }
		mpz_gcd(tmp, num, den);
		if (mpz_cmp_ui(tmp, 1) == 0)
		{
			/* verify num = c*den (mod M) */
			mpz_mul(chk, c, den);
			mpz_sub(chk, chk, num);
			mpz_mod(chk, chk, M);
			if (mpz_sgn(chk) == 0) ok = 1;
		}
	}

	mpz_clears(r0,r1,t0,t1,q,tmp,bound,chk,NULL);
	return ok;
}

void zp_symlift(mpz_t out, const mpz_t c, const mpz_t M)
{
	mpz_t half;
	mpz_init(half);
	mpz_fdiv_q_2exp(half, M, 1);
	mpz_set(out, c);
	if (mpz_cmp(out, half) > 0) mpz_sub(out, out, M);
	mpz_clear(half);
}

int zp_ratrecon_bounded(mpz_t num, mpz_t den, const mpz_t c, const mpz_t M,
                        const mpz_t denbound)
{
	mpz_t r0,r1,t0,t1,q,tmp,numbound;
	int ok = 0;
	mpz_inits(r0,r1,t0,t1,q,tmp,numbound,NULL);

	/* numbound = floor(M / (2*denbound)) */
	mpz_mul_ui(numbound, denbound, 2);
	mpz_fdiv_q(numbound, M, numbound);
	if (mpz_sgn(numbound) <= 0) { mpz_clears(r0,r1,t0,t1,q,tmp,numbound,NULL); return 0; }

	mpz_set(r0, M); mpz_set(r1, c);
	mpz_set_ui(t0, 0); mpz_set_ui(t1, 1);

	while (mpz_cmp(r1, numbound) > 0 && mpz_sgn(r1) != 0)
	{
		mpz_fdiv_q(q, r0, r1);
		mpz_mul(tmp, q, r1); mpz_sub(tmp, r0, tmp); mpz_set(r0, r1); mpz_set(r1, tmp);
		mpz_mul(tmp, q, t1); mpz_sub(tmp, t0, tmp); mpz_set(t0, t1); mpz_set(t1, tmp);
	}

	if (mpz_sgn(t1) != 0)
	{
		mpz_set(num, r1);
		mpz_set(den, t1);
		if (mpz_sgn(den) < 0) { mpz_neg(num, num); mpz_neg(den, den); }
		if (mpz_cmp(den, denbound) <= 0)
		{
			mpz_gcd(tmp, num, den);
			if (mpz_cmp_ui(tmp, 1) == 0) ok = 1;
		}
	}

	mpz_clears(r0,r1,t0,t1,q,tmp,numbound,NULL);
	return ok;
}
