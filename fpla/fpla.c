#include "fpla.h"

mpfr_prec_t FP_WPREC = 53;
mpfr_prec_t FP_RPREC = 106;

void fp_set_si(fp_t x, long num, long den)
{
	mpfr_set_si(x, num, MPFR_RNDN);
	if (den != 1) mpfr_div_si(x, x, den, MPFR_RNDN);
}

void fp_set_ui(fp_t x, unsigned long num, unsigned long den)
{
	mpfr_set_ui(x, num, MPFR_RNDN);
	if (den != 1) mpfr_div_ui(x, x, den, MPFR_RNDN);
}

fp_vec_t fp_vec(int rows, int valuenum, int valueden)
{
	int i;
	mpfr_t* vec = (mpfr_t*) calloc(rows+1, sizeof(mpfr_t));
	for (i=0;i<rows+1;i++)
	{
		mpfr_init2(vec[i], FP_WPREC);
		fp_set_si(vec[i], valuenum, valueden);
	}
	return (fp_vec_t) vec;
}

fp_mat_t fp_matzeros(int rows, int cols)
{
	int i,j;
	fp_mat_t matrix = (mpfr_t**) calloc(rows+1, sizeof(mpfr_t*));
	for (i=0;i<rows+1;i++)
		matrix[i] = (mpfr_t*) calloc(cols+1, sizeof(mpfr_t));
	for (i=0;i<rows+1;i++)
		for (j=0;j<cols+1;j++)
		{
			mpfr_init2(matrix[i][j], FP_WPREC);
			mpfr_set_zero(matrix[i][j], 1);
		}
	return matrix;
}

void fp_matfree(fp_mat_t m, int rows, int cols)
{
	int i,j;
	if (!m) return;
	for (i=0;i<rows+1;i++)
	{
		if (!m[i]) continue;
		for (j=0;j<cols+1;j++) mpfr_clear(m[i][j]);
		free(m[i]);
	}
	free(m);
}

fp_msp_t fp_msp(int rows, int cols, int numnnz)
{
	fp_msp_t msp = (fp_msp_t) calloc(1, sizeof(struct struct_fp_msp));
	msp->coeff = (fp_vec_t) fp_vec(numnnz,0,1);
	msp->pos_row = (int*) calloc(numnnz, sizeof(int));
	msp->pos_col = (int*) calloc(numnnz, sizeof(int));
	msp->nnz = numnnz;
	msp->lastnnz = -1;
	msp->rows = rows;
	msp->cols = cols;
	return msp;
}

fp_t* fp_mspget(fp_msp_t msp, int row, int col)
{
	int k;
	for (k=0;k<=msp->lastnnz;k++)
		if (msp->pos_row[k]==row && msp->pos_col[k]==col) return &msp->coeff[k];
	return NULL;
}

void fp_matset_si(fp_mat_t mat, int row, int col, int coeffnum, int coeffden)
{
	fp_set_si(mat[row][col], coeffnum, coeffden);
}

void fp_mspset_si(fp_msp_t msp, int row, int col, int coeffnum, int coeffden)
{
	fp_t* elem = fp_mspget(msp,row,col);
	if (elem == NULL)
	{
		msp->lastnnz++;
		msp->pos_row[msp->lastnnz]=row;
		msp->pos_col[msp->lastnnz]=col;
		fp_set_si(msp->coeff[msp->lastnnz], coeffnum, coeffden);
	}
	else
		fp_set_si(*elem, coeffnum, coeffden);
}

void fp_mspset_z(fp_msp_t msp, int row, int col, mpz_t coeffnum, int coeffden)
{
	fp_t* elem = fp_mspget(msp,row,col);
	if (elem == NULL)
	{
		msp->lastnnz++;
		msp->pos_row[msp->lastnnz]=row;
		msp->pos_col[msp->lastnnz]=col;
		mpfr_set_z(msp->coeff[msp->lastnnz], coeffnum, MPFR_RNDN);
		if (coeffden != 1)
			mpfr_div_ui(msp->coeff[msp->lastnnz], msp->coeff[msp->lastnnz],
			            (unsigned long)coeffden, MPFR_RNDN);
	}
	else
	{
		mpfr_set_z(*elem, coeffnum, MPFR_RNDN);
		if (coeffden != 1)
			mpfr_div_ui(*elem, *elem, (unsigned long)coeffden, MPFR_RNDN);
	}
}

void fp_mspvecmul(fp_vec_t rop, fp_msp_t op1, fp_vec_t op2)
{
	int t,i;
	mpfr_t rmul; mpfr_init2(rmul, FP_WPREC);
	for (i=0;i<op1->rows;i++)
		mpfr_set_zero(rop[i], 1);
	for (t=0;t<op1->nnz;t++)
	{
		fp_mul(rmul, op1->coeff[t], op2[op1->pos_col[t]]);
		fp_add(rop[op1->pos_row[t]], rop[op1->pos_row[t]], rmul);
	}
	mpfr_clear(rmul);
}

void fp_vecdup(fp_vec_t v1, fp_vec_t v2, int n)
{
	int t;
	for (t=0;t<n;t++) fp_set(v1[t], v2[t]);
}

void fp_vecprint(fp_vec_t vec, int n)
{
	int i;
	for (i=1;i<=n;i++) printf("%.15e\n", fp_get_d(vec[i-1]));
	printf("\n");
}

void fp_matprint(fp_mat_t mat, int n, int m)
{
	int i,j;
	for (i=1;i<=n;i++)
	{
		for (j=1;j<=m;j++) printf("%3g ", fp_get_d(mat[i-1][j-1]));
		printf("\n");
	}
	printf("\n");
}

void fp_mspprint(fp_msp_t msp)
{
	int i,j;
	fp_t* t;
	printf("max nnz=%d ", msp->nnz);
	printf("current nnz=%d\n", msp->lastnnz+1);
	for (i=0;i<msp->rows;i++)
	{
		for (j=0;j<msp->cols;j++)
		{
			t = fp_mspget(msp,i,j);
			if (t != NULL) printf("%3g ", fp_get_d(*t));
			else printf("  0 ");
		}
		printf("\n");
	}
	printf("\n");
}

/* Crout LU with implicit-scaling partial pivoting, identical in structure
 * to gmpla/mpq_ludcmp.c so that the two solvers differ only in the
 * arithmetic they use. */
int* fp_ludcmp(fp_mat_t A, int N)
{
	int i,j,k;
	int imax=0;
	int* indx = (int*) calloc(N+1, sizeof(int));
	mpfr_t max, sum, t, w;
	mpfr_init2(max,FP_WPREC); mpfr_init2(sum,FP_WPREC);
	mpfr_init2(t,FP_WPREC);   mpfr_init2(w,FP_WPREC);
	fp_vec_t vv = fp_vec(N+1,0,1);

	fp_set_d(t,1);
	for (i=1;i<=N;i++)
	{
		fp_set_d(max,0);
		for (j=1;j<=N;j++)
		{
			fp_abs(w, A[i-1][j-1]);
			if (fp_cmp(w,max)>0) fp_set(max,w);
		}
		if (fp_sgn(max)==0) /* structurally singular row */
		{
			for (i=0;i<N+1;i++) mpfr_clear(vv[i]);
			free(vv); free(indx);
			mpfr_clear(max); mpfr_clear(sum); mpfr_clear(t); mpfr_clear(w);
			return NULL;
		}
		fp_div(vv[i-1], t, max);
	}
	for (j=1;j<=N;j++)
	{
		for (i=1;i<=j-1;i++)
		{
			fp_set(sum, A[i-1][j-1]);
			for (k=1;k<=i-1;k++)
			{
				if (fp_sgn(A[i-1][k-1])!=0 && fp_sgn(A[k-1][j-1])!=0)
				{
					fp_mul(t, A[i-1][k-1], A[k-1][j-1]);
					fp_sub(sum, sum, t);
				}
			}
			fp_set(A[i-1][j-1], sum);
		}
		fp_set_d(max,0);
		imax=j;
		for (i=j;i<=N;i++)
		{
			fp_set(sum, A[i-1][j-1]);
			for (k=1;k<=j-1;k++)
			{
				if (fp_sgn(A[i-1][k-1])!=0 && fp_sgn(A[k-1][j-1])!=0)
				{
					fp_mul(t, A[i-1][k-1], A[k-1][j-1]);
					fp_sub(sum, sum, t);
				}
			}
			fp_set(A[i-1][j-1], sum);
			fp_abs(w, sum);
			fp_mul(w, w, vv[i-1]);
			if (fp_cmp(w,max)>=0) { fp_set(max,w); imax=i; }
		}
		if (j != imax)
		{
			for (k=1;k<=N;k++)
			{
				fp_set(w, A[imax-1][k-1]);
				fp_set(A[imax-1][k-1], A[j-1][k-1]);
				fp_set(A[j-1][k-1], w);
			}
			fp_set(vv[imax-1], vv[j-1]);
		}
		indx[j-1]=imax;
		if (j != N)
		{
			fp_set_d(t,1.0);
			if (fp_sgn(A[j-1][j-1])==0)
			{
				for (i=0;i<N+1;i++) mpfr_clear(vv[i]);
				free(vv); free(indx);
				mpfr_clear(max); mpfr_clear(sum); mpfr_clear(t); mpfr_clear(w);
				return NULL;
			}
			fp_div(w, t, A[j-1][j-1]);
			for (i=j+1;i<=N;i++)
				if (fp_sgn(A[i-1][j-1])!=0 && fp_sgn(w)!=0)
					fp_mul(A[i-1][j-1], A[i-1][j-1], w);
		}
	}
	if (fp_sgn(A[N-1][N-1])==0)
	{
		for (i=0;i<N+1;i++) mpfr_clear(vv[i]);
		free(vv); free(indx);
		mpfr_clear(max); mpfr_clear(sum); mpfr_clear(t); mpfr_clear(w);
		return NULL;
	}
	mpfr_clear(max); mpfr_clear(sum); mpfr_clear(t); mpfr_clear(w);
	for (i=0;i<N+1;i++) mpfr_clear(vv[i]);
	free(vv);
	return indx;
}

int fp_lubksb(fp_mat_t A, fp_vec_t b, int N, int* indx)
{
	int i,j;
	int l=0;
	mpfr_t t, sum;
	mpfr_init2(t,FP_WPREC); mpfr_init2(sum,FP_WPREC);
	/* Ly = b */
	for (i=1;i<=N;i++)
	{
		int p = indx[i-1];
		fp_set(sum, b[p-1]);
		fp_set(b[p-1], b[i-1]);
		if (l != 0)
		{
			for (j=l;j<=i-1;j++)
				if (fp_sgn(A[i-1][j-1])!=0)
				{
					fp_mul(t, A[i-1][j-1], b[j-1]);
					fp_sub(sum, sum, t);
				}
		}
		else if (fp_sgn(sum) != 0) l = i;
		fp_set(b[i-1], sum);
	}
	/* Ux = y */
	for (i=N;i>=1;i--)
	{
		fp_set(sum, b[i-1]);
		for (j=i+1;j<=N;j++)
			if (fp_sgn(A[i-1][j-1])!=0)
			{
				fp_mul(t, A[i-1][j-1], b[j-1]);
				fp_sub(sum, sum, t);
			}
		fp_set(b[i-1], sum);
		if (fp_sgn(A[i-1][i-1])==0)
		{
			mpfr_clear(t); mpfr_clear(sum);
			return -1;
		}
		fp_div(b[i-1], b[i-1], A[i-1][i-1]);
	}
	mpfr_clear(t); mpfr_clear(sum);
	return 0;
}

fp_mat_t fp_matdup(fp_mat_t A, int N)
{
	int i,j;
	fp_mat_t B = fp_matzeros(N,N);
	for (i=0;i<N;i++)
		for (j=0;j<N;j++)
			fp_set(B[i][j], A[i][j]);
	return B;
}

void fp_residual(fp_vec_t r, fp_mat_t A, fp_vec_t x, fp_vec_t b, int N)
{
	int i,j;
	mpfr_t acc, prod, xi;
	mpfr_init2(acc, FP_RPREC);
	mpfr_init2(prod, FP_RPREC);
	mpfr_init2(xi, FP_RPREC);
	for (i=0;i<N;i++)
	{
		mpfr_set(acc, b[i], MPFR_RNDN);        /* exact: FP_RPREC >= FP_WPREC */
		for (j=0;j<N;j++)
		{
			if (mpfr_sgn(A[i][j])==0) continue;
			mpfr_set(xi, x[j], MPFR_RNDN);
			mpfr_mul(prod, A[i][j], xi, MPFR_RNDN);
			mpfr_sub(acc, acc, prod, MPFR_RNDN);
		}
		mpfr_set(r[i], acc, MPFR_RNDN);        /* round back to working precision */
	}
	mpfr_clear(acc); mpfr_clear(prod); mpfr_clear(xi);
}

int fp_lurefine(fp_mat_t Aorig, fp_mat_t LU, int* indx, fp_vec_t b,
                fp_vec_t x, int N, int maxit, double tol, double* lastrel)
{
	int i, it;
	fp_vec_t r = fp_vec(N,0,1);
	mpfr_t nd, nx, w;
	mpfr_init2(nd,FP_WPREC); mpfr_init2(nx,FP_WPREC); mpfr_init2(w,FP_WPREC);

	/* initial solve: fp_lubksb overwrites its right-hand side */
	fp_vecdup(x, b, N);
	if (fp_lubksb(LU, x, N, indx) < 0)
	{
		for (i=0;i<N+1;i++) mpfr_clear(r[i]);
		free(r);
		mpfr_clear(nd); mpfr_clear(nx); mpfr_clear(w);
		return -1;
	}
	if (lastrel) *lastrel = 0.0;

	for (it=1; it<=maxit; it++)
	{
		fp_residual(r, Aorig, x, b, N);
		if (fp_lubksb(LU, r, N, indx) < 0) break;
		/* x <- x + d, and measure ||d||_inf / ||x||_inf */
		mpfr_set_zero(nd,1); mpfr_set_zero(nx,1);
		for (i=0;i<N;i++)
		{
			fp_abs(w, r[i]); if (fp_cmp(w,nd)>0) fp_set(nd,w);
			fp_add(x[i], x[i], r[i]);
			fp_abs(w, x[i]); if (fp_cmp(w,nx)>0) fp_set(nx,w);
		}
		double rel = (fp_sgn(nx)==0) ? 0.0 : fp_get_d(nd)/fp_get_d(nx);
		if (lastrel) *lastrel = rel;
		if (rel <= tol) { it++; break; }
	}

	for (i=0;i<N+1;i++) mpfr_clear(r[i]);
	free(r);
	mpfr_clear(nd); mpfr_clear(nx); mpfr_clear(w);
	return it-1;
}

/* Infinity-norm condition number of a matrix whose LU factors are known.
 * ||A||_inf is taken from the pristine matrix; ||A^{-1}||_inf is obtained
 * by solving A x = e_j for every unit vector, which costs O(N^3) but is
 * only ever used for diagnostics. */
double fp_cond_inf(fp_mat_t Aorig, fp_mat_t LU, int* indx, int N)
{
	int i,j;
	double na=0.0, ni=0.0, rowsum;
	fp_vec_t e = fp_vec(N,0,1);
	double* colabs = (double*) calloc(N, sizeof(double));

	for (i=0;i<N;i++)
	{
		rowsum = 0.0;
		for (j=0;j<N;j++) rowsum += fabs(fp_get_d(Aorig[i][j]));
		if (rowsum > na) na = rowsum;
	}
	for (j=0;j<N;j++)
	{
		for (i=0;i<N;i++) mpfr_set_zero(e[i],1);
		fp_set_si(e[j],1,1);
		if (fp_lubksb(LU,e,N,indx) < 0) { free(colabs); for(i=0;i<N+1;i++) mpfr_clear(e[i]); free(e); return -1.0; }
		for (i=0;i<N;i++) colabs[i] += fabs(fp_get_d(e[i]));
	}
	for (i=0;i<N;i++) if (colabs[i] > ni) ni = colabs[i];

	free(colabs);
	for (i=0;i<N+1;i++) mpfr_clear(e[i]);
	free(e);
	return na*ni;
}
