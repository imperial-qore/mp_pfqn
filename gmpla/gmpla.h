#ifndef  __GMP_H__
#include <gmp.h>
#endif
#ifndef MPQ_MAT
#define MPQ_MAT
#include <stdio.h>
#include <stdlib.h>

typedef mpq_t* mpq_vec_t;
typedef mpq_t** mpq_mat_t;
typedef struct struct_mpq_msp
{
	mpq_vec_t coeff;
	int* pos_row;
	int* pos_col;
	int rows;
	int cols;
	int nnz;
	int lastnnz;
} *mpq_msp_t;
extern mpq_t rmul;

/* operations on integer data structures */
int** int_mat(int rows, int cols, int initval);
long int int_matmatchrow(int** mat, long int nrows, long int ncols,int* row);
void int_matprint(int** mat, int n, int m);
int* int_vec(int n, int value);
void int_vecprint(int* vec, int n);
int int_vecsubsum(int* v, int from, int to);
int* int_vecadd(int* v1, int* v2, int n);
int* int_vecdiff(int* v1, int* v2, int n);



/* operations on multiprecision data structures */
int mpq_lubksb(mpq_mat_t A, mpq_vec_t b,int N,int* indx);
int* mpq_ludcmp(mpq_mat_t A, int N);
void mpq_matprint(mpq_mat_t mat, int n, int m);
mpq_t** mpq_matzeros(int rows, int cols);
mpq_msp_t mpq_msp(int rows, int cols, int numnnz);
mpq_t* mpq_mspget(mpq_msp_t msp, int row, int col);
void mpq_matset_si(mpq_mat_t mat, int row, int col, int coeffnum, int coeffden);
void mpq_mspset_si(mpq_msp_t msp, int row, int col, int coeffnum, int coeffden);
void mpq_mspset_z(mpq_msp_t msp, int row, int col, mpz_t coeffnum, int coeffden);
void mpq_mspvecmul(mpq_vec_t rop, mpq_msp_t op1, mpq_vec_t op2);
mpq_vec_t mpq_vec(int rows,int valuenum, int valueden);
void mpq_vecdup(mpq_vec_t v1, mpq_vec_t v2, int n);
void mpq_vecprint(mpq_vec_t, int n);
void mpq_mspprint(mpq_msp_t msp);
void mpq_matvecmul(mpq_vec_t rop, mpq_mat_t A, mpq_vec_t v, int N);
void mpq_mattransmul(mpq_mat_t rop, mpq_mat_t A, mpq_mat_t B, int nrows, int ncols);

#endif
