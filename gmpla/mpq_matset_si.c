#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "gmpla.h"

/* mpq_set_si takes an unsigned denominator: a negative coeffden became a
 * huge value and the result was not canonical. */
static void mpq_set_si_checked(mpq_t q, long num, long den)
{
	if (den == 0) { fprintf(stderr, "gmpla: zero denominator\n"); exit(1); }
	if (den < 0) { num = -num; den = -den; }
	mpq_set_si(q, num, (unsigned long) den);
	mpq_canonicalize(q);
}

void mpq_matset_si(mpq_mat_t mat, int row, int col, int coeffnum, int coeffden)
{
	mpq_set_si_checked(mat[row][col],coeffnum,coeffden);
}
