#ifndef ROUTING2VISITS_H
#define ROUTING2VISITS_H

#include <gmp.h>

typedef struct {
    mpq_t **data;
    int rows;
    int cols;
} mpq_matrix_t;

typedef struct {
    mpq_matrix_t *matrices;
    int num_classes;
} routing_matrices_t;

mpq_matrix_t* mpq_matrix_init(int rows, int cols);
void mpq_matrix_free(mpq_matrix_t *matrix);
void mpq_matrix_set(mpq_matrix_t *matrix, int row, int col, mpq_t value);
void mpq_matrix_get(mpq_t result, mpq_matrix_t *matrix, int row, int col);
void mpq_matrix_print(mpq_matrix_t *matrix);

routing_matrices_t* routing_matrices_init(int num_classes, int num_stations);
void routing_matrices_free(routing_matrices_t *routing);

mpq_matrix_t* routing2visits(routing_matrices_t *routing);

int mpq_matrix_solve(mpq_matrix_t *A, mpq_matrix_t *b, mpq_matrix_t *x);

#endif