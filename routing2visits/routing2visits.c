#include "routing2visits.h"
#include <stdio.h>
#include <stdlib.h>

mpq_matrix_t* mpq_matrix_init(int rows, int cols) {
    mpq_matrix_t *matrix = malloc(sizeof(mpq_matrix_t));
    matrix->rows = rows;
    matrix->cols = cols;
    matrix->data = malloc(rows * sizeof(mpq_t*));
    
    for (int i = 0; i < rows; i++) {
        matrix->data[i] = malloc(cols * sizeof(mpq_t));
        for (int j = 0; j < cols; j++) {
            mpq_init(matrix->data[i][j]);
        }
    }
    
    return matrix;
}

void mpq_matrix_free(mpq_matrix_t *matrix) {
    if (matrix == NULL) return;
    
    for (int i = 0; i < matrix->rows; i++) {
        for (int j = 0; j < matrix->cols; j++) {
            mpq_clear(matrix->data[i][j]);
        }
        free(matrix->data[i]);
    }
    free(matrix->data);
    free(matrix);
}

void mpq_matrix_set(mpq_matrix_t *matrix, int row, int col, mpq_t value) {
    mpq_set(matrix->data[row][col], value);
}

void mpq_matrix_get(mpq_t result, mpq_matrix_t *matrix, int row, int col) {
    mpq_set(result, matrix->data[row][col]);
}

void mpq_matrix_print(mpq_matrix_t *matrix) {
    for (int i = 0; i < matrix->rows; i++) {
        for (int j = 0; j < matrix->cols; j++) {
            gmp_printf("%Qd ", matrix->data[i][j]);
        }
        printf("\n");
    }
}

routing_matrices_t* routing_matrices_init(int num_classes, int num_stations) {
    routing_matrices_t *routing = malloc(sizeof(routing_matrices_t));
    routing->num_classes = num_classes;
    routing->matrices = malloc(num_classes * sizeof(mpq_matrix_t));
    
    for (int c = 0; c < num_classes; c++) {
        routing->matrices[c] = *mpq_matrix_init(num_stations, num_stations);
    }
    
    return routing;
}

void routing_matrices_free(routing_matrices_t *routing) {
    if (routing == NULL) return;
    
    for (int c = 0; c < routing->num_classes; c++) {
        for (int i = 0; i < routing->matrices[c].rows; i++) {
            for (int j = 0; j < routing->matrices[c].cols; j++) {
                mpq_clear(routing->matrices[c].data[i][j]);
            }
            free(routing->matrices[c].data[i]);
        }
        free(routing->matrices[c].data);
    }
    free(routing->matrices);
    free(routing);
}

int mpq_matrix_solve(mpq_matrix_t *A, mpq_matrix_t *b, mpq_matrix_t *x) {
    int n = A->rows;
    
    mpq_matrix_t *aug = mpq_matrix_init(n, n + 1);
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            mpq_set(aug->data[i][j], A->data[i][j]);
        }
        mpq_set(aug->data[i][n], b->data[i][0]);
    }
    
    mpq_t temp, pivot, factor;
    mpq_init(temp);
    mpq_init(pivot);
    mpq_init(factor);
    
    for (int i = 0; i < n; i++) {
        int max_row = i;
        for (int k = i + 1; k < n; k++) {
            mpq_t abs_k, abs_max;
            mpq_init(abs_k);
            mpq_init(abs_max);
            mpq_abs(abs_k, aug->data[k][i]);
            mpq_abs(abs_max, aug->data[max_row][i]);
            if (mpq_cmp(abs_k, abs_max) > 0) {
                max_row = k;
            }
            mpq_clear(abs_k);
            mpq_clear(abs_max);
        }
        
        if (max_row != i) {
            for (int j = 0; j <= n; j++) {
                mpq_swap(aug->data[i][j], aug->data[max_row][j]);
            }
        }
        
        if (mpq_sgn(aug->data[i][i]) == 0) {
            mpq_clear(temp);
            mpq_clear(pivot);
            mpq_clear(factor);
            mpq_matrix_free(aug);
            return -1;
        }
        
        mpq_set(pivot, aug->data[i][i]);
        
        for (int k = i + 1; k < n; k++) {
            mpq_div(factor, aug->data[k][i], pivot);
            
            for (int j = i; j <= n; j++) {
                mpq_mul(temp, factor, aug->data[i][j]);
                mpq_sub(aug->data[k][j], aug->data[k][j], temp);
            }
        }
    }
    
    for (int i = n - 1; i >= 0; i--) {
        mpq_set(x->data[i][0], aug->data[i][n]);
        
        for (int j = i + 1; j < n; j++) {
            mpq_mul(temp, aug->data[i][j], x->data[j][0]);
            mpq_sub(x->data[i][0], x->data[i][0], temp);
        }
        
        mpq_div(x->data[i][0], x->data[i][0], aug->data[i][i]);
    }
    
    mpq_clear(temp);
    mpq_clear(pivot);
    mpq_clear(factor);
    mpq_matrix_free(aug);
    
    return 0;
}

mpq_matrix_t* routing2visits(routing_matrices_t *routing) {
    int num_classes = routing->num_classes;
    int num_stations = routing->matrices[0].rows;
    int num_queues = num_stations - 1;
    
    mpq_matrix_t *visits = mpq_matrix_init(num_classes, num_stations);
    
    mpq_t one;
    mpq_init(one);
    mpq_set_ui(one, 1, 1);
    
    for (int c = 0; c < num_classes; c++) {
        mpq_set(visits->data[c][0], one);
        
        if (num_queues > 0) {
            mpq_matrix_t *A = mpq_matrix_init(num_queues, num_queues);
            mpq_matrix_t *b = mpq_matrix_init(num_queues, 1);
            mpq_matrix_t *queue_visits = mpq_matrix_init(num_queues, 1);
            
            mpq_t temp;
            mpq_init(temp);
            
            for (int i = 0; i < num_queues; i++) {
                mpq_set_ui(A->data[i][i], 1, 1);
                for (int j = 0; j < num_queues; j++) {
                    mpq_sub(A->data[i][j], A->data[i][j], routing->matrices[c].data[i + 1][j + 1]);
                }
                mpq_set(b->data[i][0], routing->matrices[c].data[0][i + 1]);
            }
            
            if (mpq_matrix_solve(A, b, queue_visits) == 0) {
                for (int i = 0; i < num_queues; i++) {
                    mpq_set(visits->data[c][i + 1], queue_visits->data[i][0]);
                }
            } else {
                fprintf(stderr, "Error: Singular matrix for class %d\n", c);
                mpq_matrix_free(visits);
                mpq_matrix_free(A);
                mpq_matrix_free(b);
                mpq_matrix_free(queue_visits);
                mpq_clear(temp);
                mpq_clear(one);
                return NULL;
            }
            
            mpq_clear(temp);
            mpq_matrix_free(A);
            mpq_matrix_free(b);
            mpq_matrix_free(queue_visits);
        }
    }
    
    mpq_clear(one);
    return visits;
}