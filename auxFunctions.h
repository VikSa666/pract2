#include "functions.h"


double **store_matrix(int rows, int cols) {
    double **A;
    int i;
    A = (double**)calloc(rows, sizeof(double*));
    if(A == NULL) {
        printf("Error storing matrix.\n");
        exit(1);
    }
    for(i = 0; i < rows; i++) {
        A[i] = (double*)calloc(cols, sizeof(double));
        if(A[i] == NULL) {
            printf("Error storing matrix.\n");
            exit(1);
        }
    }
    return A;
}

double free_matrix(double **A, int cols) {
    int i;
    for(i = 0; i < cols; i++) free(A[i]);
    free(A);
}

double special_product(double **x, double **mat, int dim) {
    int i, j;
    for(i = 0; i < dim; i++)
        for(j = 0; j < dim; j++)
            x[i][j] = 4*mat[i][j]-mat[i-1][j]-mat[i+1][j]-mat[i][j-1]-mat[i][j+1];
}

double sup_norm(double **x, double dim) {
    int i, j;
    double norm = 0.;
    for(i = 0; i < dim; i++)
        for(j = 0; j < dim; j++)
            if(fabs(x[i][j]) < norm)
                norm = fabs(x[i][j]);
    return norm;
}

double matrix_sum(double **x, double **y, int dim) {
    int i, j;
    for(i = 0; i < dim; i++)
        for(j = 0; j < dim; j++)
            x[i][j] += y[i][j];
}

double vector_product_matrix(double **x, double **y, int dim) {
    int i, j;
    double prod = 0.;
    for(i = 0; i < dim; i++)
        for(j = 0; j < dim; j++)
            prod += x[i][j]*y[i][j];
    return prod;
}

