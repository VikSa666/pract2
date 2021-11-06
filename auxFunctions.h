#include "functions.h"

// Function that stores a matrix with given dimensions
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

// Function that frees the memory of a matrix
void free_matrix(double** A, int cols) {
    int i;
    for(i = 0; i < cols; i++) free(A[i]);
    free(A);
}

// Function that implements the "special product" of the method
void special_product(double** mat, double** x, int dim) {
    int i, j;
    for(i = 1; i < dim-1; i++)
        for(j = 1; j < dim-1; j++)
            x[i][j] = 4*mat[i][j] - mat[i-1][j] - mat[i+1][j] - mat[i][j-1] - mat[i][j+1];
}

// Function that implements the ∞-norm of a matrix
double sup_norm(double** x, double dim) {
    int i, j;
    double norm = 0.;
    for(i = 0; i < dim; i++)
        for(j = 0; j < dim; j++)
            if(fabs(x[i][j]) < norm)
                norm = fabs(x[i][j]);
    return norm;
}

// Function that calculates the dot product of two vectors given in matrix-shape
double vector_product_matrix(double **x, double **y, int dim) {
    int i, j;
    double prod = 0.;
    for(i = 1; i < dim-1; i++)
        for(j = 1; j < dim-1; j++)
            prod += x[i][j]*y[i][j];
    return prod;
}

// Key function to calculate the error (notice that the result is given by the enunciate)
double error_method(double** A, double h, double t, int dim, int option) {
    int i, j;
    double error = 0.;
    for(i = 0; i < dim; i++)
        for(j = 0; j < dim; j++)
            if(error < fabs(A[i][j] - g(t, i*h, j*h, option)))
                error = fabs(A[i][j] - g(t, i*h, j*h, option));
    return error;
}

// Function to print out matrices in a fancy way
void print_matrix(double **A, int n, int m) {
    int i, j;
    for(i = 0; i < n; i++) {
        for(j = 0; j < m; j++) {
            printf("%.4e ", A[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

// Function to just fill of zeroes a matrix
void fill_zeroes(int n, double **A) {
    int i, j;
    for(i = 0; i < n; ++i) 
        for(j = 1; j < n-1; j++) 
            A[i][j] = 0;
}