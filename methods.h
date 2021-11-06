#include "auxFunctions.h"

#define T 2*PI

//
// Boundary values is a function that stores the boundary values given by g function in the 
// "boundaries" of the matrix
//
void boundary_values(double** A, double t, double h, int dim, int option) {
    int i;
    for(i = 0; i < dim; i++) {
        A[0][i] = g(t, 0, h*i, option);
        A[i][0] = g(t, h*i, 0, option);
        A[dim-1][i] = g(t, h*(dim-1), h*i, option);
        A[i][dim-1] = g(t, h*i, h*(dim-1), option);
    }
}

//
// Crank-Nicoloson method to create the system of n equations with n incognites which will be solved
// later using some method we already know
//
void crank_nicloson(double** u, double **r, double h, double dt, double t, int dim, double tol, int option) {
    int i, j;
    double alpha = dt/(h*h);
    for(i = 1; i < dim-1; i++)
        for(j = 1; j < dim-1; j++)
            r[i][j] = (1-2*alpha)*u[i][j] -r[i][j]
                + (alpha/2)*(u[i+1][j]+u[i-1][j]+u[i][j+1]+u[i][j-1])
                + (alpha/2)*(f(t,i*h, j*h, option) + f(t+dt, i*h, j*h, option));
}

//
// Conjugated Gradient method to solve the equations we talked about before
//
double conjugated_gradient(double** A, double**B, double tol, int dim, double h, double t, double dt, int option) {
    // Declaration of variables
    int i, j, it_num = 0;
    double **p, **q, **r;
    double alpha, beta, pre_norm, post_norm;

    // Store the matrices
    p = store_matrix(dim, dim);
    q = store_matrix(dim, dim);
    r = store_matrix(dim, dim);

    
    
    // Iteration number 0: we create the matrix...
    // First the boundary values
    boundary_values(A, t, h, dim, option);
    boundary_values(B, t+dt, h, dim, option);
    // Store A·r in r
    special_product(B, r, dim);
    // And then the Crank-Nicolson iterations, solution stored in r
    crank_nicloson(A, r, h, dt, dim, t, tol, option);
    for(i = 0; i < dim; i++) for(j = 0; j < dim; j++) p[i][j] = r[i][j];
    // This is to save execution time
    post_norm = vector_product_matrix(r, r, dim);
    // Iterations of CG method
    do {
        // Now q = A·p and p remains p
        special_product(p, q, dim);
        // I define this to save time on implementing the function
        pre_norm = post_norm;
        // p^T·A·p in denominator
        alpha = (-1.*pre_norm)/vector_product_matrix(p, q, dim);
        // Update values of A and r in the iteration k+1
        for(i = 1; i < dim-1; i++) {
            for(j = 1; j < dim-1; j++) {
                B[i][j] -= alpha*p[i][j];
                r[i][j] = alpha*q[i][j] + r[i][j];
            }
        }
        // Now we calculate norm_post, after the update, to calculate beta
        post_norm = vector_product_matrix(r, r, dim);
        // Now we calculate beta
        beta = post_norm/pre_norm;
        // With beta already calculated, we update the value of p in the iteration k+1
        for(i = 1; i < dim-1; i++)
            for(j = 1; j < dim-1; j++)
                p[i][j] = r[i][j] + beta * p[i][j];
        // Acumulate the number of iterations
        it_num++;
        // Finishing condition of the method itself
    } while(post_norm >= tol*tol);

    // Free all the memory
    free_matrix(p, dim);
    free_matrix(q, dim);
    free_matrix(r, dim);
    return it_num;
}