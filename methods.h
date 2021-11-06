#include "auxFunctions.h"

#define T 2*PI

#include "auxFunctions.h"

double boundary_values(double **A, double t, double h, int dim, int option) {
    int i;
    for(i = 0; i < dim; i++) {
        A[0][i] = g(t, 0, h*i, option);
        A[i][0] = g(t, h*i, 0, option);
        A[dim-1][i] = g(t, h*(dim+1), h*i, option);
        A[i][dim-1] = g(t, h*i, h*(dim+1), option);
    }
}

double crank_nicloson(double **u, double h, double dt, int dim, double tol, int option) {
    int i, j;
    double alpha = dt/(2*h*h);
    for(i = 0; i < dim; i++)
        for(j = 0; j < dim; j++)
            u[i][j] += 2*(1-2*alpha)*u[i][j] + alpha*(u[i+1][j]+u[i-1][j]+u[i][j+1]+u[i][j-1]);
}

double conjugated_gradient(double **A, double tol, int dim, double h, double t, double dt, int option) {
    // Iterator helpers
    int i, j, it_num = 0;
    // Auxiliar matrices
    double **p, **q, **r;
    // Auxiliar doubles
    double alpha, beta, pre_norm, post_norm;
    // Store the matrices
    p = store_matrix(dim+2, dim+2);
    q = store_matrix(dim+2, dim+2);
    r = store_matrix(dim+2, dim+2);
    // Store A·r in r
    special_product(r, A, dim+2);
    // Iteration number 0: we create the matrix...
    // First the boundary values
    boundary_values(A, t, h, dim, option);
    // And then the Crank-Nicolson iterations
    crank_nicloson(A, h, dt, dim, tol, option);
    // Iterations of CG method
    do {
        // Now q = A·p and p remains p
        special_product(q, p, dim+2);
        // I define this to save time on implementing the function
        pre_norm = vector_product_matrix(r, r, dim+2);
        // p^T·A·p in denominator
        alpha = (-1.*pre_norm)/vector_product_matrix(p, q, dim+2);
        // Update values of A and r in the iteration k+1
        for(i = 1; i < dim+1; i++) {
            for(j = 1; j < dim+1; j++) {
                A[i][j] -= alpha*p[i][j];
                r[i][j] += alpha*q[i][j];
            }
        }
        // Now we calculate norm_post, after the update, to calculate beta
        post_norm = vector_product_matrix(r, r, dim+2);
        // Now we calculate beta
        beta = post_norm/pre_norm;
        // With beta already calculated, we update the value of p in the iteration k+1
        for(i = 0; i < dim+1; i++)
            for(j = 0; j < dim+1; j++)
                p[i][j] = r[i][j] + beta;
        // Acumulate the number of iterations
        it_num++;
    // Finishing condition of the method itself
    } while(post_norm >= tol*tol);
    // Free all the memory
    free_matrix(p, dim+2);
    free_matrix(q, dim+2);
    free_matrix(r, dim+2);
    return it_num;
}