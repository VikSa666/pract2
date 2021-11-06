#include "methods.h"

#define T 2*PI

void main_process(double **, double, int, double, double, int, int);

int main(void) {
    // Declare the matrix which will be used, the dimension and some other variables
    int dimension, option, time_discretization;
    double **u, tol, h, max_t;
    return 0;
}

void main_process(double **A, double tol, int dim, double h, double dt, int N, int option) {
    // Our iterator and auxiliar variable to print out results
    int time_iter, it;
    // This is to iterate over time
    double dt_prev, dt_post;
    // To print the error
    double error;
    // Iterate over the time to take different times in [0,2π] in N steps
    for(time_iter = 0; time_iter < N; time_iter++) {
        // Take the t_n and t_{n+1}
        dt_prev = dt*(double)time_iter;
        dt_post = dt*((double)time_iter + 1);
        // Do the algorithm
        it = conjugated_gradient(A, tol, dim, h, dt_prev, dt, option);
        // Calculate the error
        // TODO
        // Print everything
        printf("Time = %.5le\tIterations = %d\tError = %.12le\n", dt_prev, it, error);
    }
}