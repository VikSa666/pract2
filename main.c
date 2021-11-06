#include "methods.h"


void main_process(double**, double**, double, int, double, double, int, int);

int main(void) {
    // Declare the matrix which will be used, the dimension and some other variables
    int dimension, option, time_discretization;
    double **u, **v, tol, h, delta_t, max_time;

    // Ask for some constants by terminal to the user
    printf("Choose exercise 1 or 2: ");
    scanf("%d", &option);
    if(option != 1 && option != 2) {
        printf("Those are not valid exercises...exiting program.\n");
        exit(1);
    }
    printf("Dimension of grid: ");
    scanf("%d", &dimension);
    printf("Tolerance: ");
    scanf("%le", &tol);
    printf("Time discretization: ");
    scanf("%d", &time_discretization);
    printf("\n");

    // Initialize the variables
    u = store_matrix(dimension+2, dimension+2);
    v = store_matrix(dimension+2, dimension+2);
    h = 1./(dimension+1);
    switch (option) {
    case 1:
        max_time = 2*PI;
        break;
    
    case 2:
        max_time = 10.;
        break;
    default:
        break;
    }
    delta_t = max_time/time_discretization;
    // Execute the main process
    main_process(u, v, tol, dimension+2, h, delta_t, time_discretization, option);
    // Free memory and finish
    free_matrix(u, dimension+2);
    free_matrix(v, dimension+2);
    return 0;
}

/***
 * This function implements the main algorithm and the aim of it is to iterate over a discretizated
 * interval of the shape [0,T]. If we choose N points on the interval, then the algorithm of Crank
 * Nicolson will be executed N times, so it is a very great deal of operations and processes
 **/
void main_process(double** A, double** B, double tol, int dim, double h, double dt, int N, int option) {
    // Our iterator and auxiliar variable to print out results
    int time_iter, it;
    // This is to iterate over time
    double dt_prev, dt_post;
    // To print the error
    double error;
    // Iterate over the time to take different times in [0,2π] in N steps
    for(time_iter = 0; time_iter < N+1; time_iter++) {
        // Re-initialize the matrix
        fill_zeroes(dim, A);
        // Take the t_n and t_{n+1}
        dt_prev = dt*(double)time_iter;
        dt_post = dt*((double)time_iter + 1);
        // Do the algorithm
        it = conjugated_gradient(A, B, tol, dim, h, dt_prev, dt, option);
        // Calculate the error
        error = error_method(A, h, dt_prev, dim, option);
        // Print everything
        printf("Time = %.5le\tIterations = %d\tError = %.12le\n", dt_prev, it, error);
    }
}