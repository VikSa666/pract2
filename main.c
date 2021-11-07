#include "methods.h"
#include <time.h>


void main_process(double**, double**, double, int, double, double, int, int, FILE*);

int main(void) {
    // Declare the matrix which will be used, the dimension and some other variables
    int dimension, option, time_discretization;
    double **u, **v, tol, h, delta_t, max_time;
    double begin_t, end_t;
    
    // These will be the files where I will print out the results
    FILE *f1 = fopen("results1.txt", "w"), *f2 = fopen("results2.txt", "w"), *f3 = fopen("results3.txt", "w");
    if(f1 == NULL || f2 == NULL || f3 == NULL) {
        printf("Error opening out file.\n");
        exit(-1);   
    }

    // Ask for some constants by terminal to the user
    printf("Choose exercise 1, 2 or 3: ");
    scanf("%d", &option);
    if(option < 1 || option > 3) {
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
    
    // Execute the main process in a switch to be able to choose exercise
    switch (option) {
    case 1:
        max_time = 2*PI;
        delta_t = max_time/time_discretization;
        begin_t = clock();
        main_process(u, v, tol, dimension+2, h, delta_t, time_discretization, 1, f1);
        end_t = clock();
        printf("Total time of execution: %le", (end_t-begin_t)/CLOCKS_PER_SEC);
        break;
    
    case 2:
        max_time = 10.;
        delta_t = max_time/time_discretization;
        begin_t = clock();
        main_process(u, v, tol, dimension+2, h, delta_t, time_discretization, 2, f2);
        end_t = clock();
        printf("Total time of execution: %le", (end_t-begin_t)/CLOCKS_PER_SEC);
        break;
    case 3:
        max_time = 10.;
        delta_t = max_time/time_discretization;
        begin_t = clock();
        main_process(u, v, tol, dimension+2, h, delta_t, time_discretization, 3, f3);
        end_t = clock();
        printf("Total time of execution: %le\n", (end_t-begin_t)/CLOCKS_PER_SEC);
        break;
    default:
        break;
    }
        
    // Free memory and finish
    free_matrix(u, dimension+2);
    free_matrix(v, dimension+2);
    fclose(f1);
    fclose(f2);
    fclose(f3);
    return 0;
}

//
// This function implements the main algorithm and the aim of it is to iterate over a discretizated
// interval of the shape [0,T]. If we choose N points on the interval, then the algorithm of Crank
// Nicolson will be executed N times, so it is a very great deal of operations and processes
//
void main_process(double** A, double** B, double tol, int dim, double h, double dt, int N, int option, FILE *f) {
    // Our iterator and auxiliar variable to print out results
    int time_iter, it;
    // This is to iterate over time
    double dt_prev, dt_post;
    // To print the error
    double error;
    // Print the title of the file
    fprintf(f, "========================================================================================================\n");
    fprintf(f, "Results of exercise %d with grid of %dx%d, time discretization %d and tolerance %le\n", option, dim, dim, N, tol);
    fprintf(f, "========================================================================================================\n");
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
        fprintf(f, "Time = %.5le \tIterations = %d \tError = %.12le\n", dt_prev, it, error);
    }
}