#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 4*atan(1)

double f(double t, double x, double y, int option) {
    return (cos(t)+8*PI*PI*sin(t))*sin(2*PI*x)*cos(2*PI*y);
}

double g(double t, double x, double y, int option) {
    return sin(2*PI*x)*cos(2*PI*y)*sin(t);
}