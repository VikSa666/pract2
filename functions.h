#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 4*atan(1)

double f(double t, double x, double y, int option) {
    switch (option) {
    case 1:
        return (cos(t)+8*PI*PI*sin(t))*sin(2*PI*x)*cos(2*PI*y);
        break;
    
    case 2:
        return 0;
        break;
    case 3:
        return 0;
        break;
    }
    return 0;
}

double g(double t, double x, double y, int option) {
    switch (option) {
    case 1:
        return sin(2*PI*x)*cos(2*PI*y)*sin(t);
        break;
    
    case 2:
        return 1;
        break;
    case 3:
        if(t <= 5) return 1;
        else return 0;
        break;
    }
    return 0;
}