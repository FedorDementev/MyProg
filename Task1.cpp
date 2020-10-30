#include <stdio.h>
#include <math.h>

using namespace std;

//считаем что имеется задача y' = - A * y(x), функция y(x) = y_x;

double solution(double A, double x) {
    return (exp(-A*x));
}

void schemes(double *y, double h, double A, int N, int number) {
    y[0] = 1.0;
    //printf(" %.16lf", y[0]);
    if(number == 1) {
        for(int k = 0; k < N; k++)  {
            y[k+1] = y[k] * (1 - A * h);
            //printf(" %.16lf", y[k+1]);
        }
    }
    if(number == 2) {
        for(int k = 0; k < N; k++)  {
            y[k+1] = y[k] * (1.0/(1.0 + A*h));
             //printf(" %.16lf", y[k+1]);
        }
    }

    if(number == 3) {
        for(int k = 0; k < N; k++)  y[k+1] = y[k] * ((2.0 - A * h)/(2.0 + A * h));
    }

    y[1] = 1 - A * h;

    if(number == 4) {
        for(int k = 1; k < N; k++)  y[k+1] = y[k-1] - 2 * A * h * y[k];
    }
    if(number == 5) {
        for(int k = 1; k < N; k++)  y[k+1] = (1.0/3.0) * (4 * y[k] - (1.0 + 2 * A * h) * y[k-1]);
    }
}

double errorNorm(double *y, double *solv, int N, double A, double h)
{
    double max = 0.0;

    for(int k = 0; k < N + 1; k++) {
        solv[k] = fabs(solution(A, k * h) - y[k]);
        //printf (" \n %.16lf ", solv[k]);
        //printf (" \n lalala  %.16lf ", solution(A, k * h));
        if(fabs(solv[k]) - max > 1e-13)  max = solv[k];
    }
    return max;
}

int main() {

    //double a = 0.0;
    //double b = 10.0;
    double A, res2, res1;
    double h;
    int N, number;

    printf("Введите номер схемы number = \n");
    scanf("%d", &number);

    printf("Введите A: \n");
    scanf("%lf", &A);

    res1 = 1.0;
    res2 = 1.0;

    for(int m = 0; m < 7; m++) {
        //int m = 3;
        h = pow(10, -m);
        N = pow(10, m+1);

        double * y = new double[N+1];
        double * deviation = new double[N+1];
        schemes(y, h, A, N, number);
        //printf("\n jklytftftf \n\n\n\n");
        res2 = errorNorm(y, deviation, N, A, h);
        //printf("ghjk \n");
        //printf(" \n %.16lf \n", res2);
        printf(" \n %.16lf \n", (log(res1/res2)/log(10)));
        res1 = res2;
    }

    return 0;
}

