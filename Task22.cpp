#include <stdio.h>
#include <math.h>

using namespace std;

double eigen_fun(int k, int i, double h)
{
    return(sin(M_PI*k*i*h));
    //return(sin(M_PI*k*(i+1)*h));
}

double eigen_val(int k, double h, double b)
{
    return( 2/(h*h) * (1 - cos(M_PI*k*h)) + b);
}

void f2C(double *C, double *f, double b, int n)
{
    double sum = 0;
    double sum1 = 0;
    double h = 2.0/(2*n - 1);

    for (int k = 1; k <= n-1; ++k)
    {
        for (int i = 1; i <= n-1; ++i)
            sum += f[i] * eigen_fun(k, i, h) * h;

        for (int i = 1; i <= n-1; ++i)
            sum1 += eigen_fun(k, i, h) * eigen_fun(k, i, h) * h;
        C[k] = sum/(sum1 * eigen_val(k, h, b));

        sum = 0;
        sum1 = 0;
    }
}

void C2f(double *C, double *newf, int n)
{
    double sum = 0;
    double h = 2.0/(2.0*n - 1);
    for (int i = 1; i <= n; ++i)
    {

        for (int m = 1; m <= n-1; ++m)
        {
            sum += C[m] * eigen_fun(m, i, h);
        }

        newf[i] = sum;
        sum = 0;
    }
}


void run(double * a, double *b, double *c, double *alfa, double *betta, double *newY, double *f, int N) {
    //double h = 2.0/(2.0*N - 1);
    alfa = (double *)malloc((N+1) *sizeof(double));
    betta = (double *)malloc((N+1) *sizeof(double));


    alfa[0] = 0;
    betta[0] = 0;
    alfa[1] = b[0]/c[0];
    betta[1] = f[0]/c[0];

    for(int i = 1; i < N; i++){
        alfa[i+1] = b[i] / (c[i] - a[i] * alfa[i]);
        betta[i+1] = (f[i] + a[i] * betta[i]) / (c[i] - a[i] * alfa[i]);
    }
    newY[N] = (f[N] + a[N] * betta[N]) / (c[N] - a[N] * alfa[N]);
    for(int i = (N-1); i >= 0; i--)
        newY[i] = alfa[i+1] * newY[i+1] + betta[i+1];

    free(alfa);
    free(betta);
}


double solution(double x) {     //точное решение уравнения.
    //return 3;
    return (x * x -  x);
    //return sin(3*M_PI*x);
    //return(x*x - 3 * x);
}

double funB(double x) {         //функция b(x) для прогонки. в методе Фурье мы положим её константой.
    //return (M_PI * M_PI);
    return (1);
    //return(x*x);
}

double funF(double x) {         //правая часть уравнения в задаче Коши.
    return x*x - x - 2;
    //return (x*x*x - 2*x*x - 2);
    //return (10*M_PI*M_PI*sin(3*M_PI*x));
    //return (x*x*x*x - 3*x*x*x - 2);
}

void massiv(double *a, double *b, double *c, double *B_x, double *Y, double *f, int N) {
    double h = 2.0 / (2.0 * N - 1 );
    a[0] = 0;
    for(int i = 1; i < N; i++) {
        a[i] = 1 /( h*h );
        b[i] = 1 / (h*h);
        c[i] = 2 / (h*h) + funB(i*h);
    }
    a[N] = -1;
    c[0] = 1;
    c[N] = 1;
    b[0] = 0;
    for(int k = 0; k <= N; k++) {
        B_x[k] = funB(k*h);
        Y[k] = solution(k*h);
        f[k] = funF(k*h);
    }
    f[N] = 0;
    f[0] = 0;
}

double errorNorm(double *x, double *y, double *error, int N) {   //считает вектор разности векторов x и y , выводит норму вектора разности.
    double h = 2.0 / (2.0*N - 1);
    double res = 0.0;
    for(int i = 0; i <= N; i++){
        error[i] = fabs(x[i] - y[i]);
        //printf(" %.10lf  ", error[i]);
        res += h * (x[i] - y[i]) * (x[i] - y[i]);
    }
    return sqrt(res);
}

int main() {

    int N;
    printf("введете начальное значение N:  \n");
    scanf("%d", &N);
    double h = 2.0 / (2.0 * N - 1 );
    printf(" N = %d, h =  %.16lf \n\n", N, h);

    double res1, res2, res3 = 0;
    double *a, *c, *b, *Y, *C, *f, *B_x, *newY1, *newY2, *error, *alfa, *betta;

    //double bbb = M_PI * M_PI ;   //константа для метода фурье
    double bbb = 1;

    a = (double *)malloc((N+1) *sizeof(double));
    c = (double *)malloc((N+1) *sizeof(double));
    b = (double *)malloc((N) *sizeof(double));

    C = (double *)malloc((N+1) *sizeof(double));
    Y = (double *)malloc((N+1) *sizeof(double));
    newY1 = (double *)malloc((N+1) *sizeof(double));
    newY2 = (double *)malloc((N+1) *sizeof(double));
    B_x = (double *)malloc((N+1) *sizeof(double));
    f = (double *)malloc((N+1) *sizeof(double));
    error = (double *)malloc((N+1) *sizeof(double));

    massiv(a, b, c, B_x, Y, f, N);
    run(a, b, c, alfa, betta, newY2, f, N);
    res2 = errorNorm(newY2, Y, error, N);
/*
    printf("\n  c:  \n");
    for(int j = 0; j<=N; j++) printf(" %.10lf  ", c[j]);
    printf("\n");
    printf("\n  a:  \n");
    for(int j = 0; j<=N; j++) printf(" %.10lf  ", a[j]);
    printf("\n");
    printf("\n  b:  \n");
    for(int j = 0; j<N; j++) printf(" %.10lf  ", b[j]);
*/
    printf(" для прогонки res2 =  %.16lf \n", res2);

    f2C(C, f, bbb, N);

    newY1[0] = 0;
    C2f(C, newY1, N);

    //for(int i = 0; i < N; i++) printf(" %.16lf  ", newY1[i]);
    //for(int i = 0; i <= N; i++) printf(" __%.16lf  %.16lf__ ", Y[i], newY1[i]);
    //printf("\n\n");
    res1 = errorNorm(newY1, Y, error, N);

    //for(int i = 0; i < N; i++) printf(" %.16lf  ", error[i]);
    printf("\n для Фурье res1 =  %.16lf \n", res1);


    //printf("\n  вектор разницы двух решений: \n ");
    for(int i = 0; i <= N; i++) {
        //printf(" %.16lf  \n ", fabs(newY1[i] - newY2[i]));
        if(res3 <=  fabs(newY1[i] - newY2[i])) res3 = fabs(newY1[i] - newY2[i]);
    }
    printf("\n максимальное отклонение между решениями =  %.16lf  \n \n", res3);

    free(a);
    free(b);
    free(c);
    free(Y);
    free(f);
    free(C);
    free(B_x);
    free(newY1);
    free(newY2);
    free(error);

    return 0;
}
