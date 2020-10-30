#include <stdio.h>
#include <math.h>

using namespace std;


double eigen_fun(int k, int i, double h)
{
    return(sin(M_PI*k*i*h));
    //return(sin(M_PI*k*(i+1)*h));
}

double eigen_val(int k, double h)
{
    return( 2/(h*h) * (1 - cos(M_PI*k*h)));
}

double Mult(int k, int l, int N)
{
    double h = 2.0/(2.0*N - 1);
    double res = 0, res1 = 0, res2 = 0;
    for(int i = 1; i < N; i++)
    {
        res += eigen_fun(k, i, h) * eigen_fun(l, i, h);
        res1 += eigen_fun(k, i, h) * eigen_fun(k, i, h);
        res2 += eigen_fun(l, i, h) * eigen_fun(l, i, h);
    }
    return res/(sqrt(res1 * res2));
}


double Mult_Matrix(FILE * fout, int N)
{
    double max = 0.0;
    fprintf(fout, " Матрица скалярного произведения: \n");
    for(int i = 1; i < N; i++){
        for(int j = 1; j < N; j++){
            fprintf(fout, " %.16lf ", Mult(i, j, N));
            printf("%.16lf ", Mult(i, j, N));
            if((i != j) && fabs(Mult(i, j, N)) > fabs(max)) max = Mult(i, j, N);
        }
        fprintf(fout, "\n");
        printf("\n");
    }
    printf("\n");
    return max;
}

double Er_norm(int k, int N, double *reserv)
{
    double res = 0, res1 = 0, res2 = 0;
    double h = 2.0/(2.0*N - 1);
    for(int i = 1; i <= N - 1; ++i)
    {
        reserv[i] = fabs(((eigen_fun(k, i+1, h) - 2*eigen_fun(k, i, h) + eigen_fun(k, i-1, h))/(h*h)) + eigen_val(k, h)*eigen_fun(k, i, h));
        //printf(" %.16lf  ", reserv[i]);
        res1 += reserv[i] * reserv[i];
    }
    //printf("\n %.16lf \n",  res1);
    //printf("\n %.16lf \n", eigen_val(k, h));
    return sqrt(res1)/eigen_val(k, h);
}

int main()
{
    FILE * fout;
    int N = 10;
    //double h = 2.0/(2.0*N - 1);
    double * error;
    double * reserv;
    double maxLam = 0.0;
    double maxMatr = 0.0;

    fout = fopen("output.txt", "w");
    maxMatr = Mult_Matrix(fout, N);

    error = (double *)malloc((N-1) *sizeof(double));
    reserv = (double *)malloc((N-1) *sizeof(double));

    printf(" \n ");
    for(int t = 1; t < N; t++)
    {
        error[t] = Er_norm(t, N, reserv);
        if(fabs(error[t]) > fabs(maxLam)) maxLam = error[t];
        printf(" %.16lf ",  Er_norm(t, N, reserv));
    }
    printf(" \n ");
    printf(" \n  MaxLam = %.16lf \n",  maxLam);
    printf(" \n  MaxMatr = %.16lf \n",  maxMatr);
    fprintf(fout, " \n  MaxLam = %.16lf \n",  maxLam);
    fprintf(fout, " \n  MaxMatr = %.16lf \n",  maxMatr);

    return 0;
}


