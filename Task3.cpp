#include <stdio.h>
#include <math.h>

void run(double * a, double *b, double *c, double *newY, double *f, int N) {

    double *alfa, *betta;

    alfa = new double[N+1];
    betta = new double[N+1];

    alfa[0] = 0;
    betta[0] = 0;
    alfa[1] = b[0]/c[0];
    betta[1] = f[0]/c[0];

    for(int i = 1; i < N; i++)
        alfa[i+1] = b[i] / (c[i] - a[i] * alfa[i]);
    for(int i = 1; i < N; i++)
        betta[i+1] = (f[i] + a[i] * betta[i]) / (c[i] - a[i] * alfa[i]);

    newY[N] = (f[N] + a[N] * betta[N]) / (c[N] - a[N] * alfa[N]);

    for(int i = (N-1); i >= 0; i--)
        newY[i] = alfa[i+1] * newY[i+1] + betta[i+1];

    delete []alfa;
    delete []betta;
}

double b_m(double x) {
    //return x;
    return exp(-x);
}

double exact(double t, double x) {
    //return exp(-t) * sin(3 * M_PI * x);
    return exp(t) * x * (1 - x*x);
}

double right_part(double t, double x) {
    //return exp(-t) * sin(3 * M_PI * x) * (x - 1 - 9 * M_PI * M_PI);
    return exp(t) * x * ((1-x*x) + 6 + exp(-x) * (1 - x*x));
}

double MaxFabs(FILE * inp, int N , int M) {

    double tau = 1.0 / N;
    double h = 2.0 / (2.0 * M - 1);

    double max = 0.0;
    double res;
    double error;
    for (int n = 0; n < N + 1; n++)
        for (int m = 0; m < M + 1; m++)
        {
            fscanf(inp, "%lf", &res);
            error = fabs(res - exact(n * tau, m * h));
            if (error > max)
                max = error;
        }
    return max;
}

void explicit_fun(double* prev, double* next,  int M, int N, int n) {

    double tau = 1.0 / N;
    double h = 2.0 / (2.0 * M - 1);

    for(int m = 1; m < M; m++)
        next[m] = (tau/(h*h)) * (prev[m+1] - 2*prev[m] + prev[m-1]) - tau * prev[m] * b_m(h * m) + tau * right_part(tau*n, h*m) + prev[m] ;

    next[0] = 0;
    next[M] = -next[M-1];
}

void implicit_fun(double* prev, double* next,  int M, int N, int n, double* a, double* b, double* c, double* f) {

    double tau = 1.0 / N;
    double h = 2.0 / (2.0 * M - 1);

    a[0] = 0;
    for(int i = 1; i < M; i++) {
        a[i] = 1 / (h * h);
        b[i] = 1 / (h * h);
        c[i] = 1 / tau + 2 / (h * h) + b_m(h * i);
        f[i] = right_part(tau * (n + 1), h * i) + 1/tau * prev[i];
    }
    a[M] = -1;
    c[0] = 1;
    c[M] = 1;
    b[0] = 0;
    f[0] = 0;
    f[M] = 0;

    run(a, b, c, next, f, M);
}

void krank_nikolson(double* prev, double* next,  int M, int N, int n, double* a, double* b, double* c, double* f) {

    double tau = 1.0 / N;
    double h = 2.0 / (2.0 * M - 1);

    a[0] = 0;
    for(int i = 1; i < M; i++) {
        a[i] = 1 / (2 * h * h);
        b[i] = 1 / (2 * h * h);
        c[i] = 1 / tau + 1 / (h * h) + b_m(h * i) / 2;
        f[i] = (1 / (2 * h * h)) * (prev[i+1] - 2 * prev[i] + prev[i-1]) + prev[i] * (1 / tau  - b_m(h * i) / 2) + (right_part(tau * (n + 1), h * i) +right_part(tau * n, h * i)) / 2;
    }
    a[M] = -1;
    c[0] = 1;
    c[M] = 1;
    b[0] = 0;
    f[0] = 0;
    f[M] = 0;

    run(a, b, c, next, f, M);
}

int main() {

    FILE* out;
    FILE* inp;
    FILE* result;

    int N, M, k, num;
    double max, tau, h;

    double* prev;
    double* next;
    double* step; //для перемены ролей prev и next
    double* a;
    double* b;
    double* c;
    double* f;

    printf("Явная схема - 1;\nНеявная - 2;\nСхема Кранка-Николсона - 3  \n");
    scanf("%d", &num);
    result = fopen("result.txt","w");

    for (k = 20; k < 301; k += 25) {   //будем считать сразу для разноко количества узлов по t и x

        out = fopen("output.txt", "w");

        if (num == 1)
        {
            N = (3*k*k);
            M = k;
        }
        else {
            M = k;
            N = k;
        }

        tau = 1.0 / N;
        h = 2.0 / (2.0 * M - 1);

        prev = new double[M+1];
        next = new double[M+1];
        a = new double[M+1];
        b = new double[M];
        c = new double[M+1];
        f = new double[M+1];
        //память выделил, теперь занесём начальные данные
        for(int m = 0; m <= M; m++)
        {
            prev[m] = exact(0, m * h);
            fprintf(out, " %.16lf ", prev[m]);  //создали массив сеточных значений в момент времени t = 0;
        }
        fprintf(out, "\n");

        for (int n = 0; n < N; n++)  //цикл по слоям времени
        {
            if (num == 1)
                explicit_fun(prev, next, M, N, n);
            if (num == 2)
                implicit_fun(prev, next, M, N, n, a, b, c, f);
            if (num == 3)
                krank_nikolson(prev, next, M, N, n, a, b, c, f);

            for (int m = 0; m <= M; m++)
            fprintf(out, " %.16lf ", next[m]);
            fprintf(out, "\n");
            //теперь меняю указатель чтобы считать следующий слой;
            step = next;
            next = prev;
            prev = step;
        }
        fclose(out);
        //файл output.txt готов и теперь пора считать норму;
        inp = fopen("output.txt", "r"); //открыли на чтение;

        max = MaxFabs(inp, N, M);

        if (num == 1)
            fprintf(result, "%d %.16lf\n", M, max /(tau + h * h)); // явная схема

        if (num == 2)
            fprintf(result, "%d %.16lf\n", M, max /(tau + h * h)); // неявная схема

        if (num == 3)
            fprintf(result, "%d %.16lf\n", M, max /(tau * tau + h * h)); // Кранка-Николсона


        fclose(inp);

        delete []prev;
        delete []next;
        delete []a;
        delete []b;
        delete []c;
        delete []f;
    }
    return 0;
}
