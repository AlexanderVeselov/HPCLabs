#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <iostream>
#define N 500000000

int main(int argc, char** argv)
{
    double* b = new double[N];
    double s = 0;
#pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        b[i] = i * tan(i * 3.14 / N);
        s += b[i];
    }

    printf("result: %f (%d MB used)", s, N * sizeof(double) / (1000 * 1000));

    delete[] b;
    return 0;

}
