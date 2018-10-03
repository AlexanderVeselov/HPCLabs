#include <omp.h>
#include <iostream>

int main(int argc, char** argv)
{
#pragma omp parallel
    {
        int index = omp_get_thread_num();
        int total = omp_get_num_threads();
        printf("Hello, OpenMP! I am %d of %d\n", index, total);
    }

    return 0;
}
