#include "backend.h"
#include <stdio.h>
#include <cmath>

void SendFunc(float* arr, int size)
{
    for (int i = 0; i < size; ++i)
    {
        printf("%f\n", arr[i]);
    }
}

void RecvFunc(float* arr)
{
    for (int i = 0; i < 100; ++i)
    {
        arr[i] = i;
    }
}
