#include "backend.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <stdio.h>
#include <complex>
#include <vector>
#include <algorithm>

uint32_t reverseBits(uint32_t i) {
    register uint32_t mask = 0x55555555; // 0101...
    i = ((i & mask) << 1) | ((i >> 1) & mask);
    mask = 0x33333333; // 0011...
    i = ((i & mask) << 2) | ((i >> 2) & mask);
    mask = 0x0f0f0f0f; // 00001111...
    i = ((i & mask) << 4) | ((i >> 4) & mask);
    mask = 0x00ff00ff; // 0000000011111111...
    i = ((i & mask) << 8) | ((i >> 8) & mask);
    // 00000000000000001111111111111111 no need for mask
    i = (i << 16) | (i >> 16);
    return i;
}

int lg(uint32_t i) {
    int count = -1;
    while (i) {
        i = i >> 1;
        count++;
    }
    return count;
}

// https://equilibriumofnothing.wordpress.com/2013/10/14/algorithm-iterative-fft/
// Assume that arrays sizes are even power of two
void iterativeFFT(std::vector<std::complex<double>> const& primal,
    std::vector<std::complex<double>> & dual, const int P) {
    const int N = primal.size();
    const bool inverse = P < 0;
    const int absP = inverse ? -P : P;

    // bottom level of iteration tree
    for (int i = 0; i < N; i++)
    {
        dual[i] = primal[reverseBits(i) >> (32 - absP)];
    }

    // there are absP levels above the bottom
    for (int p = 1; p <= absP; p++) {
        // complex root of unity
        const int unityStep = 0x1 << p;
        const double theta = (inverse ? -1 : 1) * 2 * M_PI / unityStep; // INVERSE
        const std::complex<double> unityRoot(cos(theta), sin(theta));

        // each higher level doubles the step size
        for (int offset = 0; offset < N; offset += unityStep)
        {
            std::complex<double> omega = 1;

            // combine within a step segment (note only iterate over half step)
            for (int k = 0; k < unityStep / 2; k++)
            {
                const std::complex<double> u = dual[offset + k];

                const std::complex<double> t = omega * dual[offset + k + unityStep / 2];
                omega *= unityRoot;

                dual[offset + k] = u + t;
                dual[offset + k + unityStep / 2] = u - t;
            }
        }
    }

    if (inverse)
    {
        for (int j = 0; j < N; j++)
        {
            dual[j] /= N;
        }
    }

}

void fft(double* arr, int size, double* out_real, double* out_imag)
{
    std::vector<std::complex<double> > array(size);
    for (int i = 0; i < size; ++i)
    {
        array[i] = std::complex<double>(arr[i], 0.0);
    }

    std::vector<std::complex<double> > result(size);

    try
    {
        iterativeFFT(array, result, lg(size));
    }
    catch (std::exception & ex)
    {
        printf("%s\n", ex.what());
    }

    for (int i = 0; i < size; ++i)
    {
        out_real[i] = result[i].real();
        out_imag[i] = result[i].imag();
    }

}
