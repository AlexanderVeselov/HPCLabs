#include "backend.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <stdio.h>
#include <complex>
#include <vector>
#include <algorithm>
#include <omp.h>

std::uint32_t reverseBits(std::uint32_t i) {
    std::uint32_t mask = 0x55555555; // 0101...
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

int lg(std::uint32_t i) {
    int count = -1;
    while (i) {
        i = i >> 1;
        count++;
    }
    return count;
}

// https://equilibriumofnothing.wordpress.com/2013/10/14/algorithm-iterative-fft/
// Assume that arrays sizes are even power of two
void fft_impl(std::vector<std::complex<double>> const& input, std::vector<std::complex<double>> & output, const bool inverse)
{
    const int N = input.size();
    const int absP = lg(N);

    // bottom level of iteration tree
#pragma loop(no_vector)
    for (int i = 0; i < N; i++)
    {
        output[i] = input[reverseBits(i) >> (32 - absP)];
    }

    // there are absP levels above the bottom
#pragma loop(no_vector)
    for (int p = 1; p <= absP; p++)
    {
        // complex root of unity
        const int unityStep = 0x1 << p;
        const double theta = (inverse ? -1 : 1) * 2 * M_PI / unityStep;
        const std::complex<double> unityRoot(cos(theta), sin(theta));

        // each higher level doubles the step size
#pragma loop(no_vector)
        for (int offset = 0; offset < N; offset += unityStep)
        {
            std::complex<double> omega = 1;

            // combine within a step segment (note only iterate over half step)
#pragma loop(no_vector)
            for (int k = 0; k < unityStep / 2; k++)
            {
                const std::complex<double> u = output[offset + k];

                const std::complex<double> t = omega * output[offset + k + unityStep / 2];
                omega *= unityRoot;

                output[offset + k] = u + t;
                output[offset + k + unityStep / 2] = u - t;
            }
        }
    }

    if (inverse)
    {
#pragma loop(no_vector)
        for (int j = 0; j < N; j++)
        {
            output[j] /= N;
        }
    }

}

void fft(double* arr, int size, double* out_real, double* out_imag)
{
    printf("call " __FUNCTION__ "\n");

    std::vector<std::complex<double> > array(size);
    for (int i = 0; i < size; ++i)
    {
        array[i] = std::complex<double>(arr[i], 0.0);
    }

    std::vector<std::complex<double> > result(size);
    fft_impl(array, result, false);

    for (int i = 0; i < size; ++i)
    {
        out_real[i] = result[i].real();
        out_imag[i] = result[i].imag();
    }

}

void fft_simd_impl(std::vector<std::complex<double>> const& input, std::vector<std::complex<double>> & output, const bool inverse)
{
    const int N = input.size();
    const int absP = lg(N);

    for (int i = 0; i < N; i++)
    {
        output[i] = input[reverseBits(i) >> (32 - absP)];
    }

    for (int p = 1; p <= absP; p++)
    {
        const int unityStep = 0x1 << p;
        const double theta = (inverse ? -1 : 1) * 2 * M_PI / unityStep;
        const std::complex<double> unityRoot(cos(theta), sin(theta));

        for (int offset = 0; offset < N; offset += unityStep)
        {
            std::complex<double> omega = 1;

            for (int k = 0; k < unityStep / 2; k++)
            {
                const std::complex<double> u = output[offset + k];

                const std::complex<double> t = omega * output[offset + k + unityStep / 2];
                omega *= unityRoot;

                output[offset + k] = u + t;
                output[offset + k + unityStep / 2] = u - t;
            }
        }
    }

    if (inverse)
    {
        for (int j = 0; j < N; j++)
        {
            output[j] /= N;
        }
    }

}

void fft_simd(double* arr, int size, double* out_real, double* out_imag)
{
    printf("call " __FUNCTION__ "\n");

    std::vector<std::complex<double> > array(size);
    for (int i = 0; i < size; ++i)
    {
        array[i] = std::complex<double>(arr[i], 0.0);
    }

    std::vector<std::complex<double> > result(size);
    fft_simd_impl(array, result, false);

    for (int i = 0; i < size; ++i)
    {
        out_real[i] = result[i].real();
        out_imag[i] = result[i].imag();
    }

}

template <typename T, std::size_t N = 16>
class AlignmentAllocator {
public:
    typedef T value_type;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;

    typedef T * pointer;
    typedef const T * const_pointer;

    typedef T & reference;
    typedef const T & const_reference;

public:
    inline AlignmentAllocator() throw () { }

    template <typename T2>
    inline AlignmentAllocator(const AlignmentAllocator<T2, N> &) throw () { }

    inline ~AlignmentAllocator() throw () { }

    inline pointer adress(reference r) {
        return &r;
    }

    inline const_pointer adress(const_reference r) const {
        return &r;
    }

    inline pointer allocate(size_type n) {
        return (pointer)_aligned_malloc(n * sizeof(value_type), N);
    }

    inline void deallocate(pointer p, size_type) {
        _aligned_free(p);
    }

    inline void construct(pointer p, const value_type & wert) {
        new (p) value_type(wert);
    }

    inline void destroy(pointer p) {
        p->~value_type();
    }

    inline size_type max_size() const throw () {
        return size_type(-1) / sizeof(value_type);
    }

    template <typename T2>
    struct rebind {
        typedef AlignmentAllocator<T2, N> other;
    };

    bool operator!=(const AlignmentAllocator<T, N>& other) const {
        return !(*this == other);
    }

    // Returns true if and only if storage allocated from *this
    // can be deallocated from other, and vice versa.
    // Always returns true for stateless allocators.
    bool operator==(const AlignmentAllocator<T, N>& other) const {
        return true;
    }
};

template <class T>
using aligned_vector = std::vector<T, AlignmentAllocator<T, 32>>;

void fft_simd_aligned_impl(aligned_vector<std::complex<double>> const& input, aligned_vector<std::complex<double>> & output, const bool inverse)
{
    const int N = input.size();
    const int absP = lg(N);

    for (int i = 0; i < N; i++)
    {
        output[i] = input[reverseBits(i) >> (32 - absP)];
    }

    for (int p = 1; p <= absP; p++)
    {
        const int unityStep = 0x1 << p;
        const double theta = (inverse ? -1 : 1) * 2 * M_PI / unityStep;
        const std::complex<double> unityRoot(cos(theta), sin(theta));

        for (int offset = 0; offset < N; offset += unityStep)
        {
            std::complex<double> omega = 1;

            for (int k = 0; k < unityStep / 2; k++)
            {
                const std::complex<double> u = output[offset + k];

                const std::complex<double> t = omega * output[offset + k + unityStep / 2];
                omega *= unityRoot;

                output[offset + k] = u + t;
                output[offset + k + unityStep / 2] = u - t;
            }
        }
    }

    if (inverse)
    {
        for (int j = 0; j < N; j++)
        {
            output[j] /= N;
        }
    }

}

void fft_simd_aligned(double* arr, int size, double* out_real, double* out_imag)
{
    printf("call " __FUNCTION__ "\n");


    aligned_vector<std::complex<double> > array(size);
    for (int i = 0; i < size; ++i)
    {
        array[i] = std::complex<double>(arr[i], 0.0);
    }

    aligned_vector<std::complex<double> > result(size);
    fft_simd_aligned_impl(array, result, false);

    for (int i = 0; i < size; ++i)
    {
        out_real[i] = result[i].real();
        out_imag[i] = result[i].imag();
    }

}

void fft_parallel_simd_aligned_impl(aligned_vector<std::complex<double>> const& input, aligned_vector<std::complex<double>> & output, const bool inverse)
{
    const int N = input.size();
    const int absP = lg(N);

#pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        output[i] = input[reverseBits(i) >> (32 - absP)];
    }

//#pragma omp parallel for
    for (int p = 1; p <= absP; p++)
    {
        const int unityStep = 0x1 << p;
        printf("unityStep: %d\n", unityStep);

        const double theta = (inverse ? -1 : 1) * 2 * M_PI / unityStep;
        printf("theta: %f\n", theta);

        const std::complex<double> unityRoot(cos(theta), sin(theta));

        for (int offset = 0; offset < N; offset += unityStep)
        {
            std::complex<double> omega = 1;

            for (int k = 0; k < unityStep / 2; k++)
            {
                const std::complex<double> u = output[offset + k];

                const std::complex<double> t = omega * output[offset + k + unityStep / 2];
                omega *= unityRoot;
                output[offset + k] = u + t;
                output[offset + k + unityStep / 2] = u - t;
            }
        }
    }

    if (inverse)
    {
        for (int j = 0; j < N; j++)
        {
            output[j] /= N;
        }
    }

}

void fft_parallel_simd_aligned(double* arr, int size, double* out_real, double* out_imag)
{
    printf("call " __FUNCTION__ "\n");

    aligned_vector<std::complex<double> > array(size);
    for (int i = 0; i < size; ++i)
    {
        array[i] = std::complex<double>(arr[i], 0.0);
    }

    aligned_vector<std::complex<double> > result(size);
    fft_parallel_simd_aligned_impl(array, result, false);

    for (int i = 0; i < size; ++i)
    {
        out_real[i] = result[i].real();
        out_imag[i] = result[i].imag();
    }

}

