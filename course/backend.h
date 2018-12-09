#pragma once

#define DLLEXPORT __declspec(dllexport)

extern "C"
{
    void DLLEXPORT fft(double* arr, int size, double* out_real, double* out_imag);
    void DLLEXPORT fft_simd(double* arr, int size, double* out_real, double* out_imag);
    void DLLEXPORT fft_simd_aligned(double* arr, int size, double* out_real, double* out_imag);
    void DLLEXPORT fft_parallel_simd_aligned(double* arr, int size, double* out_real, double* out_imag);

}
