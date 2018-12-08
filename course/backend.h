#pragma once

#define DLLEXPORT __declspec(dllexport)

extern "C"
{
    void DLLEXPORT fft(double* arr, int size, double* out_real, double* out_imag);

}
