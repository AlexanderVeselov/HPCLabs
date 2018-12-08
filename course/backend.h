#pragma once

#define DLLEXPORT __declspec(dllexport)

extern "C"
{
    void DLLEXPORT SendFunc(float* arr, int size);
    void DLLEXPORT RecvFunc(float* x);

}
