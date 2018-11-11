#include "../utils/time_profiler.hpp"
#include <cmath>
#include <omp.h>
#include <iostream>
#include <vector>
#include <stack>
#include <algorithm>

template <typename T>
class Matrix
{
public:
    Matrix(std::size_t size)
        : size_(size)
    {
        data_.resize(size_ * size_);
        std::fill(data_.begin(), data_.end(), (T)0);
    }

    Matrix(T const* data, std::size_t size)
        : size_(size)
    {
        data_.resize(size_ * size_);
        std::copy(data, data + size_ * size_, data_);
    }

    Matrix(Matrix const& other)
        : data_(other.data_)
        , size_(other.size_)
    {}

    Matrix& operator=(Matrix const& other)
    {
        data_ = other.data_;
        size_ = other.size_;
        return *this;
    }

    Matrix(Matrix && other)
        : data_(std::move(other.data_))
        , size_(std::move(other.size_))
    {}

    T& operator[](std::size_t index)
    {
        return data_[index];
    }

    T const& operator[](std::size_t index) const
    {
        return data_[index];
    }

    std::size_t GetSize() const
    {
        return size_;
    }

    T& Element(std::size_t row, std::size_t column)
    {
        return data_[row * size_ + column];
    }

    T const& Element(std::size_t row, std::size_t column) const
    {
        return data_[row * size_ + column];
    }

    Matrix<T> GetMinorSubmatrix(std::size_t row, std::size_t column) const
    {
        std::size_t submatrix_size = size_ - 1;
        Matrix<T> submatrix(submatrix_size);

        for (std::size_t i = 0; i < submatrix_size * submatrix_size; ++i)
        {
            std::size_t dst_x = i % submatrix_size;
            std::size_t dst_y = i / submatrix_size;

            std::size_t src_index = (dst_x >= column) + dst_y + (dst_y >= row) * size_ + i;

            submatrix[i] = data_[src_index];

        }

        return submatrix;
    }

    Matrix<T> GetInnerSubmatrix() const
    {
        std::size_t inner_size = size_ - 2;
        Matrix<T> inner(inner_size);

        for (std::size_t i = 0; i < inner_size * inner_size; ++i)
        {
            std::size_t dst_x = i % inner_size;
            std::size_t dst_y = i / inner_size;

            std::size_t src_index = (dst_x + 1) + (dst_y + 1) * size_;

            inner[i] = data_[src_index];

        }

        return inner;

    }

    T Determinant() const
    {
        if (size_ == 1)
        {
            return data_[0];
        }

        T result = {0};
#pragma omp parallel for reduction (+:result)
        for (int i = 0; i < GetSize(); ++i)
        {
            T value = Element(0, i) * GetMinorSubmatrix(0, i).Determinant();
            result += (i % 2) ? -value : value;
        }

        return result;
    }

private:
    std::vector<T> data_;
    std::size_t size_;

};

template <typename T>
std::ostream& operator<<(std::ostream & os, Matrix<T> const& matrix)
{
    //os << "{";
    for (std::size_t i = 0; i < matrix.GetSize(); ++i)
    {
        //os << "{";
        for (std::size_t j = 0; j < matrix.GetSize(); ++j)
        {
            os << matrix[i * matrix.GetSize() + j] << " ";//", ";
        }
        os << std::endl;
        //os << "},";
    }
   // os << "}";
    return os;
}

int main(int argc, char** argv)
{
    Matrix<float> matrix(3);
    for (int i = 0; i < matrix.GetSize(); ++i)
    {
        for (int j = 0; j < matrix.GetSize(); ++j)
        {
            matrix.Element(i, j) = (float)(std::rand() % 10);
        }
    }

    std::cout << matrix << std::endl;

    {
        PROFILE_TIME("determinant calculation");
        std::cout << matrix.Determinant() << std::endl;
    }

    return 0;
}
