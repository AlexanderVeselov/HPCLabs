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

    Matrix<T> GetCofactor(std::size_t row, std::size_t column) const
    {
        Matrix<T> cofactor(size_ - 1);
        std::size_t dst_row = 0;
        for (std::size_t src_row = 0; src_row < size_; ++src_row)
        {
            if (src_row == row)
            {
                continue;
            }

            std::size_t dst_column = 0;
            for (std::size_t src_column = 0; src_column < size_; ++src_column)
            {
                if (src_column == column)
                {
                    continue;
                }

                cofactor.Element(dst_row, dst_column) = Element(src_row, src_column);

                ++dst_column;
            }

            ++dst_row;
        }

        return cofactor;
    }

    Matrix<T> GetInner() const
    {
        Matrix<T> inner(size_ - 2);
        for (std::size_t i = 1; i < size_ - 1; ++i)
        {
            for (std::size_t j = 1; j < size_ - 1; ++j)
            {
                inner.Element(i - 1, j - 1) = Element(i, j);
            }
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
            T value = Element(0, i) * GetCofactor(0, i).Determinant();
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
    Matrix<float> matrix(5);
    for (int i = 0; i < matrix.GetSize(); ++i)
    {
        for (int j = 0; j < matrix.GetSize(); ++j)
        {
            matrix.Element(i, j) = (float)(std::rand() % 30000) + 1;
        }
    }

    std::cout << matrix << std::endl;
    {
        PROFILE_TIME("determinant calculation");
        std::cout << matrix.Determinant() << std::endl;
    }

    {
        PROFILE_TIME("condensation method");
        Matrix<float> prevprev = matrix;
        Matrix<float> prev = matrix;
        for (int iter = 0; iter < matrix.GetSize() - 1; ++iter)
        {
            Matrix<float> current(prev.GetSize() - 1);
            for (int i = 0; i < prev.GetSize() - 1; ++i)
            {
                for (int j = 0; j < prev.GetSize() - 1; ++j)
                {
                    current.Element(i, j) = prev.Element(i, j) * prev.Element(i + 1, j + 1)
                        - prev.Element(i, j + 1) * prev.Element(i + 1, j);
                }
            }
            if (iter > 0)
            {
                Matrix<float> inner = prevprev.GetInner();
                //std::cout << current << std::endl;
                //std::cout << inner << std::endl;
                //std::cout << prev << std::endl;
                //std::cout << prevprev << std::endl;
                for (int i = 0; i < inner.GetSize(); ++i)
                {
                    for (int j = 0; j < inner.GetSize(); ++j)
                    {
                        //std::cout << current.Element(i, j) << " by " << inner.Element(i, j) << std::endl;
                        current.Element(i, j) /= inner.Element(i, j);

                    }
                }
                prevprev = prev;
            }
            prev = current;
        }
        std::cout << prev << std::endl;
    }
    return 0;
}
