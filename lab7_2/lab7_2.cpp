#include <cmath>
#include <omp.h>
#include <iostream>
#include <vector>
#include <array>

template <typename T, std::size_t size>
class Matrix
{
public:
    Matrix()
    {
        std::fill(data_.begin(), data_.end(), 0);
    }

    Matrix(T const* data)
    {
        std::copy(data, data + size * size, data_);
    }

    Matrix(Matrix const& other)
        : data_(other.data_)
    {}

    Matrix(Matrix && other)
        : data_(std::move(other.data_))
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
        return size;
    }

    T& Get(std::size_t row, std::size_t column)
    {
        return data_[row * size + column];
    }

    T const& Get(std::size_t row, std::size_t column) const
    {
        return data_[row * size + column];
    }

    Matrix<T, size - 1> GetCofactor(std::size_t row, std::size_t column) const
    {
        Matrix<T, size - 1> cofactor;
        std::size_t dst_row = 0;
        for (std::size_t src_row = 0; src_row < size; ++src_row)
        {
            if (src_row == row)
            {
                continue;
            }

            std::size_t dst_column = 0;
            for (std::size_t src_column = 0; src_column < size; ++src_column)
            {
                if (src_column == column)
                {
                    continue;
                }

                cofactor.Get(dst_row, dst_column) = Get(src_row, src_column);

                ++dst_column;
            }

            ++dst_row;
        }

        return cofactor;
    }

private:
    std::array<T, size * size> data_;

};

template <typename T, std::size_t size>
std::ostream& operator<<(std::ostream & os, Matrix<T, size> const& matrix)
{
    for (std::size_t i = 0; i < size; ++i)
    {
        for (std::size_t j = 0; j < size; ++j)
        {
            os << matrix[i * size + j] << " ";
        }
        os << std::endl;
    }
    return os;
}

int main(int argc, char** argv)
{
    Matrix<std::uint32_t, 3> matrix;
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            matrix.Get(i, j) = rand();
        }
    }

    std::cout << matrix << std::endl;
    std::cout << matrix.GetCofactor(1, 0) << std::endl;

    return 0;
}
