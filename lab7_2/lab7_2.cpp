#include <cmath>
#include <omp.h>
#include <iostream>
#include <vector>

template <typename T>
class Matrix
{
public:
    Matrix(std::size_t size)
        : size_(size)
    {
        data_.resize(size_ * size_);
        std::fill(data_.begin(), data_.end(), 0);
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

    T Determinant() const
    {
        if (size_ == 1)
        {
            return data_[0];
        }

        T result = {0};

        for (std::size_t i = 0; i < GetSize(); ++i)
        {
            T element = Element(0, i);
            T det = GetCofactor(0, i).Determinant();
            T value = element * det;
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
    os << "{";
    for (std::size_t i = 0; i < matrix.GetSize(); ++i)
    {
        os << "{";
        for (std::size_t j = 0; j < matrix.GetSize(); ++j)
        {
            os << matrix[i * matrix.GetSize() + j] << ", ";
        }
        os << "},";
    }
    os << "}";
    return os;
}

int main(int argc, char** argv)
{
    Matrix<std::int64_t> matrix(4);
    for (int i = 0; i < 4; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            matrix.Element(i, j) = rand();
        }
    }

    std::cout << matrix << std::endl;
    std::cout << matrix.Determinant() << std::endl;

    return 0;
}
