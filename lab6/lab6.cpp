#include "../utils/time_profiler.hpp"

#include <omp.h>
#include <iostream>

// ВАРИАНТ 3. Дана последовательность натуральных чисел { a0, …, an–1 }. Создать
// OpenMP - приложение для вычисления общей суммы и всех промежуточных
// сумм простых чисел последовательности.

#include <vector>
#include <algorithm>

///
/// Sequential Algorithms
///
template <typename T>
T Reduce(std::vector<T> const& src)
{
    PROFILE_TIME("sequential reduce");

    T sum = 0;

    for (int i = 0; i < src.size(); ++i)
    {
        sum += src[i];
    }

    return sum;
}

template <typename T, class Pr>
std::vector<T> Filter(std::vector<T> const& src, Pr predicate)
{
    PROFILE_TIME("sequential filter");

    std::vector<T> result;

    for (int i = 0; i < src.size(); ++i)
    {
        if (predicate(src[i]))
        {
            result.push_back(src[i]);
        }
    }

    return result;
}

template <typename T>
std::vector<T> Scan(std::vector<T> const& src)
{
    PROFILE_TIME("sequential scan");

    std::vector<T> result(src.size());
    result[0] = src[0];

    for (std::size_t i = 1; i < src.size(); ++i)
    {
        result[i] += result[i - 1] + src[i];
    }

    return result;
}


///
/// Parallel algorithms
///

std::uint32_t ReduceOmp(std::vector<std::uint32_t> const& src)
{
    PROFILE_TIME("parallel reduce");

    std::uint32_t sum = 0;

#pragma omp parallel for reduction (+:sum)
    for (int i = 0; i < src.size(); ++i)
    {
        sum += src[i];
    }

    return sum;
}

template <typename T, class Pr>
std::vector<T> FilterOmpNonOrdered(std::vector<T> const& src, Pr predicate)
{
    PROFILE_TIME("parallel filter non odered");

    std::vector<T> result;

#pragma omp parallel for
    for (int i = 0; i < src.size(); ++i)
    {
        if (predicate(src[i]))
        {
#pragma omp critical
            result.push_back(src[i]);
        }
    }

    return result;
}

template <typename T, class Pr>
std::vector<T> FilterOmp(std::vector<T> const& src, Pr predicate)
{
    PROFILE_TIME("parallel filter");

    std::vector<T> result;
    std::vector<std::size_t> thread_counts;

#pragma omp parallel
    {
        std::vector<T> thread_result;
        const int curr_thread = omp_get_thread_num();
        const int num_threads = omp_get_num_threads();
#pragma omp single
        {
            thread_counts.resize(num_threads);
        }

#pragma omp for
        for (int i = 0; i < src.size(); ++i)
        {
            if (predicate(src[i]))
            {
                thread_result.push_back(src[i]);
                ++thread_counts[curr_thread];
            }
        }
#pragma omp barrier

#pragma omp single
        {
            thread_counts = Scan(thread_counts);
            result.resize(thread_counts.back());
        }

        std::size_t offset = 0;
        if (curr_thread > 0)
        {
            offset = thread_counts[curr_thread - 1];
        }
        std::copy(thread_result.begin(), thread_result.end(), result.begin() + offset);

    }

    return result;
}

template <typename T>
std::vector<T> ScanOmp(std::vector<T> const& src)
{
    PROFILE_TIME("parallel scan");

    std::vector<T> result(src.size());
    std::vector<T> thread_sums;

#pragma omp parallel
    {
        const int curr_thread = omp_get_thread_num();
        const int num_threads = omp_get_num_threads();
#pragma omp single
        {
            thread_sums.resize(num_threads);
        }

        T local_sum = 0;
#pragma omp for
        for (int i = 0; i < src.size(); i++)
        {
            local_sum += src[i];
            result[i] = local_sum;
        }

        thread_sums[curr_thread] = local_sum;
#pragma omp barrier
        T local_sum_offset = 0;
        for (int i = 0; i < curr_thread; i++)
        {
            local_sum_offset += thread_sums[i];
        }

#pragma omp for
        for (int i = 0; i < src.size(); i++)
        {
            result[i] += local_sum_offset;
        }
    }

    return result;
}

template <typename T>
void PrintVector(std::vector<T> const& data)
{
    return;
    for (auto i : data)
    {
        std::cout << i << " ";
    }
    std::cout << std::endl;
}

int main(int argc, char** argv)
{
    std::srand(0);
    constexpr std::size_t input_data_length = 20000000;
    std::vector<std::uint32_t> input_data(input_data_length);

    std::uint32_t value = 0;
    std::generate(input_data.begin(), input_data.end(), [&value]() {
        //return std::rand() + 1;
        //std::cout << value + 1 << std::endl;
        return ++value;
    });

    auto is_prime = [](std::uint32_t x)
    {
        for (std::uint32_t q = 2; q*q <= x; ++q)
        {
            if (!(x % q))
            {
                return false;
            }
        }

        return true;
    };

    // Sequential algorithm test
    std::uint32_t seq_reduced = Reduce(input_data);
    std::cout << "Sequential reduction: sum = " << seq_reduced << std::endl;

    auto seq_filtered = Filter(input_data, is_prime);
    PrintVector(seq_filtered);

    auto seq_scanned = Scan(seq_filtered);
    PrintVector(seq_scanned);

    std::cout << std::endl;

    // Parallel algorithm test
    std::uint32_t parallel_sum = ReduceOmp(input_data);
    std::cout << "Parallel reduction: sum = " << parallel_sum << std::endl;

    auto parallel_filtered = FilterOmp(input_data, is_prime);
    PrintVector(parallel_filtered);

    auto parallel_filtered_non_ordered = FilterOmpNonOrdered(input_data, is_prime);
    PrintVector(parallel_filtered_non_ordered);

    auto parallel_scanned = ScanOmp(parallel_filtered);
    PrintVector(parallel_scanned);

    return 0;
}
