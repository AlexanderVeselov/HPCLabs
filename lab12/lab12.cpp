/*

3. Задача о вычислении произведения матрицы на вектор. Дана A –
матрица размерностью m строк на n столбцов и дан вектор x размерностью n.
Исходные данные хранятся в файлах на диске, подготовленные для каждого
процессора (написать последовательную программу подготовки данных,
входным параметром которой является число процессоров). Написать про-
грамму для p процессоров (2<=p<<m,n), вычисляющую произведение Ax,
как m скалярных произведений. Распределение нагрузки на процессоры про-
вести наиболее удобным способом

*/

#include <mpi.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>

int g_rank;
int g_size;

int ThreadId()
{
    return g_rank;
}

bool MainThread()
{
    return g_rank == 0;
}

int NumThreads()
{
    return g_size;
}

void LoadMatrixFromFile(std::string const& filename, std::vector<float> & out_matrix, std::size_t & w, std::size_t & h)
{
    out_matrix.clear();
    w = 0;
    h = 0;

    std::ifstream in_file(filename);
    if (!in_file)
    {
        std::cout << "Error reading " << filename << std::endl;
        return;
    }

    while (!in_file.eof())
    {
        float value;
        std::string line;
        std::getline(in_file, line);
        std::istringstream sstr(line);

        w = 0;
        while (!sstr.eof())
        {
            sstr >> value;
            out_matrix.push_back(value);
            ++w;
        }

        ++h;

    }

}

void LoadVectorFromFile(std::string const& filename, std::vector<float> & out_vec, std::size_t & size)
{
    std::size_t w = 0;
    LoadMatrixFromFile(filename, out_vec, w, size);
    if (w != 1)
    {
        std::cout << "Error loading vector" << std::endl;
    }
}

void PrintMatrix(std::vector<float> const& matrix, std::size_t w, std::size_t h)
{
    for (std::size_t i = 0; i < h; ++i)
    {
        for (std::size_t j = 0; j < w; ++j)
        {
            std::cout << matrix[i * w + j] << " ";
        }
        std::cout << std::endl;
    }
}

void PrintVector(std::vector<float> const& vec, std::size_t size)
{
    PrintMatrix(vec, 1, size);
}

void MatVecMul(std::vector<float> const& matrix_rows, int num_rows, std::vector<float> const& vec, int vec_size, std::vector<float> & result)
{
    result.clear();
    result.resize(num_rows * vec_size);

    for (int i = 0; i < num_rows; i++)
    {
        result[i] = 0;
        for (int j = 0; j < vec_size; j++)
        {
            result[i] += matrix_rows[i * vec_size + j] * vec[j];
        }
    }
}

int GetNumRowsForProcess(int total_rows, int process_id)
{
    int rest_rows = total_rows;

    for (int i = 0; i < process_id; ++i)
    {
        rest_rows = rest_rows - rest_rows / (NumThreads() - i);
    }
    return rest_rows / (NumThreads() - process_id);

}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &g_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &g_rank);

    std::vector<float> matrix;
    std::size_t w = 0;
    std::size_t h = 0;
    std::vector<float> vec;
    std::size_t vec_size = 0;
    if (MainThread())
    {
        LoadMatrixFromFile("matrix.txt", matrix, w, h);
        //PrintMatrix(matrix, w, h);

        LoadVectorFromFile("vector.txt", vec, vec_size);
        //PrintVector(vec, vec_size);
    }

    // Broadcast width, height
    MPI_Bcast(&w, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&h, 1, MPI_INT, 0, MPI_COMM_WORLD);

    vec_size = w;

    if (!MainThread())
    {
        vec.resize(vec_size);
    }

    // Broadcast vector
    MPI_Bcast(vec.data(), h, MPI_FLOAT, 0, MPI_COMM_WORLD);

    std::vector<float> proc_rows;
    int proc_row_num = GetNumRowsForProcess(h, ThreadId());
    proc_rows.resize(w * proc_row_num);

    int* send_offset = new int[NumThreads()];
    int* send_size   = new int[NumThreads()];

    send_size[0]   = GetNumRowsForProcess(h, 0) * w;
    send_offset[0] = 0;
    for (int i = 1; i < NumThreads(); i++)
    {
        send_size[i] = GetNumRowsForProcess(h, i) * h;
        send_offset[i] = send_offset[i - 1] + send_size[i - 1];
    }

    MPI_Scatterv(matrix.data(), send_size, send_offset, MPI_FLOAT, proc_rows.data(), send_size[ThreadId()], MPI_FLOAT, 0, MPI_COMM_WORLD);

    std::vector<float> local_result;
    local_result.resize(proc_row_num);
    MatVecMul(proc_rows, proc_row_num, vec, vec_size, local_result);

    send_size[0] = GetNumRowsForProcess(h, 0);
    send_offset[0] = 0;
    for (int i = 1; i < NumThreads(); i++)
    {
        send_size[i] = GetNumRowsForProcess(h, i);
        send_offset[i] = send_offset[i - 1] + send_size[i - 1];
    }

    std::vector<float> result;
    result.resize(h);
    MPI_Gatherv(local_result.data(), send_size[ThreadId()], MPI_FLOAT, result.data(), send_size, send_offset, MPI_FLOAT, 0, MPI_COMM_WORLD);

    if (ThreadId() == 0)
    {
        PrintVector(result, h);
    }

    MPI_Finalize();

    return 0;
}
