#include <cstdlib>

// Functions for allocating/deallocating memory
double ** create_array (size_t rows, size_t cols)
{
    double ** m = new double *[rows];
    m[0] = new double[rows * cols];
    for (size_t i = 1; i != rows; ++i)
        m[i] = m[i - 1] + cols;
    return m;
}

double * create_vector (size_t a)
{
    double * v = new double[a];
    return v;
}

void free_array (double ** m)
{
    delete [] m[0];
    delete [] m;
}

void free_vector (double * v)
{
    delete [] v;
}
