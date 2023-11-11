#include <cstdlib>

// Functions for allocating/deallocating memory
double ** create_array (size_t a, size_t b)
{
    double ** m = new double *[a];
    m[0] = new double[a * b];
    for (size_t i = 1; i != a; ++i)
        m[i] = m[i - 1] + b;
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
