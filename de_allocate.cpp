#include <cstdlib>
#include "headers/reiman_functions.hpp"
#include "headers/de_allocate.hpp"

// Functions for allocating/deallocating memory
double ** create_array (size_t a, size_t b)
{
    double ** m = new double *[a];
    m[0] = new double[a * b];
    for (size_t i = 1; i != a; ++i)
        m[i] = m[i - 1] + b;
    return m;
}

void free_array (double ** m)
{
    delete [] m[0];
    delete [] m;
}
