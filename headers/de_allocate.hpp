// #include <cstdlib>

#ifndef DEALLOC_INCLUDED
#define DEALLOC_INCLUDED

double ** create_array (size_t rows, size_t cols);
double * create_vector (size_t a);
void free_array (double ** m);
void free_vector (double * v);

#endif
