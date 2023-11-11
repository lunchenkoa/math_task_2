#include <cstdlib>

#ifndef DEALLOC_INCLUDED
#define DEALLOC_INCLUDED

double ** create_array (size_t a, size_t b);
double * create_vector (size_t a);
void free_array (double ** m);
void free_vector (double * v);

#endif
