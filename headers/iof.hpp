#include <string>
#include "variables.hpp"

#ifndef IOF_INCLUDED
#define IOF_INCLUDED

bool initialization (std::string& test, primitive_variables& left, primitive_variables& right);
bool save_results (std::string test, double* x, primitive_variables* states);

#endif