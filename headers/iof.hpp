#include <string>
#include <vector>
#include "variables.hpp"

using namespace std;

#ifndef IOF_INCLUDED
#define IOF_INCLUDED

bool initialization (string& test, vector<double>& left, vector<double>& right);
bool save_results (string test, vector<double> x,  vector<double> dens, vector<double> vel, vector<double> pres, int N);

#endif