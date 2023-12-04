#include <string>
#include "variables.hpp"
#include <vector>
#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;


#ifndef IOF_INCLUDED
#define IOF_INCLUDED

bool initialization (string& test, primitive_variables& left, primitive_variables& right);
bool save_results (string test, VectorXd x, vector<primitive_variables> states, int N);

#endif