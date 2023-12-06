#include <string>
#include "variables.hpp"
#include <vector>
#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;


#ifndef IOF_INCLUDED
#define IOF_INCLUDED

bool initialization (string& test, vector<double>& left, vector<double>& right);
bool save_results (string test, VectorXd x,  vector<double> dens, vector<double> vel, vector<double> pres, int N);

#endif