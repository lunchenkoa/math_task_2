// #include <string>
#include "variables.hpp"
#include <vector>
#include <Eigen/Dense>

#ifndef GOD_FUNC
#define GOD_FUNC

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

void initialization_of_IC (vector<double> x, int N, vector<double>& dens, vector<double>& vel, vector<double>& pres, \
  vector<double>& left, vector<double>& right);
void prim2cons (int N, vector<vector<double>>& u, vector<vector<double>>& F, vector<double> dens, vector<double> vel, vector<double> pres, double adiabat);
void cons2prim (int N, vector<vector<double>> u, vector<vector<double>> F, vector<double>& dens, vector<double>& vel, vector<double>& pres, double adiabat);
void compute_sound_speed (int N, double adiabat, vector<double> dens, vector<double> pres, vector<double>& speed_res);
double compute_max_velocity (int N, double adiabat, vector<double> dens, vector<double> vel, vector<double> pres, vector<double> speed_res);
// void compute_wave_speeds (int N, double adiabat, conservative_variables cons, \
                          primitive_variables* prims, double* D_L, double* D_R);
void HLL_method (int N, double adiabat, vector<vector<double>> u, vector<vector<double>> F, double Courant);

#endif