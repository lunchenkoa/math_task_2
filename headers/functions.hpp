// #include <string>
#include "variables.hpp"
#include <vector>
#include <Eigen/Dense>

#ifndef GOD_FUNC
#define GOD_FUNC

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

void initialization_of_IC (VectorXd x, int N, vector<primitive_variables> states, \
                                  primitive_variables& left, primitive_variables& right);
void prim2cons (int N, MatrixXd u, MatrixXd F, vector<primitive_variables> prims, double adiabat);
void cons2prim (int N, MatrixXd u, MatrixXd F, vector<primitive_variables> prims, double adiabat);
void compute_sound_speed (int N, double adiabat, vector<primitive_variables> prims, VectorXd speed_res);
double compute_max_velocity (int N, double adiabat, vector<primitive_variables> prims, VectorXd speed_res);
// void compute_wave_speeds (int N, double adiabat, conservative_variables cons, \
                          primitive_variables* prims, double* D_L, double* D_R);
void HLL_method (int N, double adiabat, MatrixXd u, MatrixXd F, double Courant);

#endif