#include <vector>
#include "variables.hpp"

#ifndef GOD_FUNC
#define GOD_FUNC

using namespace std;

void initialization_of_IC (vector<double> x, int N, vector<double>& dens, vector<double>& vel, vector<double>& pres, vector<double>& left, vector<double>& right);
void prim2cons (int N, vector<vector<double>>& u, vector<vector<double>>& F, vector<double> dens, vector<double> vel, vector<double> pres, double adiabat);
void cons2prim (int N, vector<vector<double>> u, vector<vector<double>> F, vector<double>& dens, vector<double>& vel, vector<double>& pres, double adiabat);
void compute_sound_speed (int N, double adiabat, vector<double> dens, vector<double> pres, vector<double>& speed_res);
double compute_max_velocity (int N, double adiabat, vector<double> dens, vector<double> vel, vector<double> pres, vector<double> speed_res);
void compute_wave_speed (double vel_L, double vel_R, double s_vel_L, double s_vel_R, double& D_L, double& D_R);
void HLL_method (int N, double adiabat, vector<vector<double>>& u, vector<vector<double>>& F, double Courant);

#endif