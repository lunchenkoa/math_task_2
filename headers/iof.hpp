#include <string>

#ifndef IOF_INCLUDED
#define IOF_INCLUDED

bool initialization (std::string test, double& dens_L, double& vel_L, double& pres_L, \
                                       double& dens_R, double& vel_R, double& pres_R);
// bool save_results (std::string test, int N, double* x, double* pressure, double* density, double* velocity);

#endif