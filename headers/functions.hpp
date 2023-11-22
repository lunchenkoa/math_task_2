// #include <string>
#include "variables.hpp"

#ifndef GOD_FUNC
#define GOD_FUNC

void initialization_of_IC (double * x, int N, primitive_variables* states, \
                                  primitive_variables& left, primitive_variables& right);
void prim2cons (int N, conservative_variables cons, primitive_variables* prims, double adiabat);
void cons2prim (int N, conservative_variables cons, primitive_variables* prims, double adiabat);
void compute_sound_speed (int N, double adiabat, primitive_variables* prims, double* speed_res);
double compute_max_velocity (int N, double adiabat, primitive_variables* prims, double* speed_res);
// void compute_wave_speeds (int N, double adiabat, conservative_variables cons, \
                          primitive_variables* prims, double* D_L, double* D_R);
void HLL_method (int N, double adiabat, conservative_variables cons, double Courant);

#endif