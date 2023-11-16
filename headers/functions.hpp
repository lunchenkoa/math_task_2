// #include <string>
#include "variables.hpp"

#ifndef GOD_FUNC
#define GOD_FUNC

void initialization_of_IC (int N, primitive_variables* states, \
                                  primitive_variables& left, primitive_variables& right);
void prim2cons (int N, conservative_variables cons, primitive_variables* prims, double adiabat);
void cons2prim (int N, conservative_variables cons, primitive_variables* prims, double adiabat);
// double compute_sound_speed (primitive_variables& states, double adiabat);
// void compute_wave_speeds (int N, double adiabat, conservative_variables cons, \
                          primitive_variables* prims, double* D_L, double* D_R);
void HLL_method (int N, double adiabat, conservative_variables cons, primitive_variables* prims, \
                 double Courant);

#endif