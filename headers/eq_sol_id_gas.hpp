#include <string>
#ifndef GOD_FUNC
#define GOD_FUNC
void set_initial_values (double* DENS, double* VEL, double* PRES, int LR_sep, double** arr, int array_length);
void feats2vectors (double* DENS, double* VEL, double* PRES, double** arr, double adiabat, bool is_arr_u, int array_length);
void vectors2feats (double* DENS, double* VEL, double* PRES, double** arr, double adiabat, int array_length);
double sound_speed (double density, double pressure, double adiabat);
// double** speed_estimates(double adiabat, double** u_arr, int array_length);
void speed_estimates(double adiabat, double** u_arr, int array_length, double** D);
void HLL_method (double** u, double** u_init, double** F, double time, double Courant, double dx, double adiabat, int array_length);
#endif
