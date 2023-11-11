#include <string>
#ifndef GOD_FUNC
#define GOD_FUNC

void set_initial_values (double* PRES, double* DENS, double* VEL, int zero_index, double** arr);
void real2u_or_F (double* DENS, double* VEL, double* PRES, double** arr, double adiabat, bool is_arr_u);
void u2real (double* DENS, double* VEL, double* PRES, double** arr, double adiabat);
double sound_speed (double& DENS, double& PRES, double adiabat);
double vmax (double** arr, double adiabat);
void godunov (double** u, double** F, int64_t number_of_nods, double time, double Courant, double dx, double adiabat);

#endif
