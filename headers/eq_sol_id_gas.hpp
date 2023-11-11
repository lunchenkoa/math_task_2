#include <string>
#ifndef REIMAN_FUNC
#define REIMAN_FUNC

void initial_conditions(std::string test_numbr, double& ro_l,\
                        double& v_l,double& p_l,double& ro_r,\
                        double& v_r,double& p_r,double& t);

double sound_speed(double density, double pressure);

double max_speed(double density, double velocity, double pressure);

void compute_fluxes(double density_L, double velocity_L, double pressure_L,
                    double density_R, double velocity_R, double pressure_R,
                    double &flux_density, double &flux_momentum, double &flux_energy);

double der_fun_p(double x, double ro_l, double v_l,double p_l, double c_l,\
                      double ro_r, double v_r,double p_r, double c_r,\
                      double gamma, int state);

void solution_calculation (double ro_l, double v_l,double p_l, double c_l,\
                           double ro_r, double v_r,double p_r, double c_r,\
                           double p_star, double gamma, double t, int state,\
                           double *V, double *RO, double *P, double *x,\
                           int N, double xl, double xr, int& state_solution);

void configurationcalc(double ro_l, double v_l, double p_l, double c_l, double ro_r, double v_r, double p_r, double c_r, \
                        double p_star, double ro_star_l, double ro_star_r, double v_star, double x_star,\
                        double gamma, double t, int state,\
                        double x_rw_h, double x_rw_t, double x_lw_h, double x_lw_t, double x_d_l, double x_d_r, \
                        double x, double& V, double& RO, double& P);

#endif