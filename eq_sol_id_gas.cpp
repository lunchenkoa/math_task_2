#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

#include "de_allocate.hpp"

using namespace std;

// const double tol = 1e-6;

// а нужна ли эта функция вообще? я ничего не понимаю...
void set_initial_values (double* DENS, double* VEL, double* PRES, int zero_index, double** arr) // expecting vectors RHO, V, P and N_0, u0
{
    int array_length = sizeof(DENS) / sizeof(DENS[0]);

    for (int i = 0; i <= zero_index; ++i)
    {
        DENS[i] = arr[0][0];
        VEL[i] = arr[1][0];
        PRES[i] = arr[2][0];
    }

    for (int i = zero_index + 1; i < array_length; ++i)
    {
        DENS[i] = arr[0][1];
        VEL[i] = arr[1][1];
        PRES[i] = arr[2][1];
    }
}

void real2u_or_F (double* DENS, double* VEL, double* PRES, double** arr, double adiabat, bool is_arr_u)
{
    int array_length = sizeof(arr) / sizeof(arr[0]);

    if (is_arr_u)
    {
        for (int i = 0; i < array_length; ++i)
        {
            arr[i][0] = DENS[i];
            arr[i][1] = DENS[i] * VEL[i];
            arr[i][2] = PRES[i] / (adiabat - 1) + DENS[i] * pow(VEL[i], 2) / 2;
        }
    }
    else
    {
        for (int i = 0; i < array_length; ++i)
        {
            arr[i][0] = DENS[i] * VEL[i];
            arr[i][1] = DENS[i] * pow(VEL[i], 2) + PRES[i];
            arr[i][2] = VEL[i] * (PRES[i] / (adiabat - 1) + PRES[i] + DENS[i] * pow(VEL[i], 2) / 2);
        }
    }
}

void u2real (double* DENS, double* VEL, double* PRES, double** arr, double adiabat) // expect arr = u
{
    int array_length = sizeof(arr) / sizeof(arr[0]);

    for (size_t i = 1; i <= array_length; ++i)
    {
        DENS[i] = arr[i][1];
    }

    for (size_t i = 0; i < array_length; ++i)
    {
        VEL[i] = arr[i][2] / DENS[i];
        PRES[i] = (arr[i][3] - DENS[i] * pow(VEL[i], 2) / 2) * (adiabat - 1);
    }
}

// Функция для вычисления скорости звука
double sound_speed (double DENS, double PRES, double adiabat)
{
    return sqrt(adiabat * PRES / DENS);
}

double vmax (double** arr, double adiabat) // expect arr = u
{
    double v_max = -1.0;
    int array_length = sizeof(arr) / sizeof(arr[0]);

    double * RHO1 = create_vector(array_length);
    double * V1 = create_vector(array_length);
    double * P1 = create_vector(array_length);

    u2real(RHO1, V1, P1, arr, adiabat);

    for (size_t i = 1; i <= array_length - 2; ++i)
    {
        double c = sound_speed(RHO1[i], P1[i], adiabat);
        v_max = max(v_max, abs(V1[i]) + c);
    }

    free_vector(RHO1);
    free_vector(V1);
    free_vector(P1);

    return v_max;
}

// // Функция для вычисления максимальной скорости
// double* speed_estimates(double* u_l, double* u_r) 
// {   
//     /* 
//     Там в формуле на слайде 56 под максимумом стоит какое-то а...
//     Я так и не пон, что это значит. Типо, скорости D_l и D_r - трехмерные???
//     */
//     double* D = new double[2];   
//     D[0] = min(fabs(matrix_eigenvalue(u_l)), fabs(matrix_eigenvalue(u_r))); 
//     D[1] = max(fabs(matrix_eigenvalue(u_l)), fabs(matrix_eigenvalue(u_r))); 
    
//     return D;
// }

void godunov (double** u, double** F, int64_t number_of_nods, double time, double Courant, double dx, double adiabat)
{
    double ** temp_u = create_array(3, number_of_nods + 2);
    double ** F_l = create_array(3, number_of_nods + 2);
    double ** F_r = create_array(3, number_of_nods + 2);

    double * RHO = create_vector(number_of_nods + 2);
    double * V = create_vector(number_of_nods + 2);
    double * P = create_vector(number_of_nods + 2);

    int counter = 0;
    // double D_l, D_r;

    double t, dt = 0.0;
    double v_max = 0;

    while (t <= time)
    {
        v_max = vmax(u, adiabat);
        dt = Courant * dx / v_max;
        t += dt;

        for (size_t j = 1; j <= number_of_nods; ++j)
        {
            for (size_t k = 0; k < 3; ++k)
            {
                F_r[j][k] = (F[j][k] + F[j + 1][k]) / 2.0 - v_max / 2.0 * (u[j + 1][k] - u[j][k]);
                F_l[j][k] = (F[j - 1][k] + F[j][k]) / 2.0 - v_max / 2.0 * (u[j][k] - u[j - 1][k]);
            }
        }

        for (size_t j = 1; j <= number_of_nods; ++j)
        {
            for (int k = 0; k < 3; ++k)
            {
                temp_u[j][k] = u[j][k] - dt / dx * (F_r[j][k] - F_l[j][k]);
            }
        }

        for (size_t j = 1; j <= number_of_nods; ++j)
        {
            for (int k = 0; k < 3; ++k)
            {
                u[j][k] = temp_u[j][k];
            }
        }

        for (int k = 0; k < 3; ++k)
            {
                u[0][k] = u[1][k];
                u[number_of_nods + 1][k] = u[number_of_nods][k];
            }

        u2real(RHO, V, P, u, adiabat);
        real2u_or_F (RHO, V, P, F, adiabat, false);

        ++counter;
    }
    
    free_array(temp_u);
    free_array(F_l);
    free_array(F_r);

    free_vector(RHO);
    free_vector(V);
    free_vector(P);
}


// double* matrix_eigenvalue (double* u, double density, double velocity, double pressure)
// {
//     double *lambda = new double[3];
//     lambda[0] = velocity - sound_speed(density, pressure);
//     lambda[1] = velocity + sound_speed(density, pressure);
//     lambda[2] = velocity;
    
//     return lambda;
// }

// // Функция для вычисления потоков через грани ячейки
// void compute_fluxes(double density_L, double velocity_L, double pressure_L,
//                     double density_R, double velocity_R, double pressure_R,
//                     double &flux_density, double &flux_momentum, double &flux_energy) 
// {
//     double sound_speed_L = sound_speed(density_L, pressure_L);
//     double sound_speed_R = sound_speed(density_R, pressure_R);
    
//     double p_star = 0.5 * (pressure_L + pressure_R - 
//                     density_L * velocity_L * sound_speed_L - 
//                     density_R * velocity_R * sound_speed_R) / 
//                     (0.5 * (sound_speed_L + sound_speed_R));
    
//     if (p_star <= 0.0) 
//     {
//         p_star = tol;
//     }
    
//     double v_star = 0.5 * (velocity_L + velocity_R + 
//                     (pressure_L - pressure_R) / 
//                     (density_L * sound_speed_L + density_R * sound_speed_R));
    
//     double rho_star = 0.5 * (density_L + density_R);
    
//     if (v_star >= 0.0) 
//     {
//         flux_density = density_L * velocity_L;
//         flux_momentum = density_L * velocity_L * velocity_L + pressure_L;
//         flux_energy = velocity_L * (density_L * (pressure_L + density_L * velocity_L * velocity_L) / (density_L * pressure_L));
//     } else 
//     {
//         flux_density = density_R * velocity_R;
//         flux_momentum = density_R * velocity_R * velocity_R + pressure_R;
//         flux_energy = velocity_R * (density_R * (pressure_R + density_R * velocity_R * velocity_R) / (density_R * pressure_R));
//     }
    
//     if (v_star > 0.0) 
//     {
//         flux_density += p_star - pressure_L;
//         flux_momentum += (p_star + density_L * v_star * v_star) * v_star - density_L * velocity_L * (v_star - velocity_L);
//         flux_energy += v_star * (p_star + density_L * v_star * v_star);
//     } else 
//     {
//         flux_density += p_star - pressure_R;
//         flux_momentum += (p_star + density_R * v_star * v_star) * v_star - density_R * velocity_R * (v_star - velocity_R);
//         flux_energy += v_star * (p_star + density_R * v_star * v_star);
//     }
// }
