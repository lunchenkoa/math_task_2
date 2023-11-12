#include <iostream>
#include <string>
#include <cmath>
#include <algorithm>

#include "de_allocate.hpp"
#include "iof.hpp"

using namespace std;


void set_initial_values (double* DENS, double* VEL, double* PRES, int LR_sep, double** arr, int array_length) // expecting vectors RHO, V, P and N_0, u0
{
    for (size_t i = 0; i < LR_sep; ++i) // for L-part
    {
        DENS[i] = arr[0][0];
        VEL[i] = arr[1][0];
        PRES[i] = arr[2][0];
    }

    for (size_t i = LR_sep; i < array_length; ++i) // for R-part
    {
        DENS[i] = arr[0][1];
        VEL[i] = arr[1][1];
        PRES[i] = arr[2][1];
    }
}

// conversion of individual gas features into vector ones, namely u and F
void feats2vectors (double* DENS, double* VEL, double* PRES, double** arr, double adiabat, bool is_arr_u, int array_length)
{
    if (is_arr_u)
    {
        for (int i = 0; i < array_length; ++i)
        {
            arr[i][0] = DENS[i];
            arr[i][1] = DENS[i] * VEL[i];
            // equation closing the system:
            //     p = ρε(γ-1), i.e. ε = p / (ρ(γ-1)),
            // therefore, the 3d component of the vector u:
            //     ρ(ε + v^2 / 2) = p / (γ-1) + ρv^2 / 2.
            arr[i][2] = PRES[i] / (adiabat - 1) + DENS[i] * pow(VEL[i], 2) / 2;
        }
    }
    else
    {
        for (int i = 0; i < array_length; ++i)
        {
            arr[i][0] = DENS[i] * VEL[i];
            arr[i][1] = DENS[i] * pow(VEL[i], 2) + PRES[i];
            // similarly, the 3d component of the vector F:
            //     ρv(ε + v^2 / 2 + p / ρ) = pv / (γ-1) + ρv^3 / 2 + pv.
            arr[i][2] = PRES[i] * VEL[i] / (adiabat - 1) + DENS[i] * pow(VEL[i], 3) / 2 + PRES[i] * VEL[i];
        }
    }
}

// inverse transformation of vector characteristics of gas into individual features, namely rho, v, p
void vectors2feats (double* DENS, double* VEL, double* PRES, double** arr, double adiabat, int array_length)
{   // expect arr = u
    for (size_t i = 0; i < array_length; ++i)
    {
        DENS[i] = arr[i][0];
        VEL[i] = arr[i][1] / DENS[i];
        PRES[i] = (arr[i][2] - DENS[i] * pow(VEL[i], 2) / 2) * (adiabat - 1);
    }
}

double sound_speed (double density, double pressure, double adiabat)
{
    return sqrt(adiabat * pressure / density);
}

double* speed_estimates(double adiabat, double** u_init) 
{   
    double sound_vel_L = sound_speed(u_init[0][0], u_init[0][2], adiabat);
    double sound_vel_R = sound_speed(u_init[1][0], u_init[1][2], adiabat);

    // determine the propagation speed
    double * D = create_vector(2);   
    D[0] = min(u_init[0][1], u_init[1][1]) - max(sound_vel_L, sound_vel_R); 
    D[1] = max(u_init[0][1], u_init[1][1]) + max(sound_vel_L, sound_vel_R); 

    // free_vector(D);       Разве после этого в return не будет просто ничего, тк ты его удалила а потом передала на выход...?
    
    return D;
}

void HLL_method (double** u, double** u_init, double** F, double time, double Courant, double dx, double adiabat, int array_length)
{
    double u_L, u_R, F_L, F_R, F_C, F_star;
    
    double ** step_u = create_array(array_length + 1, 3);
    double * RHO = create_vector(array_length + 1);
    double * V = create_vector(array_length + 1);
    double * P = create_vector(array_length + 1);

    double D_L = speed_estimates(adiabat, u_init)[0];
    double D_R = speed_estimates(adiabat, u_init)[1];

    double t, dt = 0.0;

    while (t <= time)
    {
        dt = Courant * dx / max(abs(D_L), abs(D_R));
        t += dt;

        for (size_t j = 1; j < array_length - 1; ++j)
        {
            for (int k = 0; k < 3; ++k)
            {   
                u_L = u[j][k];
                u_R = u[j + 1][k];
                F_L = (F[j - 1][k] + F[j][k]) / 2 - D_L / 2 * (u[j][k] - u[j-1][k]);
                F_R = (F[j + 1][k] + F[j][k]) / 2 - D_R / 2 * (u[j + 1][k] - u[j][k]);
                F_C = (-D_L * F_R + D_R * F_L + D_L * D_R * (u_R - u_L)) / (D_R - D_L);
                

                if (D_L >= 0)
                {
                    F_star = F_L;
                }
                else if (D_L <= 0 && 0 <= D_R)
                {
                    F_star = F_C;
                }
                else if (D_R <= 0)
                {
                    F_star = F_R;
                }
                u[j][k + 1]
            }

        }

        for (size_t j = 1; j < array_length - 1; ++j)
        {
            for (int k = 0; k < 3; ++k)
            {
                u[j][k] = u[j][k] - dt/dx * (F_r(j, :) - F_l(j, :));
            }
        }

        for (int k = 0; k < 3; ++k)
            {
                u[0][k] = u[1][k];
                u[array_length][k] = u[array_length - 1][k]; // мб тут ошибка...
            }

        vectors2feats(RHO, V, P, u, adiabat, array_length);
        feats2vectors (RHO, V, P, F, adiabat, false, array_length);

        // ++counter;
    }
    
    free_array(u_L);
    free_array(u_R);
    free_array(F_L);
    free_array(F_R);

    free_vector(RHO);
    free_vector(V);
    free_vector(P);
}

