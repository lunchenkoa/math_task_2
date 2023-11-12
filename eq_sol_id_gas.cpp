#include <iostream>
#include <string>
#include <cmath>
#include <algorithm>

#include "de_allocate.hpp"

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
    double ** u_L = create_array(array_length, 3);
    double ** u_R = create_array(array_length, 3);
    double ** F_L = create_array(array_length, 3);
    double ** F_R = create_array(array_length, 3);

    double ** F = create_array(array_length, 3);


    double * RHO = create_vector(array_length);
    double * V = create_vector(array_length);
    double * P = create_vector(array_length);

    // int counter = 0;
    double t, dt = 0.0;

    while (t <= time)
    {
        double D_L = speed_estimates(adiabat, u_init)[0];
        double D_R = speed_estimates(adiabat, u_init)[1];

        dt = Courant * dx / max(abs(D_L), abs(D_R));
        t += dt;

        for (size_t j = 1; j < array_length - 1; ++j)
        {
            for (int k = 0; k < 3; ++k)
            {
                u_L[j][k] = (u[j][k] + u[j - 1][k]) / 2;
                u_R[j][k] = (u[j + 1][k] + u[j][k]) / 2;
                F_R[j][k] = ;       // по презе F_L = F(u_L), F_R = F(u_R) что бы это не значило
                F_L[j][k] = ;

                if (D_L >= 0)
                {
                    F[j][k] = F_L[j][k];
                }
                else if (D_L <= 0 && 0 <= D_R)
                {
                    F[j][k] = (-D_L * F_R[j][k] + D_R * F_L[j][k] + D_L * D_R * (u_R[j][k] - u_L[j][k])) / (D_R - D_L);
                }
                else if (D_R <= 0)
                {
                    F[j][k] = F_R[j][k];
                }
            }
        }

        for (size_t j = 1; j < array_length - 1; ++j)
        {
            for (int k = 0; k < 3; ++k)
            {
                u[j][k] = ;
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
