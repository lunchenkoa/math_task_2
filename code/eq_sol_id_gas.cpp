#include <iostream>
#include <string>
#include <cmath>
#include <algorithm>

#include "../headers/de_allocate.hpp"
#include "../headers/iof.hpp"

using namespace std;


// expecting vectors RHO, V, P and N_0, u0
void set_initial_values (double* DENS, double* VEL, double* PRES, int LR_sep, double** arr, int array_length) 
{
    for (size_t i = 0; i < LR_sep; ++i) // for L-part
    {
        DENS[i] = arr[0][0];
        VEL[i]  = arr[0][1];
        PRES[i] = arr[0][2];
    }

    for (size_t i = LR_sep; i < array_length; ++i) // for R-part
    {
        DENS[i] = arr[0][0];
        VEL[i]  = arr[0][1];
        PRES[i] = arr[0][2];
    }
}

// conversion of individual gas features into vector ones, namely u and F
void feats2vectors (double* DENS, double* VEL, double* PRES, double** arr, double adiabat, bool is_arr_u, int array_length)
{
    if (is_arr_u)
    {
        for (size_t i = 0; i < array_length; ++i)
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
        for (size_t i = 0; i < array_length; ++i)
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


double** speed_estimates(double adiabat, double** u_arr, int array_length, double** D)
{
    double * rho = create_vector(array_length);   //  new double [array_length];  // 
    double * v = create_vector(array_length);     //  new double [array_length];  //  
    double * p = create_vector(array_length);     //  new double [array_length];  //  
    // double ** D = create_array(2, array_length-1);

    vectors2feats(rho, v, p, u_arr, adiabat, array_length);

    double sound_vel_L, sound_vel_R = 0;
    for (size_t i = 0; i < array_length-1; ++i)
    {
        sound_vel_L = sound_speed(rho[i], p[i], adiabat);
        sound_vel_R = sound_speed(rho[i+1], p[i+1], adiabat);

        // determine the propagation speed
        D[0][i] = min(v[i], v[i + 1]) - max(sound_vel_L, sound_vel_R);
        D[1][i] = max(v[i], v[i + 1]) + max(sound_vel_L, sound_vel_R);
    }

    // delete [] p;
    // delete [] rho;
    // delete [] v;
    free_vector(rho);
    free_vector(v);
    free_vector(p);

    return D;
}

void HLL_method (double** u, double** u_init, double** F, double time, double Courant, double dx, double adiabat, int array_length)
{
    double u_L, u_R, F_L, F_R, F_C;
    
    double ** F_star = create_array(array_length - 1, 3);
    double ** step_u = create_array(array_length - 2, 3);
    double * RHO = create_vector(array_length);   // new double [array_length];  //   
    double * V = create_vector(array_length);   //   new double [array_length];  //   
    double * P = create_vector(array_length);   //   new double [array_length];  //   
    double ** D = create_array(2, array_length-1);

    double t, dt = 0.0;

    while (t <= time)
    {
        D = speed_estimates(adiabat, u, array_length, D);
        double * D_L = D[0];
        double * D_R = D[1];
        dt = Courant * dx / max(abs(D_L[0]), abs(D_R[array_length - 1]));
        t += dt;

        for (size_t j = 0; j < array_length-1; ++j)
        {
            for (size_t k = 0; k < 3; ++k)
            {   
                u_L = u[j][k];
                u_R = u[j + 1][k];
                F_L = F[j][k];
                F_R = F[j + 1][k];
                // F_L = (F[j - 1][k] + F[j][k]) / 2 - D_L[j-1] / 2 * (u[j][k] - u[j-1][k]);
                // F_R = (F[j + 1][k] + F[j][k]) / 2 - D_R[j-1] / 2 * (u[j + 1][k] - u[j][k]);
                F_C = (-D_L[j] * F_R + D_R[j] * F_L + D_L[j] * D_R[j] * (u_R - u_L)) / (D_R[j] - D_L[j]);
                

                if (D_L[j] >= 0)
                {
                    F_star[j][k] = F_L;
                }
                else if (D_L[j] <= 0 && 0 <= D_R[j])
                {
                    F_star[j][k] = F_C;
                }
                else if (D_R[j] <= 0)
                {
                    F_star[j][k] = F_R;
                }
            }
        }

        for (size_t j = 0; j < array_length - 2; ++j)
        {
            for (size_t k = 0; k < 3; ++k)
            {
                step_u[j][k] = u[j+1][k] - dt/dx * (F_star[j][k] - F_star[j+1][k]);
            }
        }

        for (size_t k = 0; k < 3; ++k)
            {
                u[0][k] = step_u[0][k];
                for (size_t g = 1; g <= array_length - 1; ++g)
                {
                    u[g][k] = step_u[g-1][k];
                }
                // u[array_length][k] = step_u[array_length - 2][k];
            }

        vectors2feats(RHO, V, P, u, adiabat, array_length);
        feats2vectors (RHO, V, P, F, adiabat, false, array_length);
        delete [] D_L;
        delete [] D_R;
        // free_vector(D_L);
        // free_vector(D_R);
        // free_array(D);

    }
    
    free_array(step_u);
    free_array(D);
    free_array(F_star);

    // delete [] P;
    // delete [] RHO;
    // delete [] V;
    free_vector(RHO);
    free_vector(V);
    free_vector(P);
}

