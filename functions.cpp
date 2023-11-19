#include <iostream>
#include <string>
#include <cmath>
#include <algorithm>

#include "headers/iof.hpp"
#include "headers/variables.hpp"

using namespace std;

void initialization_of_IC (int N, primitive_variables* states, \
                                  primitive_variables& left, primitive_variables& right)
{
    for (size_t i = 0; i < N / 2 + 1; ++i)
    {
        states[i].dens = left.dens;
        states[i].vel = left.vel;
        states[i].pres = left.pres;
    }
    for (size_t i = N / 2 + 1; i < N; ++i)
    {
        states[i].dens = right.dens;
        states[i].vel = right.vel;
        states[i].pres = right.pres;
    }
}

void prim2cons (int N, conservative_variables cons, primitive_variables* prims, double adiabat)
{
    for (size_t i = 0; i < N; ++i)
    {
        cons.u[i][0] = prims[i].dens;
        cons.u[i][1] = prims[i].dens * prims[i].vel;
        // equation closing the system:
        //     p = ρε(γ-1), i.e. ε = p / (ρ(γ-1)),
        // therefore, the 3d component of the vector u:
        //     ρ(ε + v^2 / 2) = p / (γ-1) + ρv^2 / 2.
        cons.u[i][2] = prims[i].pres / (adiabat - 1) + prims[i].dens * pow(prims[i].vel, 2) / 2;
    }
    for (size_t i = 0; i < N; ++i)
    {
        cons.F[i][0] = prims[i].dens * prims[i].vel;
        cons.F[i][1] = prims[i].dens * pow(prims[i].vel, 2) + prims[i].pres;
        // similarly, the 3d component of the vector F:
        //     ρv(ε + v^2 / 2 + p / ρ) = pv / (γ-1) + ρv^3 / 2 + pv.
        cons.F[i][2] = prims[i].pres * prims[i].vel / (adiabat - 1) + prims[i].dens * \
                       pow(prims[i].vel, 3) / 2 + prims[i].pres * prims[i].vel;
    }
}

void cons2prim (int N, conservative_variables cons, primitive_variables* prims, double adiabat)
{
    for (size_t i = 0; i < N; ++i)
    {
        prims[i].dens = cons.u[i][0];
        prims[i].vel = cons.u[i][1] / prims[i].dens;
        prims[i].pres = (cons.u[i][2] - prims[i].dens * pow(prims[i].vel, 2) / 2) * (adiabat - 1);
    }
}

// double compute_sound_speed (primitive_variables& states, double adiabat)
// {
//     return sqrt(adiabat * states.pres / states.dens);
// }

// void compute_wave_speeds (int N, double adiabat, conservative_variables cons, primitive_variables* prims, double* D_L, double* D_R)
// {
//     cons2prim (N, cons, prims, adiabat);

//     double sound_vel_L, sound_vel_R = 0;

//     for (size_t i = 0; i <= N; ++i)
//     {
//         sound_vel_L = compute_sound_speed (prims[i], adiabat);
//         sound_vel_R = compute_sound_speed (prims[i + 1], adiabat);

//         D_L[i] = min(prims[i].vel, prims[i + 1].vel) - max(sound_vel_L, sound_vel_R);
//         D_R[i] = max(prims[i].vel, prims[i + 1].vel) + max(sound_vel_L, sound_vel_R);
//     }
// }

void HLL_method (int N, double adiabat, conservative_variables cons, double Courant)
{
    double t, dt = 0.0;
    // double u_L, u_R, F_C;

    double ** F_star = create_array(N + 1, 3);
    double ** tmp_u = create_array(N, 3);

    double * D_L = create_vector(N + 1);
    double * D_R = create_vector(N + 1);
    double * s_vel = create_vector(N);
    double * F_L = create_vector(3);
    double * F_R = create_vector(3);

    primitive_variables * prims = new primitive_variables[N];

    while (t <= time)
    {
        cons2prim (N, cons, prims, adiabat);

        for (size_t i = 0; i < N; ++i)
            s_vel[i] = adiabat * prims[i].pres / prims[i].dens;

        D_L[0] = -s_vel[0];               // kinda v_{-1} = p_{-1} = ρ_{-1} = 0
        D_R[0] = prims[0].vel + s_vel[0];
        
        for (size_t i = 1; i < N; ++i)
        {      
            D_L[i] = min(prims[i - 1].vel, prims[i].vel) - max(s_vel[i - 1], s_vel[i]);
            D_R[i] = max(prims[i - 1].vel, prims[i].vel) + max(s_vel[i - 1], s_vel[i]);
        }

        D_L[N] = s_vel[N - 1];            // kinda v_N = p_N = ρ_N = 0
        D_R[N] = prims[N - 1].vel + s_vel[N - 1];

        for (size_t j = 0; j <= N; ++j)
        {
            F_L[0] = (prims[j].dens * prims[j].vel) - (prims[j - 1].dens * prims[j - 1].vel);
            F_L[1]  = (prims[j].dens * pow(prims[j].vel, 2) + prims[j].pres) - (prims[j - 1].dens * pow(prims[j - 1].vel, 2) + prims[j - 1].pres) ;
            F_L[2] = (prims[j].pres * prims[j].vel / (adiabat - 1) + prims[j].dens * pow(prims[j].vel, 3) / 2 + prims[j].pres * prims[j].vel) -\
                     (prims[j - 1].pres * prims[j - 1].vel / (adiabat - 1) + prims[j - 1].dens * pow(prims[j - 1].vel, 3) / 2 + prims[j - 1].pres * prims[j - 1].vel);

            F_R[0] = (prims[j + 1].dens * prims[j + 1].vel) - (prims[j].dens * prims[j].vel);
            F_R[1]  = (prims[j + 1].dens * pow(prims[j + 1].vel, 2) + prims[j + 1].pres) - (prims[j].dens * pow(prims[j].vel, 2) + prims[j].pres) ;
            F_R[2] = (prims[j + 1].pres * prims[j + 1].vel / (adiabat - 1) + prims[j + 1].dens * pow(prims[j + 1].vel, 3) / 2 + prims[j + 1].pres * prims[j + 1].vel) -\
                     (prims[j].pres * prims[j].vel / (adiabat - 1) + prims[j].dens * pow(prims[j].vel, 3) / 2 + prims[j].pres * prims[j].vel);


            
            for (size_t k = 0; k < 3; ++k)
            {
                if (D_L[j] > 0)
                {
                    F_star[j][k] = F_L[k];
                }
                else if (D_L[j] <= 0 && 0 <= D_R[j])
                {
                    F_star[j][k] = (-D_L[j] * F_R[k] + D_R[j] * F_L[k] + D_L[j] * D_R[j] * (cons.u[j-1][k] - cons.u[j+1][k])) / (D_R[j] - D_L[j]);
                }
                else if (D_R[j] < 0)
                {
                    F_star[j][k] = F_R[k];
                }
            }
        }

        for (size_t j = 0; j < N; ++j)
        {
            for (size_t k = 0; k < 3; ++k)
            {
                tmp_u[j][k] = cons.u[j][k] - dt / dx * (F_star[j + 1][k] - F_star[j][k]);
            }
        }
        // Надо обновить граничные условня для всех индексов, кроме первого и последнего
        for (size_t j = 0; j < 3; ++j)
        {
            cons.u[0][j] = tmp_u[0][j];
            for (size_t i = 1; i <= N; ++i)
            {
                cons.u[i][j] = tmp_u[i - 1][j];
            }
        }
        
        dt = Courant * dx / max(abs(D_L[0]), abs(D_R[N]));
        t += dt;

    }
    // Перевод в примитивные переменные из u
    
    free_array(F_star);
    free_array(tmp_u);
    free_vector(D_L);
    free_vector(D_R);
    free_vector(s_vel);
    free_vector(F_L);
    free_vector(F_R);
    delete [] prims;
}
