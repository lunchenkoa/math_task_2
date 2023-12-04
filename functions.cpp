#include <iostream>
#include <string>
#include <cmath>
#include <algorithm>

#include "headers/iof.hpp"
#include "headers/variables.hpp"
#include "headers/de_allocate.hpp"

using namespace std;

void initialization_of_IC (double * x, int N, primitive_variables* states, \
                                  primitive_variables& left, primitive_variables& right)
{
    for (size_t i = 0; i < N; ++i)
    {
        if (x[i] < 0)
        {
            states[i].dens = left.dens;
            states[i].vel = left.vel;
            states[i].pres = left.pres;
        }
        else
        {
            states[i].dens = right.dens;
            states[i].vel = right.vel;
            states[i].pres = right.pres;
        }
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

void compute_sound_speed (int N, double adiabat, primitive_variables* prims, double* speed_res)
{
    for (size_t i = 0; i < N; ++i)
        speed_res[i] = sqrt(adiabat * prims[i].pres / prims[i].dens);
}

double compute_max_velocity (int N, double adiabat, primitive_variables* prims, double* speed_res)
{
    compute_sound_speed (N, adiabat, prims, speed_res);

    double * speed_each_cell = create_vector(N);

    for (size_t j = 0; j < N; ++j)
    {
        speed_each_cell[j] = abs(prims[j].vel) + speed_res[j];
    }

    double max_vel = *max_element(speed_each_cell, speed_each_cell + N);

    free_vector(speed_each_cell);

    return max_vel;
}

void HLL_method (int N, double adiabat, conservative_variables cons, double Courant)
{
    double t, dt = 0.0;
    double v_max = 0;

    double ** F_star = create_array(N + 1, 3);
    double ** tmp_u = create_array(N, 3);

    for (size_t k = 0; k < 3; k++)
    {
        tmp_u[0][k] = cons.u[0][k];
        tmp_u[N - 1][k] = cons.u[N - 1][k];
    }

    double * D_L = create_vector(N + 1);
    double * D_R = create_vector(N + 1);
    double * s_vel = create_vector(N);
    double * F_L = create_vector(3);
    double * F_R = create_vector(3);

    primitive_variables * prims = new primitive_variables[N];
    int count = 0;

    while (t <= time_res)
    {   
        // if (count == 0)
        // {
        //     break;
        // }
        // for (size_t k = 0; k < 3; ++k)
        // {
        //     cons.u[0][k] = tmp_u[0][k];

        //     for (size_t j = 1; j < N - 1; ++j)
        //         cons.u[j][k] = tmp_u[j - 1][k];

        //     cons.u[N - 1][k] = tmp_u[N - 1][k];
        // }
        count +=1;
        // cout << "Time:" << t << endl;
        cons2prim (N, cons, prims, adiabat);

        compute_sound_speed (N, adiabat, prims, s_vel);
        v_max = compute_max_velocity (N, adiabat, prims, s_vel);
        dt = Courant * dx / v_max;
        
        for (size_t i = 1; i < N; ++i)
        {   
            // cout << " i = " << i <<" vel =" << prims[i].vel << ",  s_vel  = " << s_vel[i] << endl;
            D_L[i] = min(prims[i - 1].vel, prims[i].vel) - max(s_vel[i - 1], s_vel[i]);
            D_R[i] = max(prims[i - 1].vel, prims[i].vel) + max(s_vel[i - 1], s_vel[i]);
            // cout << "D_L = " << D_L[i] << " D_R = " << D_R[i] << endl;
            // cout << " i = " << i <<" D_L =" << D_L[i] << ",  D_R  = " << D_R[i] << endl;
        
        }

        // D_L[N] = s_vel[N - 1];            // kinda v_N = p_N = ρ_N = 0
        // D_R[N] = prims[N - 1].vel + s_vel[N - 1];

        // F_R[0] = prims[0].dens * prims[0].vel;
        // F_R[1] = prims[0].dens * pow(prims[0].vel, 2) + prims[0].pres;
        // F_R[2] = prims[0].pres * prims[0].vel / (adiabat - 1) + prims[0].dens * \
        //          pow(prims[0].vel, 3) / 2 + prims[0].pres * prims[0].vel;
            
        // for (size_t k = 0; k < 3; ++k)
        // {
        //     F_L[k] = 0;

        //     if (D_L[0] > 0)
        //     {
        //         F_star[0][k] = F_L[k];
        //     }
        //     else if (D_L[0] <= 0 && 0 <= D_R[0])
        //     {
        //         F_star[0][k] = (-D_L[0] * F_R[k] + D_R[0] * F_L[k] + D_L[0] * D_R[0] * \
        //                        (0 - cons.u[1][k])) / (D_R[0] - D_L[0]);
        //     }
        //     else if (D_R[0] < 0)
        //     {
        //         F_star[0][k] = F_R[k];
        //     }
        // }

        // F_L[0] = prims[N - 1].dens * prims[N - 1].vel;
        // F_L[1] = prims[N - 1].dens * pow(prims[N - 1].vel, 2) + prims[N - 1].pres;
        // F_L[2] = prims[N - 1].pres * prims[N - 1].vel / (adiabat - 1) + prims[N - 1].dens * \
        //          pow(prims[N - 1].vel, 3) / 2 + prims[N - 1].pres * prims[N - 1].vel;
                      
        // for (size_t k = 0; k < 3; ++k)
        // {
        //     F_R[k] = 0;

        //     if (D_L[N] > 0)
        //     {
        //         F_star[N][k] = F_L[k];
        //     }
        //     else if (D_L[N] <= 0 && 0 <= D_R[N])
        //     {
        //         F_star[N][k] = (-D_L[N] * F_R[k] + D_R[N] * F_L[k] + D_L[N] * D_R[N] * \
        //                        (cons.u[N - 2][k])) / (D_R[N] - D_L[N]);
        //     }
        //     else if (D_R[N] < 0)
        //     {
        //         F_star[N][k] = F_R[k];
        //     }
        // }
        
        // Надо что-то сделать.....

        for (size_t j = 1; j < N; ++j)
        {
            F_L[0] = prims[j - 1].dens * prims[j - 1].vel;
            F_L[1] = prims[j - 1].dens * pow(prims[j - 1].vel, 2) + prims[j - 1].pres;
            F_L[2] = prims[j - 1].pres * prims[j - 1].vel / (adiabat - 1) + prims[j - 1].dens * \
                     pow(prims[j - 1].vel, 3) / 2 + prims[j - 1].pres * prims[j - 1].vel;
                     
            F_R[0] = prims[j].dens * prims[j].vel;
            F_R[1] = prims[j].dens * pow(prims[j].vel, 2) + prims[j].pres;
            F_R[2] = prims[j].pres * prims[j].vel / (adiabat - 1) + prims[j].dens * \
                     pow(prims[j].vel, 3) / 2 + prims[j].pres * prims[j].vel;
            
            for (size_t k = 0; k < 3; ++k)
            {
                if (D_L[j] > 0)
                {
                    F_star[j][k] = F_L[k];
                }
                else if (D_L[j] <= 0 && 0 <= D_R[j])
                {
                    F_star[j][k] = (-D_L[j] * F_R[k] + D_R[j] * F_L[k] + D_L[j] * D_R[j] * \
                                   (cons.u[j][k] - cons.u[j - 1][k])) / (D_R[j] - D_L[j]);
                }
                else if (D_R[j] < 0)
                {
                    F_star[j][k] = F_R[k];
                }
            }
            // cout << " j = " << j <<" F_star =" << F_star[j][0] << ",  " << F_star[j][1] << ",  " << F_star[j][2] << endl;
        }
        
        for (size_t j = 1; j < N - 1; ++j)
            for (size_t k = 0; k < 3; ++k)
                tmp_u[j][k] = cons.u[j][k] - dt / dx * (F_star[j + 1][k] - F_star[j][k]);

        for (size_t k = 0; k < 3; ++k)
        {
            cons.u[0][k] = tmp_u[0][k];

            for (size_t j = 1; j < N - 1; ++j)
                cons.u[j][k] = tmp_u[j - 1][k];

            cons.u[N - 1][k] = tmp_u[N - 1][k];
        }
        
        t += dt;
        // if (count == 1)
        //     break;
    }
    
    free_array(F_star);
    free_array(tmp_u);
    free_vector(D_L);
    free_vector(D_R);
    free_vector(s_vel);
    free_vector(F_L);
    free_vector(F_R);
    delete [] prims;
}
