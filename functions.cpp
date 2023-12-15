#include <iostream>
#include <string>
#include <cmath>
#include <algorithm>
#include <vector>

#include "headers/iof.hpp"
#include "headers/variables.hpp"
#include "headers/de_allocate.hpp"

using namespace std;

void initialization_of_IC (vector<double> x, int N, vector<double>& dens, vector<double>& vel, \
                           vector<double>& pres, vector<double>& left, vector<double>& right)
{
    for (size_t i = 0; i < N; ++i)
    {
        if (x[i] < 0)
        {
            dens[i] = left[0];
            vel[i] = left[1];
            pres[i] = left[2];
        }
        else
        {
            dens[i] = right[0];
            vel[i] = right[1];
            pres[i] = right[2];
        }
    }
}

void prim2cons (int N, vector<vector<double>>& u, vector<vector<double>>& F, vector<double> dens, \
                vector<double> vel, vector<double> pres, double adiabat)
{
    for (size_t i = 0; i < N; ++i)
    {
        u[0][i] = dens[i];
        u[1][i] = dens[i] * vel[i];
        // equation closing the system:
        //     p = ρε(γ-1), i.e. ε = p / (ρ(γ-1)),
        // therefore, the 3d component of the vector u:
        //     ρ(ε + v^2 / 2) = p / (γ-1) + ρv^2 / 2.
        u[2][i] = pres[i] / (adiabat - 1) + dens[i] * pow(vel[i], 2) / 2;
    }

    for (size_t i = 0; i < N; ++i)
    {
        F[0][i] = dens[i] * vel[i];
        F[1][i] = dens[i] * pow(vel[i], 2) + pres[i];
        // similarly, the 3d component of the vector F:
        //     ρv(ε + v^2 / 2 + p / ρ) = pv / (γ-1) + ρv^3 / 2 + pv.
        F[2][i] = pres[i] * vel[i] / (adiabat - 1) + dens[i] * pow(vel[i], 3) / 2 + pres[i] * vel[i];
    }
}

void cons2prim (int N, vector<vector<double>> u, vector<vector<double>> F, vector<double>& dens, \
                vector<double>& vel, vector<double>& pres, double adiabat)
{
    for (size_t i = 0; i < N; ++i)
    {
        dens[i] = u[0][i];
         vel[i] = u[1][i] / dens[i];
        pres[i] = (u[2][i] - dens[i] * pow(vel[i], 2) / 2) * (adiabat - 1);
    }
}

void compute_sound_speed (int N, double adiabat, vector<double> dens, vector<double> pres, \
                          vector<double>& speed_res)
{
    for (size_t i = 0; i < N; ++i)
        speed_res[i] = sqrt(adiabat * pres[i] / dens[i]);
}

double compute_max_velocity (int N, double adiabat, vector<double> dens, vector<double> vel, \
                             vector<double> pres, vector<double> speed_res)
{
    compute_sound_speed (N, adiabat, dens, pres, speed_res);

    vector<double> speed_each_cell;
    speed_each_cell.resize(N);

    for (size_t j = 0; j < N; ++j)
    {
        speed_each_cell[j] = abs(vel[j]) + speed_res[j];
    }

    double max_vel = *max_element(speed_each_cell.begin(), speed_each_cell.end());

    return max_vel;
}

void compute_wave_speed (double vel_L, double vel_R, double s_vel_L, double s_vel_R, \
                         double& D_L, double& D_R)
{
    D_L = min(vel_L, vel_R) - max(s_vel_L, s_vel_R);
    D_R = max(vel_L, vel_R) + max(s_vel_L, s_vel_R);
}

void HLL_method (int N, double adiabat, vector<vector<double>> u, vector<vector<double>> F, \
                 double Courant)
{
    double t = 0.0, dt = 0.0;
    double v_max = 0;

// Creating auxiliary arrays

    vector<vector<double>> tmp_u, F_star1, F_star2;

    tmp_u.resize(3);
    F_star1.resize(3);
    F_star2.resize(3);

    for (size_t var = 0; var < 3; ++var)
    {
        tmp_u[var].resize(N);       // the same size as the two-dimensional vector u[][]
        F_star1[var].resize(N + 1); // (N + 1) is the number of grid nodes
        F_star2[var].resize(N + 1); 
    }

// Setting boundary conditions (BC)

    for (size_t k = 0; k < 3; ++k)
    {
        tmp_u[k][0] = u[k][0];
        tmp_u[k][N - 1] = u[k][N - 1];
    }

// Creating vectors to store wave propagation speeds (D_L, D_R), sound speed (s_vel) and flows (F_L, F_R)
    
    vector<double> D_L, D_R, s_vel, F_L, F_R;

    D_L.resize(N + 1); 
    D_R.resize(N + 1);
    s_vel.resize(N);
    F_L.resize(3);
    F_R.resize(3);

// Creating primitive variables in the form of vectors for local use

    vector<double> dens, vel, pres;

    dens.resize(N);
     vel.resize(N);
    pres.resize(N);

// Main function loop

    // int count = 0;

    while (t <= time_res)
    {
        // count += 1;

    // The input of the function was u and F, transform them into primitive variables and place them
    // in vectors created here in advance
        cons2prim (N, u, F, dens, vel, pres, adiabat);

    // Calculate the sound velocity, and then the maximum speed on the grid to find the step dt
        compute_sound_speed (N, adiabat, dens, pres, s_vel);
        v_max = compute_max_velocity (N, adiabat, dens, vel, pres, s_vel);
        dt = Courant * dx / v_max;
    
    // Calculate the speed of wave propagation on the left and right
        for (size_t i = 1; i < N; ++i)
            compute_wave_speed (vel[i - 1], vel[i], s_vel[i - 1], s_vel[i], D_L[i], D_R[i]);

    // Loop to calculate left and right fluxes to calculate flux across a cell boundary
        for (size_t j = 1; j < N; ++j)
        {
            F_L[0] = dens[j - 1] * vel[j - 1];
            F_L[1] = dens[j - 1] * pow(vel[j - 1], 2) +  pres[j - 1];
            F_L[2] = pres[j - 1] * vel[j - 1] / (adiabat - 1) + dens[j - 1] * \
                     pow( vel[j - 1], 3) / 2 + pres[j - 1] * vel[j - 1];
               
            F_R[0] = dens[j] *  vel[j];
            F_R[1] = dens[j] * pow(vel[j], 2) +  pres[j];
            F_R[2] = pres[j] * vel[j] / (adiabat - 1) + dens[j] * pow(vel[j], 3) / 2 + pres[j] * vel[j];
            
            for (size_t k = 0; k < 3; ++k)
            {
                if (D_L[j] > 0) // ((D_L[j] >= 0) && (D_R[j] > 0))
                {
                    F_star1[k][j] = F_L[k];
                }
                else if ((D_L[j] <= 0) && (D_R[j] >= 0))
                {
                    F_star1[k][j] = (-D_L[j] * F_R[k] + D_R[j] * F_L[k] + D_L[j] * D_R[j] * \
                                    (u[k][j] -  u[k][j - 1])) / (D_R[j] - D_L[j]);
                }
                else if (D_R[j] < 0) // ((D_L[j] < 0) && (D_R[j] <= 0)) 
                {
                    F_star1[k][j] = F_R[k];
                }
            }
            cout << "F*_1[0][" << j << "] = " << F_star1[0][j] << ", F*_1[1][" << j << "] = " << F_star1[1][j] << ", F*_1[2][" << j << "] = " << F_star1[2][j] << '\n';
        }

    // Similar actions for the case when the left border is [i] and the right border is [i+1]
        for (size_t i = 1; i < N; ++i)
            compute_wave_speed (vel[i], vel[i + 1], s_vel[i], s_vel[i + 1], D_L[i], D_R[i]);

        for (size_t j = 1; j < N; ++j)
        {
            F_L[0] = dens[j] * vel[j];
            F_L[1] = dens[j] * pow(vel[j], 2) +  pres[j];
            F_L[2] = pres[j] * vel[j] / (adiabat - 1) + dens[j] * pow( vel[j], 3) / 2 + pres[j] * vel[j];
               
            F_R[0] = dens[j + 1] * vel[j + 1];
            F_R[1] = dens[j + 1] * pow(vel[j + 1], 2) +  pres[j + 1];
            F_R[2] = pres[j + 1] * vel[j + 1] / (adiabat - 1) + dens[j + 1] * \
                     pow(vel[j + 1], 3) / 2 + pres[j + 1] * vel[j + 1];
            
            for (size_t k = 0; k < 3; ++k)
            {
                if (D_L[j] > 0) // ((D_L[j] >= 0) && (D_R[j] > 0)) 
                {
                    F_star2[k][j] = F_L[k];
                }
                else if ((D_L[j] <= 0) && (D_R[j] >= 0))
                {
                    F_star2[k][j] = (-D_L[j] * F_R[k] + D_R[j] * F_L[k] + D_L[j] * D_R[j] * \
                                    (u[k][j + 1] -  u[k][j])) / (D_R[j] - D_L[j]);
                }
                else if (D_R[j] < 0) // ((D_L[j] < 0) && (D_R[j] <= 0))
                {
                    F_star2[k][j] = F_R[k];
                }
            }
            // cout << "F*_2[0][" << j << "] = " << F_star2[0][j] << ", F*_2[1][" << j << "] = " << F_star2[1][j] << ", F*_2[2][" << j << "] = " << F_star2[2][j] << '\n';
        }
        
    // Application of the Harten, Lax, van Leer scheme    
        for (size_t j = 1; j < N - 1; ++j)
            for (size_t k = 0; k < 3; ++k)
                tmp_u[k][j] =  u[k][j] - dt / dx * (F_star2[k][j] - F_star1[k][j]);

    // Updating BC
        for (size_t k = 0; k < 3; ++k)
        {
            u[k][0] = tmp_u[k][0];

            for (size_t j = 1; j < N - 1; ++j)
                u[k][j] = tmp_u[k][j - 1];

            u[k][N - 1] = tmp_u[k][N - 1];
        }

    // Updating the time counter
        t += dt;
    }
}
