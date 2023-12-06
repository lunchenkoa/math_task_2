#include "de_allocate.hpp"
#include <vector>
#include <Eigen/Dense>

#ifndef VAR_FUNC
#define VAR_FUNC

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;
// a moment of time
    static double time_res = 0.01;

/*
    Since there is a gamma() in C++, let’s replace the letter γ with the third letter of the
    Phoenician alphabet 𐤂 (gimel) that generates it.
*/
    constexpr double gimel = 5.0 / 3.0; // Ratio of specific heats (adiabatic exponent)

    static int N = 320;                         // number of grid cells (40, 80, 160, 320)
    static double C = 0.9;                     // Courant number       (0.3, 0.6, 0.9)
    const double x_L = -0.5;                   // coordinate borders
    const double x_R = 0.5;
    static double dx = (x_R - x_L) / N;        // delta x (step)

struct conservative_variables
{
    double** u;    // u = {ρ, ρ*v, E}, where E = p/(γ-1)+ρ*v*v/2
    double** F;    // F = {ρ*v, ρ*v*v+p, p*v/(γ-1)+ρ*v*v*v/2+p*v}

    // conservative_variables() :  u[create_array(N, 3)), F(create_array(N, 3)) {}

    // ~conservative_variables()
    // {
    //     for (int i = 0; i < N; ++i)
    //     {
    //         delete[] u[i];
    //         delete[] F[i];
    //     }
    //     delete[] u;
    //     delete[] F;
    // }
};


struct primitive_variables     // cnhernehf bp dtrnjhjd
{
    // VectorXd dens(N);
    // VectorXd vel(N);
    // VectorXd pres(N);
    // primitive_variables.dens.resize(N);
    // primitive_variables.vel.resize(N);
    // primitive_variables.pres.resize(N);
    
    // double dens;    // ρ (gas density)
    // double vel;     // v (gas velocity)
    // double pres;    // p (gas pressure)
};

#endif
