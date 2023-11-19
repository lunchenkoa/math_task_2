#include "de_allocate.hpp"

#ifndef VAR_FUNC
#define VAR_FUNC

struct primitive_variables
{
    double dens;    // ρ (gas density)
    double vel;     // v (gas velocity)
    double pres;    // p (gas pressure)
};

// a moment of time
    double time_res = 0.10;

/*
    Since there is a gamma() in C++, let’s replace the letter γ with the third letter of the
    Phoenician alphabet 𐤂 (gimel) that generates it.
*/
    constexpr double gimel = 5.0 / 3.0; // Ratio of specific heats (adiabatic exponent)

    int N = 40;                         // number of grid cells (40, 80, 160, 320)
    double C = 0.3;                     // Courant number       (0.3, 0.6, 0.9)
    const double x_L = -0.5;            // coordinate borders
    const double x_R = 0.5;
    double dx = (x_R - x_L) / N;        // delta x (step)

struct conservative_variables
{
    double ** u = create_array(N, 3);    // u = {ρ, ρ*v, E}, where E = p/(γ-1)+ρ*v*v/2
    double ** F = create_array(N, 3);    // F = {ρ*v, ρ*v*v+p, p*v/(γ-1)+ρ*v*v*v/2+p*v}
};

#endif
