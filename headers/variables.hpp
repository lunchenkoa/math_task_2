#ifndef VAR_FUNC
#define VAR_FUNC

using namespace std;

    static double time_res = 0.1;       // a moment of time

/*
    Since there is a gamma() in C++, let’s replace the letter γ with the third letter of the
    Phoenician alphabet 𐤂 (gimel) that generates it.
*/
    constexpr double gimel = 5.0 / 3.0; // Ratio of specific heats (adiabatic exponent)

    static int N = 320;                  // number of grid cells (40, 80, 160, 320)
    static double C = 0.9;              // Courant number       (0.3, 0.6, 0.9)
    const double x_L = -0.5;            // coordinate borders
    const double x_R = 0.5;             //
    static double dx = (x_R - x_L) / N; // delta x (step)

#endif
