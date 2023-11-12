#include <iostream>
#include <string>

#include "headers/de_allocate.hpp"
#include "headers/iof.hpp"
#include "headers/eq_sol_id_gas.hpp"

using namespace std;

/*
    Numerical solution of equations one-dimensional ideal gas dynamics with adiabatic exponent
    𝛾 = 5/3 using Harten-Lax-van Leer (HLL) method.
*/

int main ()
{
// Add test selection for single compilation

    string test = "";

    while (true)
    {
        cout << "Choose a test!\nPossible input options: 1, 2, 3\n";
        cin >> test;
        cout << '\n';

        test = "input" + test + ".txt";

        if (!((test == "input1.txt") || (test == "input2.txt") || (test == "input3.txt")))
        {
            cerr << "Invalid input!\n";
        }
        else
        {
            break;
        }
    }

// Parameter initialization:
//    rho_L (rho_R) = gas density on the left (on the right),
//    v_L   (v_R)   = gas velocity on the left (on the right),
//    p_L   (p_R)   = gas pressure on the left (on the right).

    double rho_L, v_L, p_L, rho_R, v_R, p_R;
    initialization(test, rho_L, v_L, p_L, rho_R, v_R, p_R);

    double ** u_0 = create_array(2, 3); // array with test's {(rho, v, p)_L, (rho, v, p)_R}
    u_0[0][0] = rho_L;
    u_0[0][1] = v_L;
    u_0[0][2] = p_L;
    u_0[1][0] = rho_R;
    u_0[1][1] = v_R;
    u_0[1][2] = p_R;

    // moment of time
    double time = 0.10;

/*
    Since there is a gamma() in C++, let’s replace the letter γ with the third letter of the
    Phoenician alphabet 𐤂 (gimel) that generates it.
*/
    constexpr double gimel = 5.0 / 3.0; // Ratio of specific heats (adiabatic exponent)

// Continued definition of parameters

    // selecting the number of grid cells (N) and the Courant number (C)
    int N = 40;      // 40, 80, 160, 320      // Хз, как для случая A считать число Куранта (!!!)
    double C = 0.3;  // 0.3, 0.6, 0.9

    // coordinate borders
    const double x_L = -0.5;
    const double x_R = 0.5;
    // delta x (step)
    double dx = (x_R - x_L) / N;
    // origin (here we can separate Left and Right) along the x axis (left-right separator)
    int LR_sep = N / 2 + 1;
    // number of grid nodes
    int nodes = N + 1;

/*  
  + хз, надо ли проверять начальные условия для определения конфинурации и на условие вакуума.

  - ну вообще я не видела разделения на конфигурации... по описанию HLL это две УВ, разлетающиеся
    в разные стороны, двум УВ соответствует кофигурация Б
*/

// Allocation of memory to dynamic variables

    double * x = create_vector(nodes); 
    double * RHO = create_vector(nodes);
    double * V = create_vector(nodes);
    double * P = create_vector(nodes);

    double ** u = create_array(nodes, 3);
    double ** F = create_array(nodes, 3);

    for (size_t i = 1; i < nodes; ++i)
        x[i] = x_L + (i - 1) * dx;

    set_initial_values (RHO, V, P, LR_sep, u_0, nodes);
    feats2vectors (RHO, V, P, u, gimel, true, nodes);
    feats2vectors (RHO, V, P, F, gimel, false, nodes);
    HLL_method (u, u_0, F, time, C, dx, gimel, nodes);
    vectors2feats (RHO, V, P, u, gimel, nodes);
    save_results (test, x, P, RHO, V);

    free_vector(x);
    free_vector(RHO);
    free_vector(V);
    free_vector(P);

    free_array(u_0);
    free_array(u);
    free_array(F);
    
    return 0;
}
