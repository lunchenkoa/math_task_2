#include <iostream>
#include <fstream>
#include <string>
#include <vector>
 
#include "headers/variables.hpp"
#include "headers/iof.hpp"
#include "headers/functions.hpp"

using namespace std;

int main ()
{
// Add test selection for single compilation

    string test_nmbr = "";
    string test_file = "";

    while (true)
    {
        cout << "Choose a test!\nPossible input options: 1, 2, 3\n";
        cin >> test_nmbr;
        cout << '\n';

        test_file = "input/input" + test_nmbr + ".txt";
        cout << test_file << endl;
        if (!((test_file == "input/input1.txt") || (test_file == "input/input2.txt") || \
              (test_file == "input/input3.txt")))
        {
            cerr << "Invalid input!\n";
        }
        else
        {
            break;
        }
    }
    
// Allocation of memory to dynamic variables

    vector<double> l_init, r_init;
    vector<double> x;
    vector<double> dens, vel, pres;
    vector<vector<double>> u, F;
    vector<double> f_dens, f_vel, f_pres;

// Parameter initialization

    l_init.resize(3);
    r_init.resize(3);

    initialization (test_file, l_init, r_init);

    x.resize(N);

    for (size_t i = 0; i < N; ++i)
    {
        x[i] = x_L + (i + 0.5) * dx;
    }

// Initialization of initial conditions (creating primitive variables)
    
    dens.resize(N);
     vel.resize(N);
    pres.resize(N);

    initialization_of_IC (x, N, dens, vel, pres, l_init, r_init);                          

// Creating conservative variables and filling them with values

    u.resize(3);
    F.resize(3);

    for (size_t var = 0; var < 3; ++var) // 3 rows, N cols
    {
        u[var].resize(N);
        F[var].resize(N);
    }

    prim2cons (N, u, F, dens, vel, pres, gimel);

// Solution

    HLL_method (N, gimel, u, F, C);

// Storing results in primitive variables

    f_dens.resize(N);
     f_vel.resize(N);
    f_pres.resize(N);

    cons2prim (N, u, F, f_dens, f_vel, f_pres, gimel);

    save_results (test_nmbr, x, f_dens, f_vel, f_pres, N);
    
    return 0;
}
