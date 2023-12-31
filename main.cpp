#include <iostream>
#include <fstream>
#include <string>

#include "headers/variables.hpp"
#include "headers/iof.hpp"
#include "headers/functions.hpp"
#include "headers/de_allocate.hpp"

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
    
// Parameter initialization
//    {rho,v,p}_L = left.{dens,vel,pres}, {rho,v,p}_R = right.{dens,vel,pres}

    primitive_variables left, right;
    initialization (test_file, left, right);

// Allocation of memory to dynamic variables

    double * x = create_vector(N);
    // cout << "Init Values: " << endl;
    for (size_t i = 0; i < N; ++i)
    {
        x[i] = x_L + (i + 0.5) * dx;
    }

    primitive_variables * init_features = new primitive_variables[N]; // array for init {rho, v, p}
    initialization_of_IC(x, N, init_features, left, right);           // which is half filled with 
                                                                      // left characteristics and 
                                                                      // half with right ones

    conservative_variables cons_vars;
    cons_vars.u = create_array(N, 3);
    cons_vars.F = create_array(N, 3);

// Solution

    prim2cons (N, cons_vars, init_features, gimel); // fill the cons_vars

    // for (size_t i = 0; i < N; ++i)
    // {
    //     // cout << "i = " << i << " x = " << x[i] << " P = " << init_features[i].pres << " Rho = " << init_features[i].dens << " V = " << init_features[i].vel << endl;
    //     // cout << "i = " << i << " x = " << x[i] << " u = " << cons_vars.u[i][0] << ", " << cons_vars.u[i][1] << ", " << cons_vars.u[i][2] << ", " <<  " F = " <<  cons_vars.F[i][0] << ", " << cons_vars.F[i][1] << ", " << cons_vars.F[i][2] << endl;
    // }

    HLL_method (N, gimel, cons_vars, C);

    primitive_variables * final_features = new primitive_variables[N];
    cons2prim (N, cons_vars, final_features, gimel);

    save_results (test_nmbr, x, final_features, N);

// Deallocation of memory

    delete [] init_features;
    delete [] x;
    free_array(cons_vars.u);
    free_array(cons_vars.F);
    delete [] final_features;
    
    return 0;
}
