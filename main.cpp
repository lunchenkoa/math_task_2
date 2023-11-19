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

    primitive_variables * init_features = new primitive_variables[N]; // array for init {rho, v, p}
    initialization_of_IC(N, init_features, left, right);              // which is half filled with 
                                                                      // left characteristics and 
                                                                      // half with right ones

    double * x = create_vector(N);
    for (size_t i = 0; i < N; ++i)
        x[i] = x_L + i * dx;

    conservative_variables cons_vars;
    prim2cons (N, cons_vars, init_features, gimel); // здесь работа init_features по идее кончается

// Solution

    HLL_method (N, gimel, cons_vars, C);

    primitive_variables * final_features = new primitive_variables[N];
    cons2prim (N, cons_vars, final_features, gimel);

    save_results (test_nmbr, x, final_features);

// Deallocation of memory

    delete [] init_features;
    delete [] x;
    free_array(cons_vars.u);
    free_array(cons_vars.F);
    delete [] final_features;
    
    return 0;
}