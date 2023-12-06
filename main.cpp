#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <Eigen/Dense>
 
#include "headers/variables.hpp"
#include "headers/iof.hpp"
#include "headers/functions.hpp"
#include "headers/de_allocate.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;
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

    vector<double> l_init, r_init;
    l_init.resize(3);
    r_init.resize(3);
    initialization (test_file, l_init, r_init);

// Allocation of memory to dynamic variables

    // double * x = create_vector(N);
    VectorXd x(N);
    // cout << "Init Values: " << endl;1
    for (size_t i = 0; i < N; ++i)
    {
        x(i) = x_L + (i + 0.5) * dx;
    }

    // vector<primitive_variables> init_features;
    // init_features.resize(N);
    // primitive_variables init_features;
    vector<double> dens, vel, pres;
    dens.resize(N);
     vel.resize(N);
    pres.resize(N);
    // primitive_variables * init_features = new primitive_variables[N]; // array for init {rho, v, p}
    initialization_of_IC( x, N, dens, vel, pres, l_init, r_init);           // which is half filled with 
                                                                      // left characteristics and 
                                                                      // half with right ones

    // conservative_variables cons_vars;
    double ** u;
    double** F;
    u = create_array(N, 3);
    F = create_array(N, 3);
    // MatrixXd  u[N,3);
    // MatrixXd F(N,3);

// Solution

    // prim2cons (N, cons_vars, init_features, gimel); // fill the cons_vars
    prim2cons (N, u, F, dens, vel, pres, gimel); // fill the cons_vars


    // for (size_t i = 0; i < N; ++i)
    // {
    //     // cout << "i = " << i << " x = " << x[i] << " P = " << init_features[i].pres << " Rho = " << init_features[i].dens << " V = " << init_features[i].vel << endl;
    //     // cout << "i = " << i << " x = " << x[i] << " u = " << cons_vars.u[i][0] << ", " << cons_vars.u[i][1] << ", " << cons_vars.u[i][2] << ", " <<  " F = " <<  cons_vars.F[i][0] << ", " << cons_vars.F[i][1] << ", " << cons_vars.F[i][2] << endl;
    // }

    HLL_method (N, gimel, u, F, C);

    // primitive_variables * final_features = new primitive_variables[N];
    primitive_variables final_features;

    vector<double> f_dens, f_vel, f_pres;
    f_dens.resize(N);
    f_vel.resize(N);
    f_pres.resize(N);

    // final_features.resize(N);
    cons2prim (N, u, F, f_dens, f_vel, f_pres, gimel);

    save_results (test_nmbr, x, f_dens, f_vel, f_pres, N);

// Deallocation of memory

    // delete [] init_features;
    // delete [] x;
    free_array(u);
    free_array(F);
    // delete [] final_features;
    
    return 0;
}
