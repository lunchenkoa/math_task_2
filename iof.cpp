#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <Eigen/Dense>
 
using Eigen::MatrixXd;
using Eigen::VectorXd;

#include "headers/variables.hpp"

using namespace std;

// The function of reading data from a file

bool initialization (string& test, primitive_variables& left, primitive_variables& right)
{
    ifstream input(test);
    if (!input.is_open())
    {
        cerr << "Failed to open " << test << '\n';
        return false;
    }

    if (!(input >> left.dens >> left.vel >> left.pres >> right.dens >> right.vel >> right.pres))
    {
        cerr << "Error reading data from file " << test << '\n';
        return false;
    }
    input.close();

    return true;
}

// The function of writing data to a file

bool save_results (string test, VectorXd x, vector<primitive_variables> states, int N)
{
    string res_path = "./solution/output" + test + ".txt";

    ofstream output(res_path);
    if (output.is_open())
    {
        for (size_t i = 0; i < N; ++i)
            output << x[i] << " " << states[i].dens << " " << states[i].vel << " " << states[i].pres << '\n';
    }
    else
    {
        cerr << "Failed to open " << res_path << '\n';
        return false;
    }
    output << '\n';
    output.close();

    return true;
}
