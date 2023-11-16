#include <iostream>
#include <fstream>
#include <string>

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

bool save_results (string test, double* x, primitive_variables* states)
{
    string result = "";
    result = "output" + test + ".txt";

    string res_path = "../solution/" + result;

    ofstream output(res_path);
    if (output.is_open())
    {
        int array_length = sizeof(x) / sizeof(x[0]);
        for (size_t i = 0; i != array_length; ++i)
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
