#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "headers/variables.hpp"

using namespace std;

// The function of reading data from a file

bool initialization (string& test, vector<double>& left, vector<double>& right)
{
    ifstream input(test);
    if (!input.is_open())
    {
        cerr << "Failed to open " << test << '\n';
        return false;
    }

    if (!(input >> left[0] >> left[1] >> left[2] >> right[0] >> right[1] >> right[2]))
    {
        cerr << "Error reading data from file " << test << '\n';
        return false;
    }
    input.close();

    return true;
}

// The function of writing data to a file

bool save_results (string test, vector<double> x,  vector<double> dens, vector<double> vel, \
                   vector<double> pres, int N)
{
    string res_path = "./solution/output" + test + ".txt";

    ofstream output(res_path);
    if (output.is_open())
    {
        for (size_t i = 0; i < N; ++i)
            output << x[i] << " " << dens[i] << " " << vel[i] << " " << pres[i] << '\n';
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
