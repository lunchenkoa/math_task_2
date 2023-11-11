#include <iostream>
#include <fstream>
#include <string>

using namespace std;

// The function of reading data from a file

// Introduce the notation:
//    dens_L, vel_L, pres_L = gas density, gas velocity, gas pressure on the LEFT,
//                          there are left border conditions.
//    dens_R, vel_R, pres_R = gas density, gas velocity, gas pressure on the RIGHT,
//                          there are right border conditions.
//    time is a time point under study.

bool initialization (string test, double& dens_L, double& vel_L, double& pres_L, \
                                  double& dens_R, double& vel_R, double& pres_R)
{
    ifstream input(test);
    if (!input.is_open())
    {
        cerr << "Failed to open " << test << '\n';
        return false;
    }

    if (!(input >> dens_L >> vel_L >> pres_L >> dens_R >> vel_R >> pres_R))
    {
        cerr << "Error reading data from file " << test << '\n';
        return false;
    }
    input.close();

    return true;
}

bool save_results (string test, double* x, double* pressure, double* density, double* velocity)
{
    string result = "";
    if (test == "input1.txt")
    {
        result = "output1.txt";
    }
    else if (test == "input2.txt")
    {
        result = "output2.txt";
    }
    else if (test == "input3.txt")
    {
        result = "output3.txt";
    }

    string res_path = "./solution/" + result;

    ofstream output(res_path);
    if (output.is_open())
    {
        int array_length = sizeof(x) / sizeof(x[0]);
        for (size_t i = 0; i != array_length; ++i)
            output << x[i] << " " << density[i] << " " << velocity[i] << " " << pressure[i] << '\n';
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
