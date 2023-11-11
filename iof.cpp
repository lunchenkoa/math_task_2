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
