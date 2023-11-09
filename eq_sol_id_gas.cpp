
#include <iostream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

const double adiabat = 5.0 / 3.0;  // Показатель адиабаты
const double tol = 1e-6;

// Determining initial conditions for tests
void initial_conditions(string test_numbr, double& ro_l,\
                        double& v_l,double& p_l,double& ro_r,\
                        double& v_r,double& p_r,double& t)
{
    // configuration definition
    if (test_numbr.compare("1") == 0)
    {
        ro_l = 1.0; v_l = 0.0; p_l = 3.0; 
        ro_r = 1.0; v_r = 0.0; p_r = 1.0;
        t = 0.1;

    } else if (test_numbr.compare("2") == 0)
    {
        ro_l = 1.0; v_l = 1.0; p_l = 3.0; 
        ro_r = 1.0; v_r = -1.0; p_r = 1.0;
        t = 0.1;

    } else if (test_numbr.compare("3") == 0)
    {
        ro_l = 1.0; v_l = -0.1; p_l = 1.0; 
        ro_r = 1.0; v_r = 0.2; p_r = 1.0;
        t = 0.1;
    
    } else if (test_numbr.compare("4") == 0)
    {
        ro_l = 1.0; v_l = -0.1; p_l = 1.0; 
        ro_r = 1.0; v_r = 0.2; p_r = 3.0;
        t = 0.18;
    } else 
    {
        cout << "Invalid test number entered"<< endl;
    }

}

// Функция для вычисления скорости звука
double sound_speed(double density, double pressure) {
    return std::sqrt(adiabat * pressure / density);
}

// Функция для вычисления максимальной скорости
double max_speed(double density, double velocity, double pressure) {
    return std::max(std::fabs(velocity), sound_speed(density, pressure));
}

// Функция для вычисления потоков через грани ячейки
void compute_fluxes(double density_L, double velocity_L, double pressure_L,
                    double density_R, double velocity_R, double pressure_R,
                    double &flux_density, double &flux_momentum, double &flux_energy) {
    double sound_speed_L = sound_speed(density_L, pressure_L);
    double sound_speed_R = sound_speed(density_R, pressure_R);
    
    double p_star = 0.5 * (pressure_L + pressure_R - 
                    density_L * velocity_L * sound_speed_L - 
                    density_R * velocity_R * sound_speed_R) / 
                    (0.5 * (sound_speed_L + sound_speed_R));
    
    if (p_star <= 0.0) {
        p_star = tol;
    }
    
    double v_star = 0.5 * (velocity_L + velocity_R + 
                    (pressure_L - pressure_R) / 
                    (density_L * sound_speed_L + density_R * sound_speed_R));
    
    double rho_star = 0.5 * (density_L + density_R);
    
    if (v_star >= 0.0) {
        flux_density = density_L * velocity_L;
        flux_momentum = density_L * velocity_L * velocity_L + pressure_L;
        flux_energy = velocity_L * (density_L * (pressure_L + density_L * velocity_L * velocity_L) / (density_L * pressure_L));
    } else {
        flux_density = density_R * velocity_R;
        flux_momentum = density_R * velocity_R * velocity_R + pressure_R;
        flux_energy = velocity_R * (density_R * (pressure_R + density_R * velocity_R * velocity_R) / (density_R * pressure_R));
    }
    
    if (v_star > 0.0) {
        flux_density += p_star - pressure_L;
        flux_momentum += (p_star + density_L * v_star * v_star) * v_star - density_L * velocity_L * (v_star - velocity_L);
        flux_energy += v_star * (p_star + density_L * v_star * v_star);
    } else {
        flux_density += p_star - pressure_R;
        flux_momentum += (p_star + density_R * v_star * v_star) * v_star - density_R * velocity_R * (v_star - velocity_R);
        flux_energy += v_star * (p_star + density_R * v_star * v_star);
    }
}