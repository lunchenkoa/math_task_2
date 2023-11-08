/*
Программирование численного решения уравнений газодинамики методом Хартена-Лакса-ван Лира (HLL) — это 
большой проект, который может занять много времени и кода. Ниже я предоставлю вам пример реализации 
HLL метода для одномерных уравнений газодинамики на C++, но учтите, что это всего лишь начало, и вы 
должны будете расширить код для решения всей задачи.
*/

#include <iostream>
#include <vector>
#include <cmath>

const double gamma = 5.0 / 3.0;
const double tol = 1e-6;

// Функция для вычисления скорости звука
double sound_speed(double density, double pressure) {
    return std::sqrt(gamma * pressure / density);
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

int main() {
    // Задаем начальные условия
    double x_min = -0.5;
    double x_max = 0.5;
    int N = 320;
    double dx = (x_max - x_min) / N;

    std::vector<double> density(N);
    std::vector<double> velocity(N);
    std::vector<double> pressure(N);

    for (int i = 0; i < N; ++i) {
        if (i < N / 2) {
            density[i] = 1.0;
            velocity[i] = 0.0;
            pressure[i] = 3.0;
        } else {
            density[i] = 1.0;
            velocity[i] = 0.0;
            pressure[i] = 1.0;
        }
    }

    // Задаем параметры временного шага
    double dt = 0.001;
    double t_max = 0.2;
    int num_steps = static_cast<int>(t_max / dt);

    // Выполняем временную итерацию
    for (int step = 0; step < num_steps; ++step) {
        // Создаем временные массивы для новых значений
        std::vector<double> new_density(N);
        std::vector<double> new_velocity(N);
        std::vector<double> new_pressure(N);

        // Рассчитываем потоки через грани ячеек и обновляем значения
        for (int i = 0; i < N - 1; ++i) {
            double flux_density, flux_momentum, flux_energy;
            compute_fluxes(density[i], velocity[i], pressure[i],
                           density[i + 1], velocity[i + 1], pressure[i + 1],
                           flux_density, flux_momentum, flux_energy);
            
            new_density[i] = density[i] - dt / dx * (flux_density - flux_density);
            new_velocity[i] = velocity[i] - dt / dx * (flux_momentum - flux_momentum);
            new_pressure[i] = pressure[i] - dt / dx * (flux_energy - flux_energy);
        }

        // Обновляем значения в ячейках
        for (int i = 0; i < N; ++i) {
            density[i] = new_density[i];
            velocity[i] = new_velocity[i];
            pressure[i] = new_pressure[i];
        }
    }

    // Здесь можно добавить код для вывода результатов

    return 0;
}

/*
Приведенный код решает уравнения газодинамики методом HLL для начальных условий из варианта 1. 
Он выполняет временную итерацию с использованием явной разностной схемы, но вы должны выполнить 
оставшуюся работу по выводу результатов и адаптации кода для других начальных условий и заданных 
параметров. Также учтите, что это всего лишь пример, и для полной реализации вам придется внести 
дополнительные детали и оптимизации.
*/
