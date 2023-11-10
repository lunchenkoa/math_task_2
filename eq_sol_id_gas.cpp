#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include "headers/de_allocate.hpp"

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
        t = 0.1;
    } else 
    {
        cout << "Invalid test number entered"<< endl;
    }

    /*
    Каким-то образом надо сделать теперь эти начальные условия внешними переменными, чтобы каждый раз их не передавать туда сюда.
    Мб они должны быть типа global.... Хз, есть ли такое - надо гуглить.
    Наверное, стоит их объединить в вектор u?...
    */

}

// Функция для вычисления скорости звука
double sound_speed(double density, double pressure) 
{
    return sqrt(adiabat * pressure / density);
}


double* matrix_eigenvalue (double* u, double density, double velocity, double pressure)
{
    double *lambda = new double[3];
    lambda[0] = velocity - sound_speed(density, pressure);
    lambda[1] = velocity + sound_speed(density, pressure);
    lambda[2] = velocity;
    
    return lambda;
}

// Функция для вычисления максимальной скорости
double* speed_estimates(double* u_l, double* u_r) 
{   
    /* 
    Там в формуле на слайде 56 под максимумом стоит какое-то а...
    Я так и не пон, что это значит. Типо, скорости D_l и D_r - трехмерные???
    */
    double* D = new double[2];   
    D[0] = min(fabs(matrix_eigenvalue(u_l)), fabs(matrix_eigenvalue(u_r))); 
    D[1] = max(fabs(matrix_eigenvalue(u_l)), fabs(matrix_eigenvalue(u_r))); 
    
    return D;
}

// Функция для вычисления потоков через грани ячейки
void compute_fluxes(double density_L, double velocity_L, double pressure_L,
                    double density_R, double velocity_R, double pressure_R,
                    double &flux_density, double &flux_momentum, double &flux_energy) 
{
    double sound_speed_L = sound_speed(density_L, pressure_L);
    double sound_speed_R = sound_speed(density_R, pressure_R);
    
    double p_star = 0.5 * (pressure_L + pressure_R - 
                    density_L * velocity_L * sound_speed_L - 
                    density_R * velocity_R * sound_speed_R) / 
                    (0.5 * (sound_speed_L + sound_speed_R));
    
    if (p_star <= 0.0) 
    {
        p_star = tol;
    }
    
    double v_star = 0.5 * (velocity_L + velocity_R + 
                    (pressure_L - pressure_R) / 
                    (density_L * sound_speed_L + density_R * sound_speed_R));
    
    double rho_star = 0.5 * (density_L + density_R);
    
    if (v_star >= 0.0) 
    {
        flux_density = density_L * velocity_L;
        flux_momentum = density_L * velocity_L * velocity_L + pressure_L;
        flux_energy = velocity_L * (density_L * (pressure_L + density_L * velocity_L * velocity_L) / (density_L * pressure_L));
    } else 
    {
        flux_density = density_R * velocity_R;
        flux_momentum = density_R * velocity_R * velocity_R + pressure_R;
        flux_energy = velocity_R * (density_R * (pressure_R + density_R * velocity_R * velocity_R) / (density_R * pressure_R));
    }
    
    if (v_star > 0.0) 
    {
        flux_density += p_star - pressure_L;
        flux_momentum += (p_star + density_L * v_star * v_star) * v_star - density_L * velocity_L * (v_star - velocity_L);
        flux_energy += v_star * (p_star + density_L * v_star * v_star);
    } else 
    {
        flux_density += p_star - pressure_R;
        flux_momentum += (p_star + density_R * v_star * v_star) * v_star - density_R * velocity_R * (v_star - velocity_R);
        flux_energy += v_star * (p_star + density_R * v_star * v_star);
    }
}


void godunov(double *F, int N)
{

    double dt, t = 0.0;
    double ** temp_u = create_array(N, 3);
    double ** F_l = create_array(N, 3);
    double ** F_r = create_array(N, 3);

    //double * temp_u = new double[N];
    double * ro = new double[N];
    double * v = new double[N];
    double * p = new double[N];
    double D_l, D_r;
    int counter = 0;

    while (t <= time)                 // time  - это по идее всеря 0.1c
    {
        // Надо согласовать все по параметрам 
        D_l, D_r = speed_estimates(u_l, u_r);  // и да, так нельзя, но я на сегодня все...

        ! Íàõîäèì øàã ïî âðåìåíè
        dt = Curant * dx  / v_max
        ! Íàõîäèì òåêóùåå çíà÷åíèå âðåìåíè
        t = t + dt
        ! Âû÷èñëÿåì F- è F+ äëÿ òåêóùåãî ìîìåíòà âðåìåíè
        do j=1, number_of_nods
                F_r(j, :) = (F(j,:) + F(j+1,:))/2d0 - v_max/2d0*(u(j+1, :) - u(j, :))
                F_l(j, :) = (F(j-1,:) + F(j,:))/2d0 - v_max/2d0*(u(j, :) - u(j-1, :))
        end do
        ! Âû÷èñëÿåì âåêòîðíóþ âåëè÷èíó u â äàííûé ìîìåíò âðåìåíè
        do j=1, number_of_nods
            temp_u(j, :) = u(j, :) - dt/dx * (F_r(j, :) - F_l(j, :))
        end do

        u(1:number_of_nods, :) = temp_u(1:number_of_nods,:)
        u(0,:) = u(1,:)
        u(number_of_nods+1,:) = u(number_of_nods,:)
        ! Ïåðåâîäèì u â ðåàëüíûå âåëè÷èíû
        call u_to_real(u,ro,v,p)
        ! Ïî ðåàëüíûì âåëè÷èíàì îáíîâëÿåì âåêòîðíóþ âåëè÷èíó F
        call real_to_F(F,ro,v,p)

        counter = counter + 1
    end do

    deallocate(temp_u, F_l, F_r)

}

subroutine u_to_real(u, ro, v, p)
    ! Ïîäïðîãðàììà äëÿ ïåðåâîäà u-->(ro,v,p)
    real(8), dimension(0:)       :: ro, v, p
    real(8), dimension(0:,:)     :: u
    integer sizee, i

    sizee = size(u(:,1))

    ! Âû÷èñëÿåì ïëîòíîñòü
    ro = u(:, 1)
    ! Âû÷èñëÿåì ñêîðîñòü
    do i=0, sizee-1
        v(i) = u(i, 2)/ro(i)
    end do
    ! Âû÷èñëÿåì äàâëåíèå
    do i=0, sizee-1
        p(i) = (u(i, 3) - ro(i) * (v(i))**2/2) * (gamm-1)
    end do

end subroutine

subroutine real_to_u(u, ro, v, p)
    ! Ïîäïðîãðàììà äëÿ ïåðåâîäà (ro,v,p)-->u
    real(8), dimension(0:)       :: ro, v, p
    real(8), dimension(0:,:)     :: u

    ! Âû÷èñëÿåì u_1
    u(:, 1) = ro
    ! Âû÷èñëÿåì u_2
    u(:, 2) = ro * v
    ! Âû÷èñëÿåì u_3
    u(:, 3) = p / (gamm-1) + ro / 2 * v**2

end subroutine

subroutine real_to_F(F, ro, v, p)
    ! Ïîäïðîãðàììà äëÿ ïåðåâîäà (ro,v,p)--> F
    real(8), dimension(0:)       :: ro, v, p
    real(8), dimension(0:,:)     :: F


    ! Âû÷èñëÿåì F_1
    F(:, 1) = ro * v
    ! Âû÷èñëÿåì F_2
    F(:, 2) = ro * v**2 + p
    ! Âû÷èñëÿåì F_3
    F(:, 3) = v * (p/(gamm-1) + p + ro/2 * v**2)

end subroutine

function vmax(u) result(v_max)
    ! Ôóíêöèÿ äëÿ íàõîæäåíèÿ ìàêñèìàëüíîé ñêîðîñòè ïî âñåì ÿ÷åéêàì
    real(8), dimension(0:,:) :: u
    real(8), allocatable, dimension(:) :: ro1, v1, p1
    real(8)                     :: v_max, c
    integer sizee, i

    v_max = -1
    sizee = size(u(:,1))
    allocate(ro1(0:sizee-1), v1(0:sizee-1), p1(0:sizee-1))

    call u_to_real(u, ro1, v1, p1)
    do i=1, sizee-2
        c = vol_speed(ro1(i), p1(i))
        v_max = max(v_max, abs(v1(i)) + c)

    end do
    deallocate( ro1, v1, p1)
end function

function vol_speed(ro, p) result(C_s)
    ! Ôóíêöèÿ âû÷èñëåíèÿ ñêîðîñòè çâóêà
    real(8) :: ro, p, C_s
    C_s = sqrt(gamm * p/ro)

end function

subroutine set_initial_values(ro, v, p)
    ! Ïîäïðîãðàììà óñòàíîâêè íà÷àëüíûõ çíà÷åíèé
    real(8), dimension(0:) :: ro, v, p
    integer sizee, i

    sizee = size(ro)

    do i=0, N0
        ro(i) = u0(1,1)
        v(i)  = u0(2,1)
        p(i)  = u0(3,1)
    end do

    do i=N0+1, sizee-1
        ro(i) = u0(1,2)
        v(i)  = u0(2,2)
        p(i)  = u0(3,2)
    end do


end subroutine

subroutine make_dat(x, ro, v, p, suffix)
    ! Ïîäïðîãðàììà çàïèñè â ôàéë
    real(8), dimension(0:) :: x, ro, v, p
    character(16) suffix

    xsize=size(x)
    open(unit=10, file='case_'//trim(case_number)//trim(suffix)//'.dat')
    write(10,'(4x,a3,13x,a4,12x,a3,13x,a3)') '# x', '# ro', '# v', '# p'
    do i=0, xsize-1
        write(10, '(4e16.8)') x(i), ro(i), v(i), p(i)
    end do
    close(10)

end subroutine



end module