//
// Created by irina on 28.03.2020.
//

#ifndef HEAT_INTEGRO_INTERPOLATION_H
#define HEAT_INTEGRO_INTERPOLATION_H

#endif //HEAT_INTEGRO_INTERPOLATION_H

#include "function.h"

void integro_interpolation(int n, int t,  double h, double tau, int TEST_P, Date my_date)
{
    std::ofstream		fout;
    std::vector<double>	a;
    std::vector<double> B(n);
    std::vector<double> A(n+1);
    std::vector<double> C(n);
    std::vector<double> F(n);
    std::vector<double> y1(n+1);	// значения на текущем временном слое
    std::vector<double> y2(n+1);	// следующий временной слой

    // инициализация начальными данными
	for (int i = 0; i != n; i++) {
		y1[i] = u0_t(i*h, my_date);
		std::cout << y1[i] << '\t';
	}
	std::cout << std::endl;
    double sigma = 0.5;
	for (int i = 0; i != n + 1; i++) {
		a.push_back(K(i * h - 0.5 * h, my_date));
		std::cout << a[i] << '\t';
	}
	std::cout << std::endl;

    fout.open("Integtgro_interpolation_mult.txt");

    //коэффициенты прогонки
    A[0] = 0;
    for (int i = 1; i <= n; i++)
    {
        A[i] = sigma / h * a[i];
        B[i-1] = A[i];
        C[i-1] = A[i-1] + B[i-1] + my_date.c* my_date.rho * h/tau;
    }
    std::cout << "Coefficients was found\n";

    //вычисление по временным слоям
    for (int j = 0; j != t; j++)
    {

        // инициализация функции правой части
        for (int i = 1; i != n-1; i++) {
            F[i] = my_date.c * my_date.rho * y1[i] * h / tau + (1 - sigma) * a[i] * (y1[i + 1] - 2 * y1[i] + y1[i - 1]) / h;
        //   std::cout << "F = " << F[i] << std::endl;
        }
        std::cout << "F is ready" << std::endl;
        //передача значений с 1 по n-1, так как они уже определены
        y2 = progon(A, C, B, F, n, my_date.left_boarder(0, my_date));
        std::cout << "ok\n" << std::endl;
        //y2[0] = my_date.left_boarder(0, my_date);
        if ((TEST_P == 1) || (TEST_P == 3)) // TEST_P == 1 - смешанная задача
                                            // TEST_P == 3 - два потока на концах
        {
            std::cout << "progon is ready" << std::endl;
            double kappa = (sigma * a[n - 1] / h) / (my_date.c * my_date.rho * h / (2 * tau) + sigma * a[n - 2] / h);
            double mu = (my_date.c * my_date.rho * y1[n-1] * h / (2 * tau) + sigma * my_date.right_boarder(tau * j, my_date) +
                  (1 - sigma) * (my_date.right_boarder(tau * (j - 1), my_date) - a[n-1] * (y1[n-1] - y1[n-2]) / h))
                 / (my_date.c * my_date.rho * h / (2 * tau) + sigma * a[n-1] / h);
            std::cout << "mu1 = " << mu << std::endl;
            y2[n] = kappa * y2[n-1] + mu;
            if (TEST_P == 3)//из условий потока(?)
            {
                kappa = (sigma * a[0] / h) / (my_date.c * my_date.rho * h / (2 * tau) + sigma * a[0] / h);
                mu = (my_date.c * my_date.rho * y1[0] * h / (2 * tau) + sigma * my_date.left_boarder(tau * j, my_date) +
                              (1 - sigma) * (my_date.left_boarder(tau * (j - 1), my_date) - a[0]*(y1[1] - y1[0]) / h))
                             / (my_date.c * my_date.rho * h / (2 * tau) + sigma * a[0] / h);
                y2[0] = kappa * y2[1] + mu;
                std::cout << "mu2 = " << mu << std::endl;
            }
        }
        else if (TEST_P == 2)
            y2[n] = my_date.right_boarder(n*h, my_date);

            for (int i = 0; i != y2.size(); i++)
              fout << j << '\t' << i << '\t' << y2[i] << '\n';
        y1 = y2;
        std::cout << j << std::endl;
    }
    fout.close();
}