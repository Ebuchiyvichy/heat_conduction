//
// Created by irina on 28.03.2020.
//

#ifndef HEAT_INTEGRO_INTERPOLATION_H
#define HEAT_INTEGRO_INTERPOLATION_H

#endif //HEAT_INTEGRO_INTERPOLATION_H

#include "function.h"

std::vector<double> integro_interpolation(int n, double T,  double h, double tau, int TEST_P, double sigma, Date my_date)
{
    std::ofstream		fout;
    std::vector<double>	a(n+1);
    std::vector<double> B(n);
    std::vector<double> A(n+2);
    std::vector<double> C(n+1);
    std::vector<double> F(n+1);
    std::vector<double> y1(n+1);	// значения на текущем временном слое
    std::vector<double> y2(n+1);	// следующий временной слой
	std::string			str = "Integro_interpolation_mult_";

	str += std::to_string(sigma) + std::to_string(T) + ".txt";

    // инициализация начальными данными
	for (int i = 0; i <= n; i++)
		y1[i] = u0_t(i*h, my_date);
	for (int i = 0; i != n + 1; i++)
		a[i] = K(i * h - 0.5 * h, my_date);

    fout.open(str);

    //коэффициенты прогонки
    A[0] = 0;
    for (int i = 1; i <= n; i++)
    {
        A[i] = sigma / h * a[i];
        B[i-1] = A[i];
        C[i-1] = A[i-1] + B[i-1] + my_date.c * my_date.rho * h / tau;
    }
    std::cout << "Coefficients was found\n";
    std::vector<double> tmp(n);
    //вычисление по временным слоям
    for (double j = tau; j <= T; j += tau)
    {
        // инициализация функции правой части
        for (int i = 1; i != n; i++)
            F[i] = my_date.c * my_date.rho * y1[i] * h / tau + (1 - sigma) * a[i] * (y1[i + 1] - 2 * y1[i] + y1[i - 1]) / h;
        if ((TEST_P == 1) || (TEST_P == 3)) // TEST_P == 1 - смешанная задача
                                            // TEST_P == 3 - два потока на концах
        {
            double kappa = (sigma * a[n] / h) / (my_date.c * my_date.rho * h / (2 * tau) + sigma * a[n] / h);
			double mu = (my_date.c * my_date.rho * y1[n - 1] * h / (2 * tau) + sigma * my_date.right_boarder(j, my_date) +
				(1 - sigma) * (my_date.right_boarder((j - tau), my_date) - a[n] * (y1[n] - y1[n - 1]) / h))
				/ (my_date.c * my_date.rho * h / (2 * tau) + sigma * a[n] / h);
			y2 = progon(A, C, B, F, n, my_date.left_boarder(j, my_date), my_date.right_boarder(j,my_date), kappa, mu);
			y2[n] = kappa * y2[n - 1] + mu;
            if (TEST_P == 3)//из условий потока(?)
            {
                kappa = (sigma * a[0] / h) / (my_date.c * my_date.rho * h / (2 * tau) + sigma * a[0] / h);
                mu = (my_date.c * my_date.rho * y1[0] * h / (2 * tau) + sigma * my_date.left_boarder(j, my_date) +
                              (1 - sigma) * (my_date.left_boarder(j - tau, my_date) - a[0]*(y1[1] - y1[0]) / h))
                             / (my_date.c * my_date.rho * h / (2 * tau) + sigma * a[0] / h);
                y2[0] = kappa * y2[1] + mu;
            }
        }
		else if (TEST_P == 2)
		{
			y2 = progon(A, C, B, F, n, my_date.left_boarder(j, my_date),my_date.right_boarder(j,my_date), 0, 0);
			y2[n] = my_date.left_boarder(j, my_date);
		}
        // for (int i = 0; i != y1.size(); i++)
        //     fout << '\t' << i * h << '\t' << y1[i] << '\n';
        y1 = y2;
    }
    for (int i = 0; i != y1.size(); i++)
        fout << '\t' << i * h << '\t' << y1[i] << '\n';
    fout.close();
    return y1;
}