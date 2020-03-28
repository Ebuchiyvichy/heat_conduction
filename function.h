//
// Created by irina on 21.03.2020.
//

#ifndef HEAT_FUNCTION_H
#define HEAT_FUNCTION_H

#endif //HEAT_FUNCTION_H
#include "vector_new.h"

double	K(double x, Date my_date)
{
	if (my_date.test == 'a') {
		if (0 < x && x < my_date.x1)
			return my_date.k1;
		if (my_date.x1 < x && x < my_date.x2)
			return my_date.k1 * (x - my_date.x2) / (my_date.x1 - my_date.x2) + my_date.k2 * (x - my_date.x1) / (my_date.x2 - my_date.x1);
		if (my_date.x2 < x < my_date.L)
			return my_date.k2;
	}
	else
    return (my_date.alpha + my_date.beta * pow(x, my_date.gamma));
}

double P1(double t, Date my_data)
{
	if (t > EPS && t < my_data.t0)
		return my_data.Q;
	else
		return 0;
}

double P2(double t, Date my_data)
{
	if (EPS < t && t < my_data.t0)
		return 2 * my_data.Q * t;
	else
		return 0;
}

double P3(double t, Date my_data)
{
	if (EPS < t && t < my_data.t0)
		return 2 * my_data.Q * (my_data.t0 - t);
	else return 0;
}

double P4(double t, Date my_data)
{
		if (EPS < t && t < 0.5 * my_data.t0)
			return 2 * my_data.Q * t;
		else if (0.5 * my_data.t0 < t && t < my_data.t0)
			return 2 * my_data.Q * (my_data.t0 - t);
		else
			return 0;
}

double	u0(double x, Date my_data)
{
	return my_data.u0;
}

std::vector<double> progon(std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> f, int n)
{
	double m;
	std::vector<double> y(n+1);
	std::vector<double>	alpha(n);
	std::vector<double>	beta(n);

	alpha[2] = b[1] / c[1];
	beta[2] = f[1] / c[1];

    for (int i = 3; i < n; i++)
	{
		alpha[i] = b[i] / (c[i] - a[i] * alpha[i - 1]);
		beta[i] = (f[i] + a[i] * beta[i - 1]) / (b[i] - a[i] * alpha[i - 1]);
	}
	/*std::cout << "alpha:\t";
	for (int i = 2; i < n; i++)
		std::cout << alpha[i] << '\t';
	std::cout << std::endl;
	std::cout << "beta:\t";
	for (int i = 2; i < n; i++)
		std::cout << beta[i] << '\t';
	std::cout << std::endl;
	std::cout << "f14 " << f[n - 1] << std::endl;
	std::cout << "alpha14 " << alpha[n - 1] << std::endl;
	std::cout << "b14 " << b[n - 1] << std::endl;
	std::cout << "c14 " << c[n - 1] << std::endl;
	std::cout << "beta14 " << beta[n - 1] << std::endl;*/
	y[n - 1] = (f[n - 1] + a[n - 1] * beta[n - 1]) / (c[n - 1] - a[n - 1] * alpha[n - 1]);
//	std::cout << "y" << n - 1 << " = " << y[n - 1] << std::endl;
	for (int i = n - 2; i >= 1; i--) {
		y[i] = alpha[i + 1] * y[i + 1] + beta[i + 1];
//		std::cout << "y" << i << " = " << y[i] << std::endl;
	}
		return y;
}

void integro_interpolation(int n, int t,  double h, double tau, int TEST_P, Date my_date)
{
    std::ofstream		fout;
	double				mu;
	std::vector<double>	a;
	std::vector<double>	y(n+1);
	std::vector<double> A(n+1);
    std::vector<double> B(n);
    std::vector<double> C(n);
    std::vector<double> F(n);
    std::vector<double> y1(n+1);	// значения на текущем временном слое
	std::vector<double> y2(n+1);	// следующий временной слой
	
	// инициализация начальными данными
	for (int i = 0; i != n+1; i++)
		y1[i] = my_date.left_boarder(i, my_date);

	double sigma = 0.5;
	for (int i = 0; i != n+1; i++)
		a.push_back(K(my_date.u0 + i * h - 0.5 * h, my_date));
	double kappa = (sigma * a[n - 1] / h) / (my_date.c * my_date.rho * h / (2 * tau) + sigma * a[n - 1] / h);
	fout.open("Integtgro_interpolation_mult.txt");

	//коэффициенты прогонки
	A[0] = 0;
	A[1] = sigma / h * a[1];
	for (int i = 2; i <= n; i++)
	{
		A[i] = sigma / h * a[i];
		B[i-1] = A[i];
		C[i-1] = A[i-1] + B[i-1] - my_date.c* my_date.rho * h/tau;
	}
	std::cout << "Coefficients was found\n";

	//вычисление по временным слоям
	for (int j = 0; j != t; j++)
	{
		// инициализация функции правой части
		for (int i = 1; i < n; i++) {
			F[i] = -my_date.c * my_date.rho * y1[i] * h / tau - (1 - sigma) * a[i] * (y1[i + 1] - 2 * y1[i] + y[i - 1]) / h;
			std::cout << "F = " << F[i] << std::endl;
		}

		mu = (my_date.c * my_date.rho * y1[n] * h / (2 * tau) + sigma * my_date.right_boarder(tau * j, my_date) + 
			(1 - sigma) * (my_date.right_boarder(tau * (j - 1), my_date) - a[n] * (y1[n] - y1[n-1]) / h))
			/ (my_date.c * my_date.rho * h / (2 * tau) + sigma * a[n] / h);
		std::cout << "mu = " << mu << std::endl;

        //передача значений с 1 по n-1, так как они уже определены
		y2 = progon(A, C, B, F, n);
		y2[0] = my_date.left_boarder(0, my_date);
		y2[n] = kappa * y2[n-1] + mu;
		//  print_in_file(y2, fout);//??????????
		for (int i = 0; i != y2.size(); i++)
			fout << j << '\t' << i << '\t' << y2[i] << '\n';
		y1 = y2;
    }
    fout.close();
}

