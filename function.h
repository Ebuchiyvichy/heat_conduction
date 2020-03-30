//
// Created by irina on 21.03.2020.
//

#ifndef HEAT_FUNCTION_H
#define HEAT_FUNCTION_H

#endif //HEAT_FUNCTION_H
#include "vector_new.h"

double	K(double x, Date my_date)
{
		if (0 < x && x < my_date.x1)
			return my_date.k1;
		if (my_date.x1 < x && x < my_date.x2)
			return my_date.k1 * (x - my_date.x2) / (my_date.x1 - my_date.x2) + my_date.k2 * (x - my_date.x1) / (my_date.x2 - my_date.x1);
		if (my_date.x2 < x && x < my_date.L)
			return my_date.k2;
}

double	K_quasi(double x, Date my_date)
{
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
    if ((x > 0) && (x < my_data.L))
		return my_data.u0;
	return 0;
}

double u0_t(double x, Date my_data)
{
    return my_data.u0 - x * (my_data.L - x);
}

//левая прогонка
std::vector<double> progon(std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> f, int n, double y0)
{
    double m = n;
	std::vector<double> y(n+1);
	std::vector<double>	ksi(n);
	std::vector<double>	etta(n);

	ksi[m-2] = a[m-1]/b[m-1];
	etta[m-2] = - f[m-1]/b[m-1];

    for (int i = m - 2; i >= 1; i--)
	{
		ksi[i-1] = a[i] / (b[i] - c[i] * ksi[i]);
		etta[i-1] = (f[i] + c[i] * etta[i]) / (b[i] - c[i] * ksi[i]);
	}
/*	std::cout << "ksi:\t";
	for (int i = 2; i < n; i++)
		std::cout << ksi[i] << '\t';
	std::cout << std::endl;
	std::cout << "etta:\t";
	for (int i = 2; i < n; i++)
		std::cout << etta[i] << '\t';
	std::cout << std::endl;
	std::cout << "f14 " << f[n - 1] << std::endl;
	std::cout << "ksi14 " << ksi[n - 1] << std::endl;
	std::cout << "b14 " << b[n - 1] << std::endl;
	std::cout << "c14 " << c[n - 1] << std::endl;
	std::cout << "etta14 " << etta[n - 1] << std::endl;
 */
	y[0] = y0;
//	std::cout << "y" << n - 1 << " = " << y[n - 1] << std::endl;
	for (int i = 1; i != m-1; i++) {
//       std::cout << i << std::endl;
		y[i] = ksi[i - 1] * y[i - 1] + etta[i - 1];
//		std::cout << "y" << i << " = " << y[i] << std::endl;
	}
		return y;
}


