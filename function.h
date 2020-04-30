//
// Created by irina on 21.03.2020.
//

#ifndef HEAT_FUNCTION_H
#define HEAT_FUNCTION_H

#endif //HEAT_FUNCTION_H
#include "vector_new.h"

double	K(double x, Date my_date)
{
		if ((x < my_date.x1) || fabs(x-my_date.x1) <= EPS)
			return my_date.k1;
		if (my_date.x1 < x && x < my_date.x2)
			return my_date.k1 * (x - my_date.x2) / (my_date.x1 - my_date.x2) + my_date.k2 * (x - my_date.x1) / (my_date.x2 - my_date.x1);
		if ((my_date.x2 < x) || fabs(x-my_date.x2) <= EPS)
			return my_date.k2;
  //  return  1;
}

double	K_quasi(double x, Date my_date)
{
	//return (my_date.alpha + my_date.beta * pow(x, my_date.gamma));
	return (my_date.kappa * pow(x,my_date.sigma));
}

double	P1(double t, Date my_data)
{
	if (t > EPS && t < my_data.t0)
		return my_data.Q;
	else
		return 0;
}

double	P2(double t, Date my_data)
{
	if (EPS < t && t < my_data.t0)
		return 2 * my_data.Q * t;
	else
		return 0;
}

double	P3(double t, Date my_data)
{
	if (EPS < t && t < my_data.t0)
		return 2 * my_data.Q * (my_data.t0 - t);
	else return 0;
}

double	P4(double t, Date my_data)
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
	/*if (fabs(x-0) <= EPS)
		return my_data.u0;
	else if (fabs(my_data.L-x)<= EPS)
		return my_data.u0;*/

	return my_data.u0*pow(x, 1/my_data.sigma);	//метода
//	return my_data.u0;	//Ирин вариант
// 	return 0;
}

double	u0_t(double x, Date my_data)
{
    //return my_data.u0 + x * (my_data.L - x);
    //return sin(x);
//	return my_data.u0*pow(x, 1/my_data.sigma);

	return 0;	//метода
//	return my_data.u0;	//Ирин вариант
}

// // правая прогонка
// std::vector<double> progon(std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> f, int n, double kappa, double  mu, DATE my_date)
// {
// 	std::vector<double> y(n+1);
// 	std::vector<double>	alpha(n);
// 	std::vector<double>	beta(n);

// 	alpha[0] = c[0]/b[0];
// 	beta[0] = f[0]/b[0];
// 	f[1]+=a[1] * my_date.left_boarder(0, my_date);
// 	a[1] = 0;
// 	y[0] = my_date.left_boarder(0, my_date);
// 	for (int i = 1; i <= n-1; i++)
// 	{
// 		alpha[i] = c[i-1]/(b[i-1]-alpha[i-1]*a[i-1]);
// 		beta[i] = (f[i-1]+a[i-1]*beta[i-1])/(b[i-1]-a[i-1]*alpha[i-1]);

// 	}
// 	for (int i = 1; i <= n; i++)
// 	{
// 		y[i] = alpha[i+1] x[]
// 	}



// }


//левая прогонка
std::vector<double> progon(std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> f, int n, double y0, double yn, double kappa, double  mu)
{
	std::vector<double> y(n+1);
	std::vector<double>	ksi(n);
	std::vector<double>	etta(n);

	f[1]+=a[1]*y0;
	a[1] = 0;
	f[n-1] += a[n-1]*yn;
	ksi[n - 1] = a[n - 1] / (b[n - 1] - c[n - 1] * kappa);
	etta[n - 1] = (f[n - 1] + c[n - 1] * mu) / (b[n - 1] - c[n - 1] * kappa);

    for (int i = n - 2; i >= 1; i--)
	{
		ksi[i] = a[i] / (b[i] - c[i] * ksi[i+1]);
		etta[i] = (f[i] + c[i] * etta[i+1]) / (b[i] - c[i] * ksi[i+1]);
	}
	y[0] = y0;

	for (int i = 1; i != n-1; i++)
		y[i] = ksi[i] * y[i-1] + etta[i];
	return y;
}


