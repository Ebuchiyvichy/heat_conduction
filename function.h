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
	return  1;
}

double	K_quasi(double x, Date my_date)
{
	//return (my_date.alpha + my_date.beta * pow(x, my_date.gamma));
	return (my_date.kappa * pow(x,my_date.sigma));
}

double	P1(double t, Date my_data)
{
	// if (t > EPS && t < my_data.t0)
	// 	return my_data.Q;
	// else
		return 0;
}

double	P2(double t, Date my_data)
{
	if (EPS <= t && t < my_data.t0)
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

//	return my_data.u0*pow(x, 1/my_data.sigma);	//метода
	return my_data.u0;	//Ирин вариант
//	return 0;
}

double	u0_t(double x, Date my_data)
{
    //return my_data.u0 + x * (my_data.L - x);
    //return sin(x);

	return 0;	//метода
	//return my_data.u0;	//Ирин вариант
}

// правая прогонка
std::vector<double> progon(int n, std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> d, Date my_data)

{
	std::vector<double> y(n+1);
	std::vector<double> alfa(n+1);
	std::vector<double> betta(n+1);
	
	a[1] = 0; c[n-1] = 0;
	alfa[2] = c[1] / b[1]; betta[2] = d[1] / b[1];


	for (int i = 2; i < n; i++) {
		alfa[i+1] = c[i] / (b[i] - a[i] * alfa[i]);
		betta[i+1] = (a[i] * betta[i] + d[i]) / (-a[i] * alfa[i] + b[i]);
	}

	for (int i = n - 1; i > 0; i--) 
		y[i] = alfa[i+1] * y[i + 1] + betta[i+1];
	
	return  y;
}

double intergate(std::vector<double> y, double h){
    double sum = 0.0;
	int n = y.size();

    sum += 0.5*y[0] * h;
    for (int i = 1; i < n - 1; i++)
        sum += y[i]*h;
	sum += 0.5* y[n-1] * h;

    return sum;
}