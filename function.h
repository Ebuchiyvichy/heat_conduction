//
// Created by irina on 21.03.2020.
//

#ifndef HEAT_FUNCTION_H
#define HEAT_FUNCTION_H

#endif //HEAT_FUNCTION_H
#include "vector_new.h"

double K (double x, double x1, double x2, double L)
{
    double k1 =1;
    double k2 = 0.1;
    if (fabs(x - 0) <= EPS && fabs(x-x1)<= EPS)
        return k1;
    if (fabs(x-x1) > EPS && fabs(x-x2) > EPS)
        return (k1*(x-x2)/(x1-x2)+k2*(x-x1)/(x2-x1));
	if (fabs(x - x2) <= EPS && fabs(x - L) <= EPS)
		return k2;
}

double K (double u)
{
    double a = 2;
    double b = 0.5;
    double gamma = 3;
    return (a + b* pow(u, gamma));
}

double P(double t, double t0)
{
    double Q = 10;
    if ( t > EPS && t < t0)
        return Q;
}

double u0(double x, double L)
{
    double u0 = 0.04;
    return u0;
}

std::vector<double> progon(std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> f, int n)
{
    double m;
    std::vector<double> y(n);
    for (int i = 1; i < n; i++)
    {
        m = c[i] / b[i-1];
        b[i] -= m*a[i-1];
        f[i] = f[i] -f[i-1];
    }
    y[n-1] = f[n-1]/c[n-1];

    for (int i = n-2; i >= 1; i--)
        y[i] = (f[i] - c[i]*y[i+1])/b[i];
    return y;
}

void intergo_interpolation(int n, int t,  double h, double tau, double c, double p, double L)
{
    std::ofstream		fout;
	double	x(n);
	std::vector<double>	a(n);
	std::vector<double>	y(n);
	std::vector<double> A(n);
    std::vector<double> B(n-1);
    std::vector<double> C(n-1);
    std::vector<double> F(n);
    std::vector<double> y1(n); // значения на первом временном слое
    // инициализация начальными данными
    for (int i = 0; i != n; i++)
        y1[i] = u0(i*h, L);

    double sigma = 0.5;
    double kappa = (sigma*a[n-1]/h)/(c*p*h/(2*tau)+sigma*a[n-1]/n);
    for (int i = 0; i != n; i++)
        a.push_back(K(x + i*h - 0.5*h));
    fout.open("Integtgro_interpolation_mult.txt");

    for (int i = 1; i != n-1; i++)
    {
        A[i] = sigma/h * a[i];
        B[i] = sigma/h * a[i+1];
        C[i] = -sigma/h * (a[i] + a[i+1]) - c* p * h/tau;
    }
    A[n-1] = sigma / h *a[n-1]; A[0] = 0;

    std::vector<double> y2(n); // следующий временной слой
    for (int j = 0; j != t; j++)
    {
        // инициализация функции правой части
        for (int i = 1; i != n-1; i++)
            F[i] = -c*p*y1[i]*h/tau - (1-sigma)*a[i]*(y1[i+1] - 2 * y1[i] + y[i-1])/h;
        double mu = (c*p*y1[n-1] *h/(2*tau) + sigma * P(tau*j,tau)+(1-sigma)*(P(tau * (j-1), tau)-(y1[n-1] - y1[n-2])/h))/(c*p*h/(2*tau)+sigma * a[n]/h);

        //передача значений с 1 по n-1, так как они уже определены
        y2 = progon(A, C, B, F, n - 1);
        y2[n-1] = kappa * y2[n-2] + mu;
        y2[0] = u0(h, L);
     //   print_in_file(y2, fout);//??????????
        y1 = y2;
    }
    fout.close();
}

