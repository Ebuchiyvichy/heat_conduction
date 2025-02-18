//
// Created by irina on 21.03.2020.
//

#ifndef HEAT_VECTOR_NEW_H
#define HEAT_VECTOR_NEW_H

#endif //HEAT_VECTOR_NEW_H
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <functional>
#include <string>;

double	EPS = 1.e-3;
double PI = 3.1415926535;

struct Date
{
	//double	t0 = 1; //из методы
	//double	Q = 10;
	//double	L = 10;
	//double	u0 = 0.1;
	//double	alpha = 0;
	//double	beta = 0.5;
	//double	gamma = 2;
	//double	c = 1;
	//double	rho = 1;

	double	t0 = 0.5; //Ирин вариант
	double	Q = 10;
	double	L = 1;
	double	alpha = 0.5;
	double	beta = 2;
	double	gamma = 2;
	double	c = 1;
	double	rho = 1;
    double u0 = 0.1;

	double	k1 = 1;
	double	k2 = 0.1;
	double	x1 = 1.0 / 3;
	double	x2 = 2.0 / 3;
    double kappa = 0.5;
    double sigma = 2;
    //double u0 = pow(sigma*c*c/kappa, 1/sigma);
	std::function<double(double, Date)>	left_boarder;
	std::function<double(double, Date)>	right_boarder;
};

void print(std::vector<double> x)
{
    for (int i = 0; i!= x.size(); i++)
        std::cout << x[i] << '\t';
    std::cout << std::endl;
}

// переодпределение операций под вектора
std::vector<double> operator * (double a, std::vector<double> b)
{
    std::vector<double> c(b);
    for (int i = 0; i != c.size(); i++)
        c[i] = a * b[i];
    return c;
}
std::vector<double> operator + (std::vector<double> a, std::vector<double> b)
{
    std::vector<double> c(b);
    for (int i = 0; i != c.size(); i++)
        c[i] = a[i] + b[i];
    return c;
}
std::vector<double> operator - (std::vector<double> a, std::vector<double> b)
{
    std::vector<double> c(b);
    for (int i = 0; i != c.size(); i++)
        c[i] = a[i] - b[i];
    return c;
}
std::vector<double> operator / (std::vector<double> a, double b)
{
    std::vector<double> c(a);
    for (int i = 0; i != c.size(); i++)
        c[i] = a[i] / b;
    return c;
}

// нормы векторов
double norm(std::vector<double> x)
{
    double sum = 0;
    for(int i = 0; i != x.size(); i++)
        sum += (x[i]*x[i]);
    return sqrt(sum);
}

double	norm(std::vector<double> a, std::vector<double> b)
{
    double	max;
    max = fabs(a[0] - b[0]);
    for (int i = 1; i < a.size(); i++)
        if (fabs(a[i] - b[i]) > max)
            max = fabs(a[i] - b[i]);
    return max;
}