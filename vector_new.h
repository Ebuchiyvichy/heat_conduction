//
// Created by irina on 21.03.2020.
//

#ifndef HEAT_VECTOR_NEW_H
#define HEAT_VECTOR_NEW_H

#endif //HEAT_VECTOR_NEW_H
#include <vector>
#include <math>
#include <iostream>
#include <ofstream>
void print(std::vector<double> x)
{
    for (int i = 0; i!= x.size(); i++)
        cout << x[i] << '\t';
    cout << std::endl;
}

void print_in_file(std::vector<double> x, ofstream fout)
{
    for (int i = 0; i!= x.size(); i++)
        fout << x[i] << '\t';
    fout << std::endl;
}

std::vector<double>	cpy_vector(std::vector<double> tmp, std::vector<double> x)
{

    for (int i = 0; i < x.size(); i++)
        tmp[i] = x[i];
    return (tmp);
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
    for (int i = 1; i < dim; i++)
        if (fabs(a[i] - b[i]) > max)
            max = fabs(a[i] - b[i]);
    return max;
}