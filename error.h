//
// Created by irina on 04.04.2020.
//

#ifndef HEAT_ERROR_H
#define HEAT_ERROR_H

#endif //HEAT_ERROR_H

//
// Created by irina on 27.02.2020.
//

#ifndef DIFF_SOLVE_ERROR_H
#define DIFF_SOLVE_ERROR_H

#endif //DIFF_SOLVE_ERROR_H
#include "quasilinear.h"

void error_check(double h, double tau, double sigma, int TEST_P,  Date DATA)
{
    std::ofstream fout;
    double	temp;
    fout.open("Error_file.txt");
    fout << "Integrointerpolation \n";
    int dim = 2;
    double T = 0.05;
    std::cout << "Error is ready\n";
    std::vector<double> tmp1(DATA.L / h);
    std::vector<double> tmp2(DATA.L / (h/2));
    std::vector<double> tmp3(DATA.L / (h/4));
    tmp1 = integro_interpolation(DATA.L / h , T, h, tau, 2, 0.75, DATA);
    tmp2 = integro_interpolation(DATA.L * 2 / h , T, h / 2, tau / 4, 2, 0.75, DATA);
    double n = fabs(tmp1[10] -tmp2[10]);
    fout << n <<'\t';
    tmp3 = integro_interpolation(DATA.L * 4 / h , T, h / 4, tau / 16, 2, 0.75, DATA);
    fout << n / fabs(tmp2[10] - tmp3[10]) << '\t' << log( n/fabs(tmp2[10] -tmp3[10]))/log(2) << '\n';
    fout << fabs(tmp2[10] -tmp3[10]);
    fout.close();
}
