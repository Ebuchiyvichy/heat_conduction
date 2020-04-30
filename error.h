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
    fout << "\n Integrointerpolation \n";
    int dim = 2;
    double T = tau;
    std::cout << "Error is ready\n";
    double f0 = fabs(integro_interpolation(DATA.L / h , T, h, tau, TEST_P, sigma, DATA)[10] -
                     integro_interpolation(DATA.L * 2 / h , T, h / (2), tau / (4), TEST_P, sigma, DATA)[10 * 2]);

    if (true)//(fabs(sigma - 0.5) >= EPS)
    {
        for (int i = 1; i <= 3; i++)
        {
            fout <<  f0 << '\t';

            double f1 = fabs(integro_interpolation(DATA.L * pow(2, i) / h , T, h / pow(2, i), tau / pow(4, i), TEST_P, sigma, DATA)[10 * pow(2, i)] -
                    integro_interpolation(DATA.L * pow(2, i + 1) / h , T, h / pow(2, i + 1), tau / pow(4, i + 1), TEST_P, sigma, DATA)[10 * pow(2, i + 1)]);
            fout << f0 / f1 << '\t' << log(f0 / f1) << '\n';
            f0 = f1;
        }
    }
 /*   fout << "\n Quasi \n";
    std::cout << "Error is ready\n";
    f0 = fabs(quasilinear(DATA.L   / h , T, h, tau, TEST_P, DATA)[10] -
                      quasilinear(DATA.L * 2 / h , T, h / (2), tau / (2), TEST_P, DATA)[10 * 2]);

    if (true)//(fabs(sigma - 0.5) >= EPS)
    {
        for (int i = 1; i <= 8; i += 2)
        {
            fout <<  f0 << '\t';

            double f1 = fabs(quasilinear(DATA.L * 2 * i / h , T, h / (2 * i), tau / (2 * i), TEST_P, DATA)[10 * 2 * i] -
                                     quasilinear(DATA.L * 4 * i / h , T, h / (4 * i), tau / (4 * i), TEST_P, DATA)[10 * 4 * i]);
            fout << f0/f1 << '\t' << log(f0)/log(f1) << '\n';
            f0 = f1;
        }
    }
    fout << "\n Nonlinear \n";
    std::cout << "Error is ready\n";
    f0 = fabs(non_linear(DATA.L   / h , T, h, tau, TEST_P, DATA)[10] -
                      non_linear(DATA.L * 2 / h , T, h / (2), tau / (2), TEST_P, DATA)[10 * 2]);

    if (true)//(fabs(sigma - 0.5) >= EPS)
    {
        for (int i = 1; i <= 8; i += 2)
        {
            fout <<  f0 << '\t';

            double f1 = fabs(non_linear(DATA.L * 2 * i / h , T, h / (2 * i), tau / (2 * i), TEST_P, DATA)[10 * 2 * i] -
                                     non_linear(DATA.L * 4 * i / h , T, h / (4 * i), tau / (4 * i), TEST_P, DATA)[10 * 4 * i]);
            fout << f0/f1 << '\t' << log(f0)/log(f1) << '\n';
            f0 = f1;
        }
    }
	*/

    fout.close();

}
