//
// Created by irina on 28.03.2020.
//

#ifndef HEAT_INTEGRO_INTERPOLATION_H
#define HEAT_INTEGRO_INTERPOLATION_H

#endif //HEAT_INTEGRO_INTERPOLATION_H

#include "function.h"

std::vector<double> integro_interpolation(int n, double T,  double h, double tau, int TEST_P, double sigma, Date my_date)
{
    int N = n-1;
    std::ofstream		fout;
    std::vector<double>	a(N);
    std::vector<double> B(N-1);
    std::vector<double> A(N-1);
    std::vector<double> C(N-1);
    std::vector<double> F(N-1);
    std::vector<double> y1(n);	// значения на текущем временном слое
    std::vector<double> y2(n);	// следующий временной слой
	std::string			str = "Integro_interpolation_mult_";

	str += std::to_string(sigma) + ".txt";

    // инициализация начальными данными
	for (int i = 0; i < n; i++)
		y1[i] = u0_t(i*h, my_date);
	for (int i = 0; i < N; i++)
		a[i] = K(i * h - 0.5 * h, my_date);

    fout.open(str);

    //коэффициенты прогонки
    for (int i = 0; i < N-1; i++)
    {
        A[i] = sigma / h * a[i];
        B[i] = sigma / h * a[i+1];
        C[i] = - (A[i] + B[i] + my_date.c * my_date.rho * h / tau);
    }
    std::cout << "Coefficients was found\n";
    //вычисление по временным слоям
    for (double j = tau; j <= 2*tau; j += tau)
    {
        // инициализация функции правой части
        for (int i = 0; i < N-1; i++)
            F[i] = -(my_date.c * my_date.rho * y1[i+1] * h / tau + (1 - sigma) * a[i] * (y1[i + 2] - 2 * y1[i+1] + y1[i]) / h);
        if ((TEST_P == 1) || (TEST_P == 3)) // TEST_P == 1 - смешанная задача
                                            // TEST_P == 3 - два потока на концах
        {
            double kappa2 = (sigma * a[N-1] / h) / (my_date.c * my_date.rho * h / (2 * tau) + sigma * a[N-1] / h);
			double mu2 = (my_date.c * my_date.rho * y1[N] * h / (2 * tau) + sigma * my_date.right_boarder(j, my_date) +
				(1 - sigma) * (my_date.right_boarder((j - tau), my_date) - (y1[N] - y1[N - 1]) / h))
				/ (my_date.c * my_date.rho * h / (2 * tau) + sigma * a[N-1] / h);
			
            
            //y2 = progon(A, C, B, F, n, my_date.left_boarder(j, my_date), kappa, mu);
			
            if (TEST_P == 3)//из условий потока(?)
            {
                
                double kappa1 = (sigma * a[0] / h) / (my_date.c * my_date.rho * h / (2 * tau) + sigma * a[0] / h);
                double mu1 = (my_date.c * my_date.rho * y1[0] * h / (2 * tau) + sigma * my_date.left_boarder(j, my_date) +
                              (1 - sigma) * (my_date.left_boarder(j - tau, my_date) + (y1[1] - y1[0]) / h))
                             / (my_date.c * my_date.rho * h / (2 * tau) + sigma * a[0] / h);
                C[N-2] += B[N-2]*kappa2; F[N-2] -= B[N-2]*mu2;
                C[0] += A[0]*kappa1; F[0] -= A[0]*mu1;
                A[0] = 0; B[N-2] = 0; 
                y2 = convert_dim(y2,progon(A,C,B,F,N-1));
                y2[0] = kappa1 * y2[1] + mu1;
                y2[N] = kappa2 * y2[N-1] + mu2;
            }
            else if(TEST_P == 1)
            {
                F[0] -= A[0]* my_date.left_boarder(j, my_date);
                F[N-2] -= B[N-2]*mu2; C[N-2] += B[N-2]*kappa2;
                A[0] = 0; B[N-2] = 0; 
                y2 = convert_dim(y2,progon(A,C,B,F,N-1));
                y2[0] = my_date.left_boarder(j, my_date);
                y2[N] = kappa2 * y2[N-1] + mu2;
            }
        }
		else if (TEST_P == 2)
		{
            F[0] -= A[0]* my_date.left_boarder(j, my_date);
            F[N-2] -= B[N-2]*my_date.right_boarder(j,my_date);
            A[0] = 0; B[N-2] = 0;
			y2 = convert_dim(y2,progon(A,C,B,F,N-1));
            y2[0] = my_date.left_boarder(j, my_date);
			y2[N] = my_date.right_boarder(j, my_date);
		}
        for (int i = 0; i < n; i++)
            fout << j-tau << '\t' << i * h << '\t' << y1[i] << '\n';
        y1 = y2;
    }
    fout.close();
    return y1;
}

std::vector<double> new_test(int n, double t, double h, double tau, int TEST_P, Date my_date)
{
	std::ofstream		fout;
	std::vector<double>	a(n+1);
	std::vector<double> B(n);
	std::vector<double> A(n + 1);
	std::vector<double> C(n);
	std::vector<double> F(n);
	std::vector<double> y1(n + 1);	// значения на текущем временном слое
	std::vector<double> y2(n + 1);	// следующий временной слой

    double sigma = 0.5;
	// инициализация начальными данными
	for (int i = 0; i != n + 1; i++)
		y1[i] = u0_t(i * h, my_date);
	fout.open("test.txt");
	//вычисление по временным слоям
	for (double j = tau; j <= t; j += tau)
	{
		for (int i = 1; i != n + 1; i++)
			a[i] = 0.5 * (K_quasi(y1[i], my_date) + K_quasi(y1[i - 1], my_date));
		for (int i = 1; i <= n; i++)
		{
			A[i] = sigma / h * a[i];;
			B[i - 1] = A[i];
			C[i - 1] = A[i - 1] + B[i - 1] + my_date.c * my_date.rho * h / tau;
		}
		for (int i = 1; i != n; i++)
			F[i] = my_date.c * my_date.rho * y1[i+1] * h / tau + (1 - sigma) * a[i] * (y1[i + 2] - 2 * y1[i+1] + y1[i]) / h;

		if (TEST_P == 1) // TEST_P == 1 - смешанная задача
											// TEST_P == 3 - два потока на концах
		{
           double kappa = (sigma * a[n-1] / h) / (my_date.c * my_date.rho * h / (2 * tau) + sigma * a[n-1] / h);
			double mu= (my_date.c * my_date.rho * y1[n] * h / (2 * tau) + sigma * my_date.right_boarder(j, my_date) +
				(1 - sigma) * (my_date.right_boarder((j - tau), my_date) - (y1[n] - y1[n- 1]) / h))
				/ (my_date.c * my_date.rho * h / (2 * tau) + sigma * a[n-1] / h);
            y2 = convert_dim(y2, progon(A, C, B, F, n-1, my_date.left_boarder(j, my_date), 0, 0));
            y2[n] = kappa* y2[n-1] + mu;
			
		}

		else if (TEST_P == 2)
		{
			y2 = progon(A, C, B, F, n, my_date.left_boarder(j, my_date), 0, 0);
			y2[n] = my_date.right_boarder(j, my_date);
		}

		for (int i = 0; i != y1.size(); i++)
			fout << j << '\t' << i * h << '\t' << y1[i] << '\n';
		y1 = y2;
	}
	fout.close();
	return y1;
}