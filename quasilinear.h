#pragma once
#include "integro_interpolation.h"

int M = 4;
//правильно для 2го теста
std::vector<double> quasilinear(int n, double t, double h, double tau, int TEST_P, Date my_date)
{
	std::ofstream		fout;
	std::vector<double>	a(n+1);
	std::vector<double> B(n);
	std::vector<double> A(n + 1);
	std::vector<double> C(n);
	std::vector<double> F(n);
	std::vector<double> y1(n + 1);	// значения на текущем временном слое
	std::vector<double> y2(n + 1);	// следующий временной слой

	// инициализация начальными данными
	for (int i = 0; i != n + 1; i++)
		y1[i] = u0_t(i * h, my_date);
	fout.open("quasilinear.txt");
	//вычисление по временным слоям
	for (double j = 0; j <= t; j += tau)
	{
		for (int i = 1; i != n + 1; i++)
			a[i] = 0.5 * (K_quasi(y1[i], my_date) + K_quasi(y1[i - 1], my_date));
		for (int i = 1; i <= n; i++)
		{
			A[i] = a[i] / h;
			B[i - 1] = A[i];
			C[i - 1] = A[i - 1] + B[i - 1] + my_date.c * my_date.rho * h / tau;
		}
		for (int i = 1; i != n; i++)
			F[i] = my_date.c * my_date.rho * y1[i] * h / tau;

		if ((TEST_P == 1) || (TEST_P == 3)) // TEST_P == 1 - смешанная задача
											// TEST_P == 3 - два потока на концах
		{
			double kappa = (a[n] / h) / (my_date.c * my_date.rho * h / (2 * tau) + a[n] / h);
			double mu = (my_date.c * my_date.rho * y1[n] * h / (2 * tau) + my_date.right_boarder(j, my_date))
				/ (my_date.c * my_date.rho * h / (2 * tau) + a[n] / h);
			y2 = progon(A, C, B, F, n, my_date.left_boarder(j, my_date), kappa, mu);
			y2[n] = kappa * y2[n - 1] + mu;
			if (TEST_P == 3 )
			{
				kappa = (a[0] / h) / (my_date.c * my_date.rho * h / tau + a[1] / h + a[0] / h);
				mu = (my_date.c * my_date.rho * my_date.right_boarder(j - tau, my_date) * h / tau + my_date.left_boarder(j, my_date))
					/ (my_date.c * my_date.rho * h / tau + a[n + 1] / h + a[n] / h);
				y2[0] = kappa * y2[1] + mu;
			}
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

std::vector<double> non_linear(int n, double t, double h, double tau, int TEST_P, Date my_date)
{
	std::ofstream		fout;
	std::vector<double>	a(n + 1);
	std::vector<double> B(n);
	std::vector<double> A(n + 1);
	std::vector<double> C(n);
	std::vector<double> F(n);
	std::vector<double> y1(n + 1);	// значения на текущем временном слое
	std::vector<double> y2(n + 1);	// значения на следующем временном слое
	std::vector<double> y_tmp(n + 1);	// вспомогательный вектор для выполнения итераций

	// инициализация начальными данными
	for (int i = 0; i <= n; i++)
		y2[i] = u0_t(i * h, my_date);
	y1 = y2;

	fout.open("Non_linear.txt");
	//вычисление по временным слоям
	for (double j = 0; j <= t; j += tau)
	{
		for (int k = 0; k != M; k++) {
			for (int i = 1; i != n + 1; i++)
				a[i] = (0.5 * (K_quasi(y1[i], my_date) + K_quasi(y1[i - 1], my_date)));
			for (int i = 1; i <= n; i++)
			{
				A[i] = a[i] / h;
				B[i - 1] = A[i];
				C[i - 1] = A[i - 1] + B[i - 1] + my_date.c * my_date.rho * h / tau;
			}
			// инициализация функции правой части
			for (int i = 1; i != n; i++)
				F[i] = my_date.c * my_date.rho * y2[i] * h / tau;
			if ((TEST_P == 1) || (TEST_P == 3)) // TEST_P == 1 - смешанная задача
												// TEST_P == 3 - два потока на концах
			{
				double kappa = (a[n] / h) / (my_date.c * my_date.rho * h / (2 * tau) + a[n] / h);
				double mu = (my_date.c * my_date.rho * y2[n] * h / (2 * tau) + my_date.right_boarder(j, my_date))
					/ (my_date.c * my_date.rho * h / (2 * tau) + a[n] / h);
				y_tmp = progon(A, C, B, F, n, my_date.left_boarder(j, my_date), kappa, mu);
				y_tmp[n] = kappa * y1[n - 1] + mu;
				if (TEST_P == 3)//ИЗ УСЛОВИЙ ПОТОКА, А НЕ ИЗ ПРОГОНКИ???????
				{
					kappa = (a[0] / h) / (my_date.c * my_date.rho * h / tau + a[1] / h + a[0] / h);
					mu = (my_date.c * my_date.rho * my_date.right_boarder(tau * (j - 1), my_date) * h / tau + my_date.right_boarder(tau * j, my_date))
						/ (my_date.c * my_date.rho * h / tau + a[n + 1] / h + a[n] / h);
					y1[0] = kappa * y1[1] + mu;
				}
			}
			else if (TEST_P == 2)
			{
				y_tmp = progon(A, C, B, F, n, my_date.left_boarder(j, my_date), 0, 0);
				y_tmp[n] = my_date.right_boarder(j, my_date);
			}
			y1 = y_tmp;
		}
		for (int i = 0; i != y2.size(); i++)
			fout << j << '\t' << i * h << '\t' << y1[i] << '\n';
		y2 = y1;
	}
	fout.close();
	return y1;
}

int non_linear_iter(int n, double t, double h, double tau, int TEST_P, Date my_date)
{
	std::ofstream		fout;
	std::vector<double>	a(n + 1);
	std::vector<double> B(n);
	std::vector<double> A(n + 1);
	std::vector<double> C(n);
	std::vector<double> F(n);
	std::vector<double> y1(n + 1);	// значения на текущем временном слое
	std::vector<double> y2(n + 1);	// значения на следующем временном слое
	std::vector<double> y_tmp(n + 1);	// вспомогательный вектор для выполнения итераций

	// инициализация начальными данными
	for (int i = 0; i <= n; i++)
		y2[i] = u0_t(i * h, my_date);
	y1 = y2;

	fout.open("Non_linear.txt");
	//вычисление по временным слоям
	int iter = 0;
	for (double j = 0; j <= t; j += tau)
	{
		do {
			y1 = y2;
			for (int i = 1; i != n + 1; i++)
				a[i] = (0.5 * (K_quasi(y1[i], my_date) + K_quasi(y1[i - 1], my_date)));
			for (int i = 1; i <= n; i++)
			{
				A[i] = a[i] / h;
				B[i - 1] = A[i];
				C[i - 1] = A[i - 1] + B[i - 1] + my_date.c * my_date.rho * h / tau;
			}
			// инициализация функции правой части
			for (int i = 1; i != n; i++)
				F[i] = my_date.c * my_date.rho * y2[i] * h / tau;
			if (TEST_P == 2)
			{
				F[0] -= A[0]*my_date.left_boarder(j, my_date);
				F[n-1] -= B[n-1]*my_date.right_boarder(j, my_date);
				y2 = progon(A, C, B, F, n, my_date.left_boarder(j, my_date), 0, 0);
				y2[n] = my_date.right_boarder(j, my_date);
			}
			else {
				std::cout << "Неправильный номер теста для подсчета итераций";
				break;
			}
			iter++;
		} while (norm(y1,y2) > EPS);
	}
	fout.close();
	return iter;
}
