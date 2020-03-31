#pragma once
#include "integro_interpolation.h"

int M = 3;

void quasilinear(int n, int t, double h, double tau, int TEST_P, Date my_date)
{
	std::ofstream		fout;
	std::vector<double>	a(n+1);
	std::vector<double> B(n);
	std::vector<double> A(n + 1);
	std::vector<double> C(n);
	std::vector<double> F(n);
	std::vector<double> y1(n + 1);	// �������� �� ������� ��������� ����
	std::vector<double> y2(n + 1);	// ��������� ��������� ����

	// ������������� ���������� �������
	for (int i = 0; i != n+1; i++)
		y1[i] = u0_t(i*h, my_date);

	for (int i = 1; i != n; i++)
		a[i] = (0.5*(K_quasi(y1[i], my_date) + K_quasi(y1[i-1], my_date)));

	fout.open("Quasilinear.txt");
	//������������ ��������
	A[0] = 0;
	for (int i = 1; i <= n; i++)
	{
		A[i] = a[i] / h;
		B[i - 1] = A[i];
		C[i - 1] = A[i - 1] + B[i - 1] + my_date.c* my_date.rho * h / tau;
	}
	std::cout << "Coefficients was found\n";

	//���������� �� ��������� �����
	for (double j = 0; j <= t; j += tau)
	{
		for (int i = 1; i != n; i++)
			a[i] = (0.5*(K_quasi(y1[i], my_date) + K_quasi(y1[i - 1], my_date)));
		//������������ ��������
		A[0] = 0;
		for (int i = 1; i <= n; i++)
		{
			A[i] = a[i] / h;
			B[i - 1] = A[i];
			C[i - 1] = A[i - 1] + B[i - 1] + my_date.c* my_date.rho * h / tau;
		}
		std::cout << "Coefficients was found\n";

		// ������������� ������� ������ �����
		for (int i = 1; i != n - 1; i++)
			F[i] = my_date.c * my_date.rho * y1[i] * h / tau;
		//�������� �������� � 1 �� n-1, ��� ��� ��� ��� ����������
		y2 = progon(A, C, B, F, n, my_date.left_boarder(0, my_date));
		//y2[0] = my_date.left_boarder(0, my_date);
		if ((TEST_P == 1) || (TEST_P == 3)) // TEST_P == 1 - ��������� ������
											// TEST_P == 3 - ��� ������ �� ������
		{
			std::cout << "progon is ready" << std::endl;
			double kappa = (a[n-1] / h) / (my_date.c * my_date.rho * h / tau + a[n] / h + a[n-1] / h);
			double mu = (my_date.c * my_date.rho * y1[n - 1] * h / tau)
				/ (my_date.c * my_date.rho * h / tau + a[n] / h + a[n-1] / h);
			std::cout << "mu1 = " << mu << std::endl;
			y2[n] = kappa * y2[n - 1] + mu;
			if (TEST_P == 3)//�� ������� ������, � �� �� ��������???????
			{
				kappa = (a[0] / h) / (my_date.c * my_date.rho * h / tau + a[1] / h + a[0] / h);
				mu = (my_date.c * my_date.rho * my_date.right_boarder(tau * (j - 1), my_date) * h / tau)
					/ (my_date.c * my_date.rho * h / tau + a[n + 1] / h + a[n] / h);
				y2[0] = kappa * y2[1] + mu;
				std::cout << "mu2 = " << mu << std::endl;
			}
		}
		else if (TEST_P == 2)
			y2[n] = my_date.left_boarder(n*h, my_date);

		for (int i = 0; i != y1.size(); i++)
			fout << j << '\t' << i * h << '\t' << y1[i] << '\n';
		y1 = y2;
		std::cout << j << std::endl;
	}
	fout.close();
}

void non_linear(int n, int t, double h, double tau, int TEST_P, Date my_date)
{
	std::ofstream		fout;
	std::vector<double>	a(n + 1);
	std::vector<double> B(n);
	std::vector<double> A(n + 1);
	std::vector<double> C(n);
	std::vector<double> F(n);
	std::vector<double> y1(n + 1);	// �������� �� ������� ��������� ����
	std::vector<double> y2(n + 1);	// �������� � ��������� ����
	std::vector<double> y_tmp(n + 1);	// ��������������� ������ ��� ���������� ��������

	// ������������� ���������� �������
	for (int i = 0; i <= n; i++)
		y1[i] = u0_t(i*h, my_date);

	for (int i = 1; i != n+1; i++)
		a[i] = (0.5*(K_quasi(y1[i], my_date) + K_quasi(y1[i - 1], my_date)));

	fout.open("Non_linear.txt");

	//������������ ��������
	A[0] = 0;
	for (int i = 1; i <= n; i++)
	{
		A[i] = a[i] / h;
		B[i - 1] = A[i];
		C[i - 1] = A[i - 1] + B[i - 1] + my_date.c* my_date.rho * h / tau;
	}
	std::cout << "Coefficients was found\n";

	//���������� �� ��������� �����
	for (double j = 0; j <= t; j += tau)
	{
		for (int k = 0; k != M; k++) {
			for (int i = 1; i != n+1; i++)
				a[i] = (0.5*(K_quasi(y1[i], my_date) + K_quasi(y1[i - 1], my_date)));

			//������������ ��������
			A[0] = 0;
			for (int i = 1; i <= n; i++)
			{
				A[i] = a[i] / h;
				B[i - 1] = A[i];
				C[i - 1] = A[i - 1] + B[i - 1] + my_date.c* my_date.rho * h / tau;
			}

			// ������������� ������� ������ �����
			for (int i = 1; i != n; i++)
				F[i] = my_date.c * my_date.rho * y1[i] * h / tau;
			//�������� �������� � 1 �� n-1, ��� ��� ��� ��� ����������
			y_tmp = progon(A, C, B, F, n, my_date.left_boarder(0, my_date));
			//y2[0] = my_date.left_boarder(0, my_date);
			if ((TEST_P == 1) || (TEST_P == 3)) // TEST_P == 1 - ��������� ������
												// TEST_P == 3 - ��� ������ �� ������
			{
				std::cout << "progon is ready" << std::endl;
				double kappa = (a[n - 1] / h) / (my_date.c * my_date.rho * h / tau + a[n] / h + a[n - 1] / h);
				double mu = (my_date.c * my_date.rho * y_tmp[n - 1] * h / tau)
					/ (my_date.c * my_date.rho * h / tau + a[n] / h + a[n - 1] / h);
				std::cout << "mu1 = " << mu << std::endl;
				y1[n] = kappa * y1[n - 1] + mu;
				if (TEST_P == 3)//�� ������� ������, � �� �� ��������???????
				{
					kappa = (a[0] / h) / (my_date.c * my_date.rho * h / tau + a[1] / h + a[0] / h);
					mu = (my_date.c * my_date.rho * my_date.right_boarder(tau * (j - 1), my_date) * h / tau)
						/ (my_date.c * my_date.rho * h / tau + a[n + 1] / h + a[n] / h);
					y1[0] = kappa * y1[1] + mu;
					std::cout << "mu2 = " << mu << std::endl;
				}
			}
			else if (TEST_P == 2)
				y_tmp[n] = my_date.left_boarder(n*h, my_date);
			y1 = y_tmp;
		}
		for (int i = 0; i != y1.size(); i++)
			fout << j << '\t' << i * h << '\t' << y1[i] << '\n';
	}
	fout.close();
}