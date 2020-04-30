#pragma once
#include "integro_interpolation.h"

int M = 4;
//ïðàâèëüíî äëÿ 2ãî òåñòà
std::vector<double> quasilinear(int n, double t, double h, double tau, int TEST_P, Date my_date)
{
	std::ofstream		fout;
	std::vector<double>	a(n+1);
	std::vector<double> B(n);
	std::vector<double> A(n + 1);
	std::vector<double> C(n);
	std::vector<double> F(n);
	std::vector<double> y1(n + 1);	// çíà÷åíèÿ íà òåêóùåì âðåìåííîì ñëîå
	std::vector<double> y2(n + 1);	// ñëåäóþùèé âðåìåííîé ñëîé

	// èíèöèàëèçàöèÿ íà÷àëüíûìè äàííûìè
	for (int i = 0; i != n + 1; i++)
		y1[i] = u0_t(i * h, my_date);

	fout.open("quasilinear.txt");
	//âû÷èñëåíèå ïî âðåìåííûì ñëîÿì
	for (double j = tau; j <= t; j += tau)
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

		if ((TEST_P == 1) || (TEST_P == 3)) // TEST_P == 1 - ñìåøàííàÿ çàäà÷à
											// TEST_P == 3 - äâà ïîòîêà íà êîíöàõ
		{
			double kappa = (a[n] / h) / (my_date.c * my_date.rho * h / (2 * tau) + a[n] / h);
			double mu = (my_date.c * my_date.rho * y1[n] * h / (2 * tau) + my_date.right_boarder(j, my_date))
				/ (my_date.c * my_date.rho * h / (2 * tau) + a[n] / h);
			y2 = progon(A, C, B, F, n, my_date.left_boarder(j, my_date), 0, kappa, mu);
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
			y2 = progon(A, C, B, F, n, my_date.left_boarder(j, my_date), my_date.right_boarder(j, my_date), 0, 0);
			y2[n] = my_date.right_boarder(j, my_date);
		}
		y1 = y2;
	}
	for (int i = 0; i != y2.size()-1; i++)
			fout << i * h << '\t' << y1[i] << '\n';
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
	std::vector<double> y1(n + 1);	// çíà÷åíèÿ íà òåêóùåì âðåìåííîì ñëîå
	std::vector<double> y2(n + 1);	// çíà÷åíèÿ íà ñëåäóþùåì âðåìåííîì ñëîå
	std::vector<double> y_tmp(n + 1);	// âñïîìîãàòåëüíûé âåêòîð äëÿ âûïîëíåíèÿ èòåðàöèé

	// èíèöèàëèçàöèÿ íà÷àëüíûìè äàííûìè
	for (int i = 0; i <= n; i++)
		y2[i] = u0_t(i * h, my_date);
	y1 = y2;

	fout.open("Non_linear.txt");
	//âû÷èñëåíèå ïî âðåìåííûì ñëîÿì
	for (double j = tau; j <= t; j += tau)
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
			// èíèöèàëèçàöèÿ ôóíêöèè ïðàâîé ÷àñòè
			for (int i = 1; i != n; i++)
				F[i] = my_date.c * my_date.rho * y2[i] * h / tau;
			if ((TEST_P == 1) || (TEST_P == 3)) // TEST_P == 1 - ñìåøàííàÿ çàäà÷à
												// TEST_P == 3 - äâà ïîòîêà íà êîíöàõ
			{
				double kappa = (a[n] / h) / (my_date.c * my_date.rho * h / (2 * tau) + a[n] / h);
				double mu = (my_date.c * my_date.rho * y2[n] * h / (2 * tau) + my_date.right_boarder(j + tau, my_date))
					/ (my_date.c * my_date.rho * h / (2 * tau) + a[n] / h);
				y_tmp = progon(A, C, B, F, n-1, my_date.left_boarder(j, my_date), 0, kappa, mu);
				y_tmp[n] = kappa * y_tmp[n - 1] + mu;
				if (TEST_P == 3)//ÈÇ ÓÑËÎÂÈÉ ÏÎÒÎÊÀ, À ÍÅ ÈÇ ÏÐÎÃÎÍÊÈ???????
				{
					kappa = (a[0] / h) / (my_date.c * my_date.rho * h / tau + a[1] / h + a[0] / h);
					mu = (my_date.c * my_date.rho * my_date.right_boarder(tau * (j - 1), my_date) * h / tau + my_date.right_boarder(tau * j, my_date))
						/ (my_date.c * my_date.rho * h / tau + a[n + 1] / h + a[n] / h);
					y1[0] = kappa * y1[1] + mu;
				}
			}
			else if (TEST_P == 2)
			{
				y_tmp = progon(A, C, B, F, n, my_date.left_boarder(j, my_date), my_date.right_boarder(j, my_date), 0, 0);
				y_tmp[n] = my_date.right_boarder(j, my_date);
			}
			y1 = y_tmp;
		}
		
		y2 = y1;
	}
	for (int i = 0; i != y2.size(); i++)
			fout << i * h << '\t' << y1[i] << '\n';
	fout.close();
	return y1;
}

std::vector<double> non_linear_iter(int n, double t, double h, double tau, int TEST_P, Date my_date)
{
	std::ofstream		fout;
	std::vector<double>	a(n + 1);
	std::vector<double> B(n);
	std::vector<double> A(n + 1);
	std::vector<double> C(n);
	std::vector<double> F(n);
	std::vector<double> y1(n + 1);	// �������� �� ������� ��������� ����
	std::vector<double> y2(n + 1);	// �������� �� ��������� ��������� ����
	std::vector<double> y_tmp(n + 1);	// ��������������� ������ ��� ���������� ��������
	int iter;
	// ������������� ���������� �������
	for (int i = 0; i <= n; i++)
		y2[i] = u0_t(i * h, my_date);
	y1 = y2;
	y_tmp = y1;

	fout.open("Non_linear_iter.txt");
	//���������� �� ��������� �����
	for (double j = tau; j <= t; j += tau)
	{
		iter = 0;
		do {
			iter++;
			y1 = y_tmp;
			for (int i = 1; i != n + 1; i++)
				a[i] = (0.5 * (K_quasi(y1[i], my_date) + K_quasi(y1[i - 1], my_date)));
			for (int i = 1; i <= n; i++)
			{
				A[i] = a[i] / h;
				B[i - 1] = A[i];
				C[i - 1] = A[i - 1] + B[i - 1] + my_date.c * my_date.rho * h / tau;
			}
			// ������������� ������� ������ �����
			for (int i = 1; i != n; i++)
				F[i] = my_date.c * my_date.rho * y2[i] * h / tau;
			if ((TEST_P == 1) || (TEST_P == 3)) // TEST_P == 1 - ��������� ������
												// TEST_P == 3 - ��� ������ �� ������
			{
				double kappa = (a[n] / h) / (my_date.c * my_date.rho * h / (2 * tau) + a[n] / h);
				double mu = (my_date.c * my_date.rho * y2[n] * h / (2 * tau) + my_date.right_boarder(j, my_date))
					/ (my_date.c * my_date.rho * h / (2 * tau) + a[n] / h);
				y_tmp = progon(A, C, B, F, n, my_date.left_boarder(j, my_date),my_date.right_boarder(j,my_date), kappa, mu);
				y_tmp[n] = kappa * y1[n - 1] + mu;
				if (TEST_P == 3)//�� ������� ������, � �� �� ��������???????
				{
					kappa = (a[0] / h) / (my_date.c * my_date.rho * h / tau + a[1] / h + a[0] / h);
					mu = (my_date.c * my_date.rho * my_date.right_boarder(tau * (j - 1), my_date) * h / tau + my_date.right_boarder(tau * j, my_date))
						/ (my_date.c * my_date.rho * h / tau + a[n + 1] / h + a[n] / h);
					y1[0] = kappa * y1[1] + mu;
				}
			}
			else if (TEST_P == 2)
			{
				y_tmp = progon(A, C, B, F, n, my_date.left_boarder(j, my_date),0, 0, 0);
				y_tmp[n] = my_date.right_boarder(j, my_date);
			}
		} while (norm(y1,y_tmp)> EPS);
		std::cout << iter << "\n";
		fout << j << '\t' <<  iter << "\n";
		y2 = y1;
	}

	fout.close();
	return y1;
}