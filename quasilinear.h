#pragma once
#include "integro_interpolation.h"

// количество внутренних итераций
int M = 4;
std::vector<double> quasilinear(int n, double t, double h, double tau, int TEST_P, Date my_date)
{
	std::ofstream		fout;
	std::vector<double>	a(n+1);
	std::vector<double> B(n);
	std::vector<double> A(n + 1);
	std::vector<double> C(n);
	std::vector<double> F(n);
	std::vector<double> y1(n + 1);	
	std::vector<double> y2(n + 1);	
	std::vector<double> kappa;
	std::vector<double> mu;

	//заполнение первого временного слоя
	for (int i = 0; i != n + 1; i++)
		y1[i] = u0_t(i * h, my_date);

	fout.open("quasilinear.txt");
	
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

	if (TEST_P == 1 || TEST_P == 3){
        double kappa2 =  (a[n] / h) / (my_date.c * my_date.rho * h / (2 * tau) + a[n] / h);
        kappa.push_back(kappa2);
        mu.push_back(0);
    }
    if (TEST_P == 3)
    {
        double kappa1 = (a[0] / h) / (my_date.c * my_date.rho * h / tau + a[1] / h + a[0] / h);
        kappa.push_back(kappa1);
        mu.push_back(0);
    }
	if ((TEST_P == 1) || (TEST_P == 3)) 
	{
		mu[0] = (my_date.c * my_date.rho * y1[n] * h / (2 * tau) + my_date.right_boarder(j, my_date))
				/ (my_date.c * my_date.rho * h / (2 * tau) + a[n] / h);
			if (TEST_P == 1)
			{
				 F[n-1] += B[n-1] * mu[0]; 
	            C[n-1] -= B[n-1] * kappa[0];    
				y2[n] = my_date.right_boarder(0, my_date);
				y2 = progon(n,A, C, B, F, my_date);
            	y2[0] = kappa[0]* y2[n-1] + mu[0];
			}

			if (TEST_P == 3)
			{
				mu[1] = (my_date.c * my_date.rho * my_date.right_boarder(j - tau, my_date) * h / tau + my_date.left_boarder(j, my_date))
					/ (my_date.c * my_date.rho * h / tau + a[n + 1] / h + a[n] / h);
				C[1] -= A[1]*kappa[1];
	            C[n-1] -= B[n-1]*kappa[0];
	            F[1] += A[1]*mu[1]; F[n-1] += B[n-1]*mu[0];
                y2 = progon(n, A, C, B, F, my_date);
                y2[0] = kappa[1]*y2[1] + mu[1];
	            y2[n] = kappa[0]* y2[n-1] + mu[0];
			}
		}
		else if (TEST_P == 2)
		{
			F[1] += A[1]* my_date.left_boarder(0, my_date);
		    F[n-1] += B[n-1] * my_date.right_boarder(my_date.L, my_date);
			y2[0] = my_date.left_boarder(0, my_date);
            y2[n] = my_date.right_boarder(my_date.L, my_date);
			y2 = progon(n,A, C, B, F,my_date);
		}
		for (int i = 0; i != y2.size()-1; i++)
			fout << j-tau<< '\t' << i * h << '\t' << y1[i] << '\n';
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
	std::vector<double> y1(n + 1);
	std::vector<double> y2(n + 1);
	std::vector<double> y_tmp(n + 1);
	std::vector<double> kappa, mu;

	for (int i = 0; i <= n; i++)
		y2[i] = u0_t(i * h, my_date);
	y1 = y2;

	fout.open("Non_linear.txt");

	for (double j = tau; j <= t; j += tau)
	{
		// внутренние итерации
		for (int k = 0; k != M; k++) {
			for (int i = 1; i != n + 1; i++)
				a[i] = (0.5 * (K_quasi(y1[i], my_date) + K_quasi(y1[i - 1], my_date)));
			for (int i = 1; i <= n; i++)
			{
				A[i] = a[i] / h;
				B[i - 1] = A[i];
				C[i - 1] = A[i - 1] + B[i - 1] + my_date.c * my_date.rho * h / tau;
			}

			for (int i = 1; i != n; i++)
				F[i] = my_date.c * my_date.rho * y2[i] * h / tau;
			if (TEST_P == 1 || TEST_P == 3){
       			double kappa2 =  (a[n] / h) / (my_date.c * my_date.rho * h / (2 * tau) + a[n] / h);
        		kappa.push_back(kappa2);
        		mu.push_back(0);
    		}
    		if (TEST_P == 3)
    		{
       			double kappa1 = (a[0] / h) / (my_date.c * my_date.rho * h / tau + a[1] / h + a[0] / h);
        		kappa.push_back(kappa1);
        		mu.push_back(0);
    		}

			if ((TEST_P == 1) || (TEST_P == 3)) 
			{
				mu[0] = (my_date.c * my_date.rho * y1[n] * h / (2 * tau) + my_date.right_boarder(j, my_date))
				/ (my_date.c * my_date.rho * h / (2 * tau) + a[n] / h);
				if (TEST_P == 1)
				{
					F[n-1] += B[n-1] * mu[0]; 
	            	C[n-1] -= B[n-1] * kappa[0]; 
					y2[n] = my_date.right_boarder(0, my_date);
					y_tmp = progon(n, A,C, B, F, my_date);
					y_tmp[n] = kappa[0] * y_tmp[n - 1] + mu[0];
				}
				if (TEST_P == 3){	
					mu[1] = (my_date.c * my_date.rho * my_date.right_boarder(j - tau, my_date) * h / tau + my_date.left_boarder(j, my_date))
					/ (my_date.c * my_date.rho * h / tau + a[n + 1] / h + a[n] / h);
					C[1] -= A[1]*kappa[1];
	            	C[n-1] -= B[n-1]*kappa[0];
	            	F[1] += A[1]*mu[1]; F[n-1] += B[n-1]*mu[0];
					y_tmp = progon(n, A,C, B, F, my_date);
					y_tmp[0] = kappa[1] * y_tmp[1] + mu[1];
					y_tmp[n] = kappa[0]*  y_tmp[n-1] + mu[0];
				}
			}
			else if (TEST_P == 2)
			{
				F[1] += A[1]* my_date.left_boarder(0, my_date);
		   		F[n-1] += B[n-1] * my_date.right_boarder(my_date.L, my_date);
				y_tmp[0] = my_date.left_boarder(0, my_date);
           		y_tmp[n] = my_date.right_boarder(my_date.L, my_date);
				y_tmp = progon(n,A, C, B, F,my_date);
			}
			y1 = y_tmp;
		}
		for (int i = 0; i != y2.size()-1; i++)
			fout << j-tau << "\t" << i * h << '\t' << y2[i] << '\n';
		y2 = y1;
	}
	fout.close();
	return y1;
}

//функция подсчета количества итераций
std::vector<double> non_linear_iter(int n, double t, double h, double tau, int TEST_P, Date my_date)
{
	std::ofstream		fout;
	std::vector<double>	a(n + 1);
	std::vector<double> B(n);
	std::vector<double> A(n + 1);
	std::vector<double> C(n);
	std::vector<double> F(n);
	std::vector<double> y1(n + 1);	
	std::vector<double> y2(n + 1);
	std::vector<double> y_tmp(n + 1);
	std::vector<double> kappa, mu;
	int iter;

	for (int i = 0; i <= n; i++)
		y2[i] = u0_t(i * h, my_date);
	y1 = y2;
	y_tmp = y1;

	fout.open("Non_linear_iter.txt");
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

			for (int i = 1; i != n; i++)
				F[i] = my_date.c * my_date.rho * y2[i] * h / tau;
			if (TEST_P == 1 || TEST_P == 3){
       			double kappa2 =  (a[n] / h) / (my_date.c * my_date.rho * h / (2 * tau) + a[n] / h);
        		kappa.push_back(kappa2);
        		mu.push_back(0);
    		}
    		if (TEST_P == 3)
    		{
       			double kappa1 = (a[0] / h) / (my_date.c * my_date.rho * h / tau + a[1] / h + a[0] / h);
        		kappa.push_back(kappa1);
        		mu.push_back(0);
    		}

			if ((TEST_P == 1) || (TEST_P == 3)) 
			{
				mu[0] = (my_date.c * my_date.rho * y1[n] * h / (2 * tau) + my_date.right_boarder(j, my_date))
				/ (my_date.c * my_date.rho * h / (2 * tau) + a[n] / h);
				if (TEST_P == 1)
				{
					F[n-1] += B[n-1] * mu[0]; 
	            	C[n-1] -= B[n-1] * kappa[0]; 
					y2[n] = my_date.right_boarder(0, my_date);
					y_tmp = progon(n, A,C, B, F, my_date);
					y_tmp[n] = kappa[0] * y_tmp[n - 1] + mu[0];
				}
				if (TEST_P == 3){	
					mu[1] = (my_date.c * my_date.rho * my_date.right_boarder(j - tau, my_date) * h / tau + my_date.left_boarder(j, my_date))
					/ (my_date.c * my_date.rho * h / tau + a[n + 1] / h + a[n] / h);
					C[1] -= A[1]*kappa[1];
	            	C[n-1] -= B[n-1]*kappa[0];
	            	F[1] += A[1]*mu[1]; F[n-1] += B[n-1]*mu[0];
					y_tmp = progon(n, A,C, B, F, my_date);
					y_tmp[0] = kappa[1] * y_tmp[1] + mu[1];
					y_tmp[n] = kappa[0]*  y_tmp[n-1] + mu[0];
				}
			}
			else if (TEST_P == 2)
			{
				F[1] += A[1]* my_date.left_boarder(0, my_date);
		   		F[n-1] += B[n-1] * my_date.right_boarder(my_date.L, my_date);
				y_tmp[0] = my_date.left_boarder(0, my_date);
           		y_tmp[n] = my_date.right_boarder(my_date.L, my_date);
				y_tmp = progon(n,A, C, B, F,my_date);
			}
		} while (norm(y1,y_tmp)> EPS);

		fout << j << '\t' <<  iter << "\n";
		y2 = y1;
	}

	fout.close();
	return y1;
}
