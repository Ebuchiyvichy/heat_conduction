//
// Created by irina on 17.02.2020.
//

#ifndef DIFF_SOLVE_MATRIXCLASS_H
#define DIFF_SOLVE_MATRIXCLASS_H

#endif //DIFF_SOLVE_MATRIXCLASS_H

#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <functional>
#include <string>
double			EPS = 10.e-3;
int	const		dim = 2;
double const    t0 = 0;
double const	PI = 3.14159265358979;

std::function<std::vector<double>(std::vector<double>, double tau)> func; // указатель на функцию для конкретного теста

class Matrix
{
public:
	std::vector<std::vector<double>> value;
	size_t size;

public:
	Matrix(size_t count)
	{
		size = count;
		for (int i = 0; i != size; i++)
		{
			value.push_back(std::vector<double>(count));
		}
	}
	~Matrix() {

	}

	void init()
	{
		std::ifstream file;
		file.open("matrixexample.txt");

		for (int i = 0; i != size; i++)
			for (int j = 0; j != size; j++)
				file >> value[i][j];
		file.close();
	}
	void trunc()
	{
		double temp;
		for (int i = 0; i != size; i++)
		{
			for (int j = i + 1; j != size; j++)
			{
				temp = value[i][j];
				value[i][j] = value[j][i];
				value[j][i] = temp;
			}
		}
	}
	void print()
	{
		for (int n = 0; n != size; n++)
		{
			std::cout << "\n";
			for (int j = 0; j != size; j++)
			{
				std::cout << std::setw(8) << value[n][j] << '\t';
			}
		}
		std::cout << "\n";
	}

	void onebyone()
	{
		for (int i = 0; i != size; i++)
			for (int j = 0; j != size; j++)
			{
				if (i == j)
					value[i][i] = 1;
				else
					value[i][j] = 0;
			}
	}
	friend void cpy(Matrix &A, Matrix& B)
	{
		for (int i = 0; i != A.size; i++)
			for (int j = 0; j != A.size; j++)
				B.value[i][j] = A.value[i][j];
	}

	friend Matrix operator + (const Matrix &A, const Matrix &B)
	{
		Matrix C(A.size);
		if (A.size != B.size)
			std::cout << "Not correct size A and B";
		else
		{
			for (int i = 0; i != A.size; i++) {
				for (int j = 0; j != A.size; j++)
					C.value[i][j] = A.value[i][j] + B.value[i][j];
			}
		}
		return C;
	}
	friend Matrix operator - (const Matrix &A, const Matrix &B)
	{
		Matrix C(A.size);
		if (A.size != B.size && C.size != B.size)
			std::cout << "Not correct size A and B";
		else
		{
			for (int i = 0; i != A.size; i++) {
				for (int j = 0; j != A.size; j++)
					C.value[i][j] = A.value[i][j] - B.value[i][j];
			}
		}
		return C;
	}
	friend Matrix operator *(const Matrix &A, const Matrix &B)
	{
		Matrix C(A.size);
		for (int i = 0; i != A.size; i++)
			for (int j = 0; j != B.size; j++)
			{
				double sum = 0.0;
				for (int k = 0; k != A.size; k++)
					sum += A.value[i][k] * B.value[k][j];
				C.value[i][j] = sum;
			}
		return C;
	}
	friend Matrix operator * (Matrix &A, double c)
	{
		for (int i = 0; i != A.size; i++) {
			for (int j = 0; j != A.size; j++)
				A.value[i][j] = A.value[i][j] * c;
		}
		return A;
	}
	friend std::vector<double> operator * (const Matrix &A, const std::vector<double> &b) {
		std::vector<double> c(b.size());
		if (A.size != b.size())
			std::cout << "Wrong size in vector or matrix!\n";
		else {
			for (int i = 0; i != A.size; i++)
				for (int j = 0; j != A.size; j++) {
					c[i] += A.value[i][j] * b[j];
				}
		}
		return c;
	}
	friend std::vector<double> operator * (const std::vector<double> &b, const Matrix &A) {
		std::vector<double> c(b.size());
		if (A.size != b.size())
			std::cout << "Wrong size in vector or matrix!\n";
		else {
			for (int i = 0; i != A.size; i++)
				for (int j = 0; j != A.size; j++) {
					c[i] += A.value[j][i] * b[j];
				}
		}
		return c;
	}
	//     Matrix& operator = (const Matrix &C)
	//    {
	//        for (int i = 0; i != size; i++)
	//            for (int j = 0; j != size; j++)
	//                value[i][j] = C.value[i][j];
	//    }
//	void	QR(size_t n, Matrix &A);
	void QR_find_x(Matrix& A)
	{
		double	c;
		double	s;
		double	tmp;

		for (int k = 0; k < size; k++)
		{
			for (int i = k + 1; i < size; i++)
			{
				if (fabs(A.value[k][i]) > 10e-8)
				{
					c = A.value[k][k] / sqrt(A.value[k][k] * A.value[k][k] + A.value[i][k] * A.value[i][k]);
					s = A.value[i][k] / sqrt(A.value[k][k] * A.value[k][k] + A.value[i][k] * A.value[i][k]);
					for (int j = 0; j <= i; j++) // change T-matrix
					{
						tmp = value[k][j];
						value[k][j] = value[k][j] * c + value[i][j] * s;
						value[i][j] = c * value[i][j] - s * tmp;
					}
					for (int j = k; j < size; j++) // change A-matrix
					{
						tmp = A.value[k][j];
						A.value[k][j] = c * A.value[k][j] + s * A.value[i][j];
						A.value[i][j] = c * A.value[i][j] - s * tmp;
					}
				}
			}
		}
	}
	void inverse_matrix(const Matrix& R, Matrix& T)
	{
		//with Hausse
		std::vector<double> x(T.size);
		std::vector<double> y(T.size);
		Matrix E(R.size);

		E.onebyone();
		for (int i = 0; i < T.size; i++)
		{
			for (int j = 0; j < T.size; j++)
			{
				if (j == i)
					y[j] = 1;
				else
					y[j] = 0;
			}
			x = find_x(R, multi_vect(y, T));
			for (int j = 0; j < T.size; j++)
			{
				value[i][j] = x[j];
			}
		}
		trunc();
	}
	friend std::vector<double>	find_x(const Matrix& A, std::vector<double> b_)
	{
		std::vector<double>	x(A.size);

		for (int i = A.size - 1; i >= 0; i--)
		{
			for (int j = i + 1; j < A.size; j++)
				b_[i] = b_[i] - A.value[i][j] * b_[j];
			b_[i] = b_[i] / A.value[i][i];
		}
		x = cpy_vector(x, b_, A.size);
		return (x);
	}
	friend std::vector<double>	multi_vect(std::vector<double> I, const Matrix& T)
	{
		std::vector<double> tmp(dim);

		for (int j = 0; j < T.size; j++)
			for (int i = 0; i < T.size; i++)
				tmp[j] += T.value[j][i] * I[i];
		return (tmp);
	}
};

