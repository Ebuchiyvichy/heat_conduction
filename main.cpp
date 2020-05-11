
//
// Created by irina on 21.03.2020.
//
#include "error.h"

int main()
{
	Date	DATA;
	double	h = 0.01;
	double	tau = 0.01;
	double	T = 1.0;
	/*
		TEST_P == 1 - right edge with flux, left with temperature
	  	TEST_P == 2 - temperature
      	TEST_P == 3 - flux on edges
	*/

	std::cout << "Variant 2" << std::endl;
	DATA.left_boarder = &u0;
	DATA.right_boarder = &P2;	//���� �������

	integro_interpolation(DATA.L / h, T, h, tau, 4, 1, DATA);
	quasilinear(DATA.L / h, T, h, tau, 2, DATA);
	non_linear_iter(DATA.L / h, T, h, tau, 2, DATA);

	system("pause");
	return (0);
}