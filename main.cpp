
//
// Created by irina on 21.03.2020.
//
#include "error.h"

int main()
{
	Date	DATA;
	double	h = 0.01;
	double	tau = 0.01;
	double	T = 1;

	std::cout << "Variant 2 (a)" << std::endl;
	DATA.left_boarder = &u0;
	DATA.right_boarder = &u0;	//���� �������

	integro_interpolation(DATA.L / h, T, h, tau, 2, 1, DATA);
//	integro_interpolation(DATA.L / h, T, h, tau, 2, 0.5, DATA);
//	integro_interpolation(DATA.L / h, 0.3, h, tau, 2, 0.2, DATA);
// 	quasilinear(DATA.L / h, T, h, tau, 2, DATA);
//	non_linear(DATA.L / h, T, h, tau, 2, DATA);
	/*TEST_P == 1 - right edge with flux, left with temporary temperature
	  TEST_P == 2 - temporary temperature
      TEST_P == 3 - flux on edges*/
//	error_check(h, tau, 0, 1, DATA);
	std::cout << "End";
	system("pause");
	return (0);
}