
//
// Created by irina on 21.03.2020.
//
#include "quasilinear.h"

int main()
{
	Date	DATA;
	double	h = 0.03;
	double	tau = 0.1;
	double	T = 10;

	std::cout << "Variant 2 (a)" << std::endl;
	DATA.left_boarder = &u0;
	DATA.right_boarder = &P2;
//	integro_interpolation(DATA.L / h, T, h, tau, 2, DATA);
//	quasilinear(DATA.L / h, T, h, tau, 2, DATA);
	non_linear(DATA.L / h, T, h, tau, 2, DATA);
	/*TEST_P == 1 - смешанная задача
	  TEST_P == 2 - постоянные значения на концах
      TEST_P == 3 - два потока на концах*/
	system("pause");
	return (0);
}