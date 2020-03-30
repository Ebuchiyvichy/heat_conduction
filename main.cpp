
//
// Created by irina on 21.03.2020.
//
#include "quasilinear.h"

int main()
{
	Date	DATA;
	double	h = 0.03;
	double	tau = 0.01;
	double	T = 1;

	std::cout << "Variant 2 (a)" << std::endl;
	DATA.left_boarder = &u0;
	DATA.right_boarder = &P2;
//	integro_interpolation(DATA.L / h, T, h, tau, 2, 0, DATA);
	integro_interpolation(DATA.L / h, T, h, tau, 2, 0.5, DATA);
	integro_interpolation(DATA.L / h, T, h, tau, 2, 1, DATA);

	//	quasilinear(DATA.L / h, T, h, tau, 2, DATA);
//	non_linear(DATA.L / h, T, h, tau, 2, DATA);
	/*TEST_P == 1 - ��������� ������
	  TEST_P == 2 - ���������� �������� �� ������
      TEST_P == 3 - ��� ������ �� ������*/
	system("pause");
	return (0);
}