//
// Created by irina on 21.03.2020.
//
#include "function.h"

int main()
{
	Date	DATA;

	std::cout << "Variant 2 (a)" << std::endl;
	DATA.left_boarder = &u0;
	DATA.right_boarder = &P2;
	integro_interpolation(15, 10, 0.1, 0.1, 1, DATA);

	DATA.test = 'b';
	std::cout << "Variant 2 (b)" << std::endl;
	DATA.left_boarder = &u0;
	DATA.right_boarder = &P2;
	integro_interpolation(15, 10, 0.1, 0.1, 1, DATA);
	system("pause");
	return (0);
}
