#include<iostream>
#include<cmath>

import Math;
using namespace Math;

void funC()
{

	Complex c1;
	std::cout << c1 << std::endl;
	Complex c2(2);
	std::cout << c2 << std::endl;
	Complex c3(2, 2);
	std::cout << c3 << std::endl;
	c1 = Complex::FromExponentialForm(2, 2);
	std::cout << c1 << std::endl;
	c1 = Complex::FromAlgebraicForm(3, 2);
	std::cout << c1 << std::endl;
	std::cout << c1.Re() << " " << c1.Im() << " " << c1.Mod() << " " << c1.Arg() << std::endl;
	std::cout << (double)c1 << std::endl;
	std::cout << -c1 << std::endl;
	std::cout << ++c1 << std::endl;
	std::cout << c1++ << std::endl;
	std::cout << --c1 << std::endl;
	std::cout << c1-- << std::endl;
	std::cout << (c1 -= c2) << std::endl;
	std::cout << (c1 += c2) << std::endl;
	std::cout << (c1 *= c2) << std::endl;
	std::cout << (c1 /= c2) << std::endl;
	std::cout << (c1 - c2 + c3 * c1 / c2) << std::endl;
	std::cout << FindGreatestCommonDivisor(8, 12) << std::endl;
	std::cout << FindLeastCommonMultiple(12, 8) << std::endl;
}

void funR()
{
	Rational c1;
	std::cout << c1 << std::endl;
	Rational c2(5);
	std::cout << c2 << std::endl;
	Rational c3(6, 4);
	std::cout << c3 << std::endl;
	std::cout << c3.Nominator() << " " << c3.Denominator() << std::endl;
	std::cout << (double)c3 << std::endl;
	std::cout << -c1 << std::endl;
	std::cout << ++c1 << std::endl;
	std::cout << c1++ << std::endl;
	std::cout << --c1 << std::endl;
	std::cout << c1-- << std::endl;
	std::cout << (c1 -= c2) << std::endl;
	std::cout << (c1 += c2) << std::endl;
	std::cout << (c1 *= c2) << std::endl;
	std::cout << (c1 /= c2) << std::endl;
	std::cout << (c1 - c2 + c3 * c1 / c2) << std::endl;
	std::cout << (c1 == c2) << std::endl;
	std::cout << c1 << " " << c2 << std::endl;

	std::cout << (c1 < c2) << std::endl;
	std::cout << (c1 > c2) << std::endl;
	std::cout << (c1 <= c2) << std::endl;
	std::cout << (c1 >= c2) << std::endl;


}

void main()
{
	funC();
	funR();
	double u = 3;

	Complex f1 = f((Complex)u);
	Rational f2 = f((Rational)u);
	double f3 = f(u);
	std::cout << f1 << std::endl;
	std::cout << f2 << std::endl;
	std::cout << f3 << std::endl;
}