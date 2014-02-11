#pragma once
#include <cmath>


namespace Const
{
#ifdef _MSC_VER
	static double pi = 4 * std::atan(1);
#else
	static constexpr double pi = 4 * std::atan(1.0);
#endif
	static const double c	= 3.5E+4;	// speed of sound in air
	static const double rho	= 1.14E-3;	// density of air
	static const double mu	= 1.86E-4;	// viscosity of air
}
