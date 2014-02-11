#include <cassert>
#include <cmath>
#include <complex>
#include <vector>
#include <Eigen/Core>

#include "Const.h"
#include "Tract.h"

using namespace std;
using namespace Eigen;


Tract::Tract(double c1, double dl, const double S[], int length) :
	a	(130 * Const::pi),
	b	(std::pow(30 * Const::pi, 2)),
	om0	(pow(406 * Const::pi, 2)),
	c1	(c1),
	dl	(dl),
	S	(length)
{
	for (int i = 0; i < length; i++)
		this->S[i] = S[i];
}
Tract::Tract(double c1, double dl, vector<double> &S) :
	a	(130 * Const::pi),
	b	(std::pow(30 * Const::pi, 2)),
	om0	(pow(406 * Const::pi, 2)),
	c1	(c1),
	dl	(dl),
	S	(S)
{
}

Tract::~Tract(void)
{
}

complex<double> Tract::GetBeta(double f)
{
	double omg = 2.0 * Const::pi * f;
	complex<double> s(0, omg);
	return s * om0 / ((s + a) * s + b) + sqrt(s * c1);
}

vector<double> Tract::GetArea(void)
{
	return S;
}

double Tract::GetStartArea(void)
{
	return S.front();
}

double Tract::GetEndArea(void)
{
	return S.back();
}

double Tract::GetElemLength(void)
{
	return dl;
}

Matrix2cd Tract::ChainMatrix(double freq, int from)
{
	return ChainMatrix(freq, from, S.size() - 1);
}
Matrix2cd Tract::ChainMatrix(double freq, int from, int to)
{
	double n = to - from + 1;

	double omg = 2.0 * Const::pi * freq;
	complex<double> s(0, omg);

	complex<double> loss = s * om0 / ((s + a) * s + b) + sqrt(s * c1);
	complex<double> gamma = sqrt((a + s) / (loss + s));
	complex<double> sigma = gamma * (loss + s) / Const::c;

	Matrix2cd r = Matrix2cd::Identity();
	complex<double> sh = sinh(sigma * dl);
	Matrix2cd K;
    K(0, 0) = K(1, 1) = cosh(sigma * dl);

    for (int i = 0; i < n; i++)
    {
        complex<double> coeff = -Const::rho * Const::c / S[from + i] * gamma;

		K(0, 1) = coeff * sh;
        K(1, 0) = 1.0 / coeff * sh;

        r = K * r;
    }

    return r;
}
