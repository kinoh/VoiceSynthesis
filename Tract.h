#pragma once
#include <cmath>
#include <complex>
#include <vector>
#include <Eigen/Core>


class Tract
{
private:
	const double a;		// ratio of wall resistance to mass
	const double b;		// squared angular frequency of mechanical resonance
	const double om0;	// lowest squared angular frequency of acoustic resonance
	const double c1;	// correction for thermal conductivity and viscosity
	const double dl;	// length of elementary section
    std::vector<double> S;	// area functions

public:
	Tract(double c1, double dl, const double S[], int length);
	Tract(double c1, double dl, std::vector<double> &S);
	~Tract(void);
	Eigen::Matrix2cd ChainMatrix(double freq, int from);
	Eigen::Matrix2cd ChainMatrix(double freq, int from, int to);
	std::vector<double> GetArea(void);
	double GetStartArea(void);
	double GetEndArea(void);
	double GetElemLength(void);
	std::complex<double> GetBeta(double f);
};
