#pragma once
#include <vector>
#include <chrono>
#include <complex>
#include <boost/circular_buffer.hpp>
#include <gsl/gsl_odeiv.h>
#include "Tract.h"


class SpeechSynthesizer
{
public:
	double h_min;

private:
	// vocal chord
	const double m1;	// mass 1 of vocal cords
	const double m2;	// mass 2 of vocal cords
	const double d1;	// thickness of m1
	const double d2;	// thickness of m2
	const double lg;	// effective length of vocal chords
	const double es;	// nonlinear spring coeff.
    const double eh;
    const double ks1;	// linear spring coeff.
    const double ks2;
    const double kh1;
    const double kh2;
    const double kc;	// coupling spring coeff.
    const double om1;
	double Ag0;
	double p_s;
	double x1, x2, v1, v2;
	boost::circular_buffer<double> ug_buf;	// front : current (latest) value
	double dydt[5];
	// vocal tract
    const double toVelum;
    const double toSinus;
    const double toConstr;
    const double R_sin;
    const double L_sin;
    const double C_sin;
	Tract oral, nasal;
	double *z_in_curr, *h_out_curr;
	double *z_in_prev, *h_out_prev;
	// other params
	const double r_N;
	const double r_vib;
	// simulation params
	const double dt, df;
	const int N;
	const double calcSpan;
	double lastCalcSpan;
	double lastCalcTime;
	gsl_odeiv_system system;
	gsl_odeiv_step *step;
	double t;

public:
	SpeechSynthesizer(int sampling, int freqResolution, std::vector<double> &tractArea, double dl);
	~SpeechSynthesizer(void);
	double Step(void);
	double GetTime(void);

private:
	static int Differentiate(double t, const double y[], double dydt[], void *params);
	void CalcTract(void);
	void CalcResponse(double res[], const std::complex<double> spec[]);
	void GetTractImpedance(double res[], double t);
	void GetOutputImpedance(double res[], double t);
};
