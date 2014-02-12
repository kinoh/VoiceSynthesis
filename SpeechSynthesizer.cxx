#include <cmath>
#include <complex>
#include <cstring>
#include <algorithm>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv.h>
#include <Eigen/Core>

#include "Const.h"
#include "SpeechSynthesizer.h"

using namespace std;
using namespace Eigen;


static const double NasalCavityArea[] = { 1, 2, 3, 4, 6, 8, 8, 7, 4, 2, 2 };

SpeechSynthesizer::SpeechSynthesizer(int sampling, int freqResolution, vector<double> &tractArea, double dl) :
	m1	(0.125),
	m2	(0.025),
	d1	(0.25),
	d2	(0.05),
	lg	(1.4),
	es	(100),
	eh	(500),
	ks1	(8.0E+4),
	ks2	(8.0E+3),
	kh1	(3 * ks1),
	kh2	(3 * ks2),
	kc	(2.5E+4),
	om1	(sqrt(ks1 / m1)),
	Ag0	(0.05),
	p_s	(7.85E+3),
	x1	(0),
	x2	(0),
	v1	(0),
	v2	(0),
	ug_buf	(freqResolution),
	toVelum		(8.0),
    toSinus		(7.0),
    toConstr	(3.0),
    R_sin		(1.0),
    L_sin		(5.94E-3),
    C_sin		(15.8E-6),
	oral		(4.0, dl, tractArea),
	nasal		(72.0, 1.0, NasalCavityArea, 11),
	r_N(0.2),
	r_vib(3.0),
	dt	(1.0 / sampling),
	df	(1.0 / (freqResolution * dt)),
	N	(freqResolution),
	calcSpan(0.020),
	t(0)
{
	z_in_curr = new double[N];
	h_out_curr = new double[N];
	z_in_prev = new double[N];
	h_out_prev = new double[N];

	for (int i = 0; i < N; i++)
		ug_buf.push_front(0);

	lastCalcTime = -2 * calcSpan;
	CalcTract();
	memcpy(z_in_prev, z_in_curr, N);
	memcpy(h_out_prev, h_out_curr, N);
	lastCalcSpan = calcSpan;

	system.function = Differentiate;
	system.jacobian = nullptr;
	system.dimension = 5;
	system.params = this;
	step = gsl_odeiv_step_alloc(gsl_odeiv_step_rk4, system.dimension);
}

SpeechSynthesizer::~SpeechSynthesizer(void)
{
	delete [] z_in_curr;
	delete [] h_out_curr;
	delete [] z_in_prev;
	delete [] h_out_prev;
}

int SpeechSynthesizer::Differentiate(double t, const double y[], double dydt[], void *params)
{
	SpeechSynthesizer *s = (SpeechSynthesizer*) params;
	const double h_min = s->h_min;
	const double x1 = y[0];
	const double v1 = y[1];
	const double x2 = y[2];
	const double v2 = y[3];
	const double ug = y[4];
	double *dx1 = &dydt[0];
	double *dv1 = &dydt[1];
	double *dx2 = &dydt[2];
	double *dv2 = &dydt[3];
	double *dug = &dydt[4];

	const int n = s->N;
	double *z_in = new double[n];
	double V = 0;

	s->GetTractImpedance(z_in, t);
	for (int i = 0; i < n; i++)
	{
		V += z_in[i] * s->ug_buf[i];
	}

	double xc = -s->Ag0 / (2 * s->lg);
	double h1 = x1 - xc;
	double h2 = x2 - xc;
	if (h1 > 0 && h1 < h_min) h1 = h_min;
	if (h2 > 0 && h2 < h_min) h2 = h_min;
	double Rv = 1.5 * Const::mu / s->lg;
	double Rv1 = Rv * s->d1 / gsl_pow_3(h1);
	double Rv2 = Rv * s->d2 / gsl_pow_3(h2);
	double Rbb = Const::rho / (8 * gsl_pow_2(s->lg));
	double Rc = Rbb * 1.37 / gsl_pow_2(h1) * (ug);
	double R12 = Rbb * (1 / gsl_pow_2(h2) - 1 / gsl_pow_2(h1)) * (ug);
	double Re = -Rbb * 0.5 / gsl_pow_2(h2) * (ug);
	double Lg = Const::rho / (2 * s->lg);
	double Lg1 = Lg * s->d1 / h1;
	double Lg2 = Lg * s->d2 / h2;

	*dug = (h1 <= 0 || h2 <= 0
			? -ug / (t - s->t)
			: (s->p_s - V - (Rc + Rv1 + R12 + Rv2 + Re) * ug) / (Lg1 + Lg2));

	double P1, P2;
	if (h1 <= 0)
	{
		P1 = s->p_s;
		P2 = 0;
	}
	else if (h2 <= 0)
	{
		P1 = s->p_s;
		P2 = P1;
	}
	else
	{
		P1 = s->p_s - (Rc + Rv1 / 2) * ug - Lg1 / 2 * *dug;
		P2 = P1 - (R12 + (Rv1 + Rv2) / 2) * ug - (Lg1 + Lg2) / 2 * *dug;
	}

	double s1 = s->ks1 * x1 * (1 + s->es * gsl_pow_2(x1));
	if (h1 <= 0) s1 += s->kh1 * h1 * (1 + s->eh * gsl_pow_2(h1));
	double s2 = s->ks2 * x2 * (1 + s->es * gsl_pow_2(x2));
	if (h2 <= 0) s2 += s->kh2 * h2 * (1 + s->eh * gsl_pow_2(h2));
	double r1 = 2 * (h1 > 0 ? 0.2 : 1.1) * sqrt(s->ks1 * s->m1);
	double r2 = 2 * (h2 > 0 ? 0.6 : 1.9) * sqrt(s->ks2 * s->m2);
	
	*dv1 = (-r1 * v1 - s1 - s->kc * (x1 - x2) + P1 * s->lg * s->d1) / s->m1;
	*dv2 = (-r2 * v2 - s2 - s->kc * (x2 - x1) + P2 * s->lg * s->d2) / s->m2;

	*dx1 = v1;
	*dx2 = v2;

	delete [] z_in;

	return GSL_SUCCESS;
}

double SpeechSynthesizer::Step(void)
{
	CalcTract();

	double y[5] = { x1, v1, x2, v2, ug_buf.front() };
	double yerr[5];

	gsl_odeiv_step_apply(step, t, dt, y, yerr, t == 0 ? nullptr : dydt, dydt, &system);
	t += dt;

	if (yerr[0] > 1E-3)
		std::cerr << "large error" << std::endl;

	x1 = y[0];
	v1 = y[1];
	x2 = y[2];
	v2 = y[3];
	ug_buf.pop_back();
	ug_buf.push_front(y[4]);

	double *h_out = new double[N];
	double p = 0;

	GetOutputImpedance(h_out, t);
	for (int i = 0; i < N; i++)
		p += h_out[i] * ug_buf[i];

	delete [] h_out;

	std::printf("%f %.3e %.3e %.3e %.3e %.3e %.3e\n", t, x1 * 2 * lg + Ag0, v1, x2 * 2 * lg + Ag0, v2, ug_buf.front(), p);

	return p;
}

double SpeechSynthesizer::GetTime(void)
{
	return t;
}

void SpeechSynthesizer::CalcTract(void)
{
	double dur = t - lastCalcTime;

	if (dur < calcSpan)
		return;

	lastCalcSpan = dur;
	lastCalcTime = t;

#ifdef _MSC_VER
	complex<double> *Z_in = new complex<double>[N - 1];
	complex<double> *H_out = new complex<double>[N - 1];
#else
	complex<double> Z_in[N - 1];
	complex<double> H_out[N - 1];
#endif

	int i_v = (int) (toVelum / oral.GetElemLength());
	int i_c = i_v + (int) (toConstr / oral.GetElemLength());
	int i_s = (int) (toSinus / nasal.GetElemLength());

    for (int i = 1; i < N; i++)
    {
		double f = i * df;
        double omg = 2 * Const::pi * f;
        complex<double> s(0, omg);
		
		Matrix2cd K_G = oral.ChainMatrix(f, 0, i_v);
        Matrix2cd K_C = oral.ChainMatrix(f, i_v, i_c);
        Matrix2cd K_L = oral.ChainMatrix(f, i_c);
		double r_L = sqrt(oral.GetEndArea() / Const::pi);
		double ka = omg / Const::c * r_L;
		complex<double> Z_L = Const::rho * Const::c / oral.GetEndArea() * complex<double>(ka * ka / 2, 8.0 * ka / (3.0 * Const::pi));

        complex<double> Z_sin = R_sin + L_sin * s + 1.0 / (C_sin * s);
        Matrix2cd K_N1 = nasal.ChainMatrix(f, 0, i_s);
        Matrix2cd K_N2 = nasal.ChainMatrix(f, i_s);
        Matrix2cd K_sin;
		K_sin << 1, 0, -1.0 / Z_sin, 1;
        Matrix2cd K_N = K_N2 * K_sin * K_N1;
        complex<double> Z_N = 4 * Const::pi * gsl_pow_2(r_N) * Const::rho * Const::c / (gsl_pow_2(Const::c) + gsl_pow_2(r_N * omg)) * (gsl_pow_2(r_N * omg) + Const::c * r_N * s);

        Matrix2cd K_oral = K_L * K_C;
		complex<double> Z_VN = (K_N(1, 1) * Z_N - K_N(0, 1)) / (K_N(0, 0) - K_N(1, 0) * Z_N);

		complex<double> Z_VT = (K_oral(1, 1) * Z_L - K_oral(0, 1)) / (K_oral(0, 0) - K_oral(1, 0) * Z_L);
        Matrix2cd K_cN, K_cT;
		K_cN << 1, 0, -1.0 / Z_VN, 1;
        K_cT << 1, 0, -1.0 / Z_VT, 1;

        Matrix2cd K_fric = K_C * K_cN * K_G;
        Matrix2cd K_tract = K_L * K_fric;
        Matrix2cd K_nasal = K_N * K_cT * K_G;

		//double Z_0 = Const::rho * Const::c / oral.GetEndArea();
        Z_in[i - 1] = (K_tract(1, 1) * Z_L - K_tract(0, 1)) / (K_tract(0, 0) - K_tract(1, 0) * Z_L);

		complex<double> H_vib = oral.GetStartArea() / Const::c * s * r_vib / (Const::c + s * r_vib) * Z_in[i - 1] * oral.GetBeta(f);

        H_out[i - 1] = Z_L / (K_tract(0, 0) - K_tract(1, 0) * Z_L)
            + Z_N / (K_nasal(0, 0) - K_nasal(1, 0) * Z_N)
            + H_vib;
    }

	swap(z_in_curr, z_in_prev);
	swap(h_out_curr, h_out_prev);
	CalcResponse(z_in_curr, Z_in);
	CalcResponse(h_out_curr, H_out);

#ifdef _MSC_VER
	delete [] Z_in;
	delete [] H_out;
#endif
}

inline static complex<double> Exp(const complex<double> z)
{
    double re = real(z);
    double im = imag(z);
    return exp(re) * complex<double>(cos(im), sin(im));
}

void SpeechSynthesizer::CalcResponse(double res[], const complex<double> spec[])
{
    for (int i = 0; i < N; i++)
    {
		double t = i * dt;
		double c_t = 0.54 - 0.46 * cos(Const::pi * (1 + (double)i / N));
		complex<double> r(0, 0);

		for (int j = 1; j < N; j++)
		{
	        double omg = 2 * Const::pi * j * df;
			complex<double> s(0, omg);
			double c_f = 1 / (1 + exp(4 * (8.0 * j / N - 5)));	// designed to decrease around 5/8

			r += c_t * c_f / (2 * Const::pi * N) * spec[j - 1] * Exp(s * t);
		}

		res[i] = real(r);
	}
}

void SpeechSynthesizer::GetTractImpedance(double res[], double t)
{
	double dur = t - lastCalcTime;

	for (int i = 0; i < N; i++)
	{
		res[i] = z_in_curr[i] + (z_in_curr[i] - z_in_prev[i]) / calcSpan * dur;
	}
}
void SpeechSynthesizer::GetOutputImpedance(double res[], double t)
{
	double dur = t - lastCalcTime;

	for (int i = 0; i < N; i++)
	{
		res[i] = h_out_curr[i] + (h_out_curr[i] - h_out_prev[i]) / calcSpan * dur;
	}
}
