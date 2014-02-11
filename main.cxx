#include <iostream>
#include <cstdlib>
#include <vector>
#include <Eigen/Core>

#include "Tract.h"
#include "SpeechSynthesizer.h"

#define VOWEL_A 0.45, 0.20, 0.26, 0.21, 0.32, 0.30, 0.33, 1.05, 1.12, 0.85, 0.63, 0.39, 0.26, 0.28, 0.23, 0.32, 0.29, 0.28, 0.40, 0.66, 1.20, 1.05, 1.62, 2.09, 2.56, 2.78, 2.86, 3.02, 3.75, 4.60, 5.09, 6.02, 6.55, 6.29, 6.27, 5.94, 5.28, 4.70, 3.87, 4.13, 4.25, 4.27, 4.69, 5.03


static void TestSynth(int freq, int resol, double sec, double h_min)
{
#ifdef _MSC_VER
	double s[] = { VOWEL_A };
	int n = sizeof(s) / sizeof(s[0]);
	std::vector<double> a(n);
	for (int i = 0; i < n; i++)
		a[i] = s[i];

#else
	std::vector<double> a{ VOWEL_A };
#endif

	SpeechSynthesizer synthesizer(freq, resol, a, 0.396825);
	synthesizer.h_min = h_min;

	std::printf("#t       Ag1       v1        Ag2       v2        ug        p\n");

	for (int i = 0; i < sec * freq; i++)
	{
		synthesizer.Step();
	}
}

int main(int argc, char *argv[])
{
	int freq, resol;
	double sec;
	double h_min;

	if (argc - 1 < 4)
	{
		freq = 20000;
		resol = 256;
		sec = 1.0;
		h_min = 2.0E-3;
	}
	else
	{
		freq = std::atoi(argv[1]);
		resol = std::atoi(argv[2]);
		sec = std::atof(argv[3]);
		h_min = std::atof(argv[4]);
	}

	std::cout << "# sampling   : " << freq << std::endl;
	std::cout << "# resolution : " << resol << std::endl;
	std::cout << "# time       : " << sec << std::endl;
	std::cout << "# h_min      : " << h_min << std::endl;

	TestSynth(freq, resol, sec, h_min);

	return 0;
}
