/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"
#define PI 3.14159265358979323846
// function to compute the impulse response "h" based on the sinc function
void impulseResponseLPF(float Fs, float Fc, unsigned short int num_taps, std::vector<float> &h)
{
	// allocate memory for the impulse response
	h.resize(num_taps, 0.0);
	auto norm_cutoff = Fc / (Fs / 2);

	// the rest of the code in this function is to be completed by you
	// based on your understanding and the Python code from the first lab
	for (auto i = 0; i < num_taps; i++)
	{
		if (i == (num_taps - 1) / 2)
		{
			h[i] = norm_cutoff;
		}
		else
		{
			h[i] = h[i] = norm_cutoff * ((sin(PI * norm_cutoff * (i - (num_taps - 1)/2))) / (PI * norm_cutoff * (i-(num_taps - 1)/2)));
		}
		h[i] = h[i] * pow(sin((i * PI) / (num_taps)), 2);
	}
}


// function to compute the filtered output "y" by doing the convolution
// of the input data "x" with the impulse response "h"
void convolveFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h)
{
	// allocate memory for the output (filtered) data
	y.resize(x.size()+h.size()-1, 0.0);

	// the rest of the code in this function is to be completed by you
	// based on your understanding and the Python code from the first lab
	auto max_size = x.size();
	if (h.size() > max_size)
	{
		max_size = h.size();
	}
	for (auto n = 0; n < y.size(); n++)
	{
		for (auto m = 0; m < h.size(); m++)
		{
			if ((n-m) >= 0 || (n-m) < max_size)
			{
				y[n] += x[n-m] * h[m];
			}
		}
	}
}
