/*
Comp Eng 3DY4 (Computer Systems Integration Project)
Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"
#include <unistd.h>
#include <math.h>
#include <valarray>
#define PI 3.14159265358979323846
// function to compute the impulse response "h" based on the sinc function
void resampler(const int step_size, const int upsample_size, std::vector<double> &y, const std::vector<double> &x, const std::vector<double> &h, std::vector<double> &state)
{
	auto max_size = x.size();
	if (h.size() > max_size)
	{
		max_size = h.size();
	}
	long special = 0;
	for(auto phase = 0; phase < upsample_size; phase++)
	{

		for (auto n = phase; n < y.size(); n = n + upsample_size)
		{
			int x_count = n;
			// std::cout << "phase : " <<  phase << " \n";
			special = 0;
			y[n] = 0;
			for (auto m = 0; m < h.size(); m = m + 1)
			{
				if (((step_size%upsample_size)*n-m) >= 0 && ((step_size%upsample_size)*n-m) < max_size)
				{
					y[n] += x[(step_size%upsample_size)*n-m] * h[m];
				} 
				else if(((step_size%upsample_size)*n-m) < 0 && state.size() > 0)
				{
					y[n] += state[state.size() - 1 - special] * h[m];
					special++;
				}
			}

			// sleep(1);
			// std::cout << " n: " << n <<" y[n] : " <<  y[n] << " \n";
		}

	}
}


void CDR(std::vector<double> rds_data_RRC, int k, int sample_point)
{
	std::vector<double> sample_point_array(10, 0);

	for(auto i = 0; i < 10; i++)
	{
		for(auto j = 1; j < k; j++)
		{
			if(abs(rds_data_RRC[i*k+j]) < abs(rds_data_RRC[i*k+j -1]) ){
				sample_point_array[i] = j;
			}

		}
	}

	// take the avarage of the array to find the avarage max value
	int total = 0;
	for(auto i = 0; i < sample_point_array.size();i++)
	{
		total = total + sample_point_array[i];
	}
	k = int(total/sample_point_array.size());
}

void Manchester_and_differntial(std::vector<double> rds_data_RRC, int sample_point, std::vector<double> logicdata, int k)
{
	// since we want to get the middle sample lets offset
	// 1 represnts a high , 0 represents a below
	// MAnchester coding HL -> 1. LH -> 0

	std::vector<double> prediff(rds_data_RRC.size(), 0.0);
	auto counter = 0;
	for(auto i = sample_point; i < rds_data_RRC.size();i = i + 2*sample_point)
	{
		if(rds_data_RRC[i] > 0 && rds_data_RRC[i+sample_point] < 0) {prediff[counter] = 1;}
		else {prediff[counter] = 0;}
		counter++;
	}

	logicdata.resize(prediff.size());
	logicdata[0] = prediff[0];
	for(auto i = 1; i < logicdata.size(); i++)
	{
		logicdata[i] = (!prediff[i] != !prediff[i-1]) ? 1 : 0; 
	}
}

void impulseResponseRootRaisedCosine(double Fs, int N_taps, std::vector<double> impulseResponseRRC)
{
	double T_s = 1.0/2375.0;
	float beta = 0.9;
	impulseResponseRRC.resize(N_taps);
	auto t = 0;
	for(int k = 0; k < N_taps;k++)
	{
		t = float((k-N_taps/2))/Fs;

		if(t == 0.0)
		{
			impulseResponseRRC[k] = 1.0 + beta*((4/PI)-1);
		}
		else if(t == -T_s/(4*beta) || t == T_s/(4*beta))
		{
			impulseResponseRRC[k] = (beta/sqrt(2))*(((1+2/PI)*(sin(PI/(4*beta)))) + ((1-2/PI)*(cos(PI/(4*beta)))));
		}
		else
		{
			impulseResponseRRC[k] = (sin(PI*t*(1-beta)/T_s) +  4*beta*(t/T_s)*cos(PI*t*(1+beta)/T_s))/ (PI*t*(1-(4*beta*t/T_s)*(4*beta*t/T_s))/T_s);
		}
	}
}
