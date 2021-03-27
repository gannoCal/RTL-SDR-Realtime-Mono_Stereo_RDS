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
	for(auto phase = 0; phase < upsample_size; phase++){

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
						}else if(((step_size%upsample_size)*n-m) < 0 && state.size() > 0){
							y[n] += state[state.size() - 1 - special] * h[m];
							special++;
						}

					}

				// sleep(1);
				// std::cout << " n: " << n <<" y[n] : " <<  y[n] << " \n";
			}

	}
}

void Manchester_and_differntial(rds_data_RRC,logicdata)
{
	// since we want to get the middle sample lets offset
// 1 represnts a high , 0 represents a below
// MAnchester coding HL -> 1. LH -> 0
	int v = k/2
	for(auto i = 0; i < rds_data_RRC.size();i = i + 2*v){
		if(rds_data_RRC[i] > 0 && rds_data_RRC[i+v] < 0){prediff[counter] = 1;}
		else{prediff[counter] = 0;}
	}

	logicdata.resize(prediff.size());
	logicdata[0] = prediff[0];
	for(auto i = 1; i < logicdata.size();i = i + 1){
		logicdata[i] = prediff[i]^prediff[i-1];
	}
}

void impulseResponseRootRaisedCosine(Fs, N_taps,impulseResponseRRC)
{
	double T_s = 1/2375;
	float beta = 0.9;
	impulseResponseRRC.resize(N_taps);
	for(int k = 0; k < N_taps;k++){
			t = float((k-N_taps/2))/Fs;

			if(t == 0.0){
				impulseResponseRRC[k] = 1.0 + beta*((4/PI)-1);
			}
			else if(t == -T_symbol/(4*beta) || t == T_symbol/(4*beta)){
				impulseResponseRRC[k] = (beta/np.sqrt(2))*(((1+2/PI)*(sin(PI/(4*beta)))) + ((1-2/PI)*(cos(PI/(4*beta)))));
			}
			else{
				impulseResponseRRC[k] = (sin(PI*t*(1-beta)/T_symbol) +  4*beta*(t/T_symbol)*cos(PI*t*(1+beta)/T_symbol))/ (PI*t*(1-(4*beta*t/T_symbol)*(4*beta*t/T_symbol))/T_symbol)

			}
	}
}
