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
void resampler(double Fs, double Fc, unsigned short int num_taps, std::vector<double> &h, double decim)
{

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
				impulseResponseRRC[k] = (sin(PI*t*(1-beta)/T_symbol) +  4*beta*(t/T_symbol)*cos(PI*t*(1+beta)/T_symbol))/(PI*t*(1-(4*beta*t/T_symbol)*(4*beta*t/T_symbol))/T_symbol)

			}
	}
}
