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
void impulseResponseLPF(double Fs, double Fc, unsigned short int num_taps, std::vector<double> &h, double decim)
{

	double cutoff = Fc/((Fs/decim)/2);
	// allocate memory for the impulse response
	h.resize(num_taps, 0.0);
	for(auto i = 0 ; i < num_taps ; i++){
		if(i == (num_taps-1)/2){
			h[i] = cutoff;
		}else{
			h[i] = cutoff * sin( PI * cutoff *( i-(num_taps-1)/2 ) ) / ( PI * cutoff *( i-(num_taps-1)/2 ) );
		}
		h[i] = h[i] * (sin(i * PI / num_taps)*sin(i * PI / num_taps));
		//printf("h[%d] = %f\n",i,h[i]);
	}

}


void impulseResponseBPF(double Fs, double Fb,double Fe, unsigned short int num_taps, std::vector<double> &h, double decim)
{

	double center = ((Fe+Fb)/2.0)/((Fs/decim)/2.0);
    double pass = (Fe-Fb)/((Fs/decim)/2.0);
	// allocate memory for the impulse response
	h.resize(num_taps, 0.0);
	for(auto i = 0 ; i < num_taps ; i++){
		if(i == (num_taps-1)/2){
			h[i] = pass;
		}else{
			h[i] = pass * sin( PI * (pass/2) *( i-(num_taps-1)/2 ) ) / ( PI * (pass/2) *( i-(num_taps-1)/2 ) );
		}
        h[i] = h[i] * cos(i*PI*center);
		h[i] = h[i] * (sin(i * PI / num_taps)*sin(i * PI / num_taps));
		//printf("h[%d] = %f\n",i,h[i]);
	}

}


// function to compute the filtered output "y" by doing the convolution
// of the input data "x" with the impulse response "h"
/*void convolveFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h)
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
}*/


void convolveFIR_N_dec(const int step_size, std::vector<double> &y, const std::vector<double> &x, const std::vector<double> &h, std::vector<double> &state )
{
	auto max_size = x.size();
	if (h.size() > max_size)
	{
		max_size = h.size();
	}
	long special = 0;
	for (auto n = 0; n < y.size(); n++)
	{
		special = 0;
		y[n] = 0;
		for (auto m = 0; m < h.size(); m++)
		{
			if ((step_size*n-m) >= 0 && (step_size*n-m) < max_size)
			{
				y[n] += x[step_size*n-m] * h[m];
			}else if((step_size*n-m) < 0 && state.size() > 0){
				y[n] += state[state.size() - 1 - special] * h[m];
				special++;
			}


		}
	}
	for(auto ii = 0 ; ii < state.size(); ii++){
		state[ii] = x[(x.size()) - state.size() + ii];
	}


}

void convolve_UPSAMPLE_N_dec(const int step_size, const int upsample_size, std::vector<double> &y, const std::vector<double> &x, const std::vector<double> &h, std::vector<double> &state)
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

	// for(auto ii = 0 ; ii < y.size(); ii++){
	// 	y[ii] =y[ii]*upsample_size;
	// }

//	for(auto i = 0 ; i < y.size(); i++){y[i] = y[i]*upsample_size;}


	for(auto ii = 0 ; ii < state.size(); ii++){
		state[ii] = x[(x.size()) - state.size() + ii];
	}

}

void fmDemodArctanBlock(std::vector<double> &fm_demod,std::vector<double> &I, std::vector<double> &Q,std::vector<double> &prev_phase){
	fm_demod.resize(I.size(), 0.0);
	double thetadelta = 0, a, b, c, current_phase;
	for(auto n = 0; n < I.size(); n++){
		//std::cout << "Bad Samples -> I :" << I[0] << " Q : " << Q[0] << " \n";
		a = b =c = current_phase = 0;
		if(n == 0){
			a = I[n]*(Q[n]-prev_phase[0]);		//prev phase is never being stored
			b = Q[n]*(I[n]-prev_phase[1]);
			c = (I[n]*I[n] + Q[n]*Q[n]);
			thetadelta = (a-b)/c;
		}
		else{
			a = I[n]*(Q[n]-Q[n-1]);
			b = Q[n]*(I[n]-I[n-1]);
			c = (I[n]*I[n] + Q[n]*Q[n]);

			thetadelta = (a-b)/c;
		}
        //thetadelta = fmod(thetadelta,PI);
        // int saftey = 0;
        // while(thetadelta > PI || thetadelta < -PI){
        //     if(saftey == 3){
        //         break;
        //     }
        //     if(thetadelta > PI){
        //         thetadelta = thetadelta - 2*PI;
        //     }
        //     else if(thetadelta < -PI){
        //         thetadelta = thetadelta + 2*PI;
        //     }
        //     saftey++;
        // }                                               //Unwrap

		if(!std::isnan(thetadelta)){
		fm_demod[n] = thetadelta;
		}else{
		fm_demod[n] = (a-b)*2;
		}
	}
	prev_phase.resize(2);
	prev_phase[0] = Q[Q.size() - 1];
	prev_phase[1] = I[I.size() - 1];
}
