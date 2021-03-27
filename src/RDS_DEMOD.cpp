/*
Comp Eng 3DY4 (Computer Systems Integration Project)
Copyright by Nicola Nicolici
Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"


int main(int argc,char** argv)
{

	std::vector<int> rds_ds;
	std::vector<int> rds_data;


/////// MIXER SHIT - ENDp
std::vector<double> tone_ds((block_size)/20, 0); // get from before
double PLL_freq = 3e3;
double Fs = 3e3;
if(mode == 1){Fs = 250e3;}
else{Fs = 250e3;}
double PLLNCOscale = 2.0;
double phaseAdjust = 0.0;
double normBandwidth = 0.01;
std::vector<double> ncoOut((block_size)/20, 0);
std::vector<double> prevstate(6,0);

	fmPll(tone_ds,PLLfreq,PLLfs,PLLNCOscale,phaseAdjust,normBandwidth,ncoOut,prevstate);
	for(auto i = 0; i < rds_ds.size(),i++){
		RDS_data[i] = ncoOut[i] * rds_ds[i] * 2;
	}
	/////// MIXER SHIT - END


// low pass + resampler - begin
int Fc = 3000;
int k = 24;


	impulseResponseLPF( Fs, Fc, num_taps,RDS_coeff, 1);
resampler(80,19, rds_data_ds, rds_data, rds_coeff,state_stereo_data);

	// low pass + resampler - end

// RRC begin
impulseResponseRootRaisedCosine(Fs, N_taps,impulseResponseRRC)
convolveFIR_N_dec(1, rds_data_RRC, rds_data_ds, impulseResponseRRC, std::vector<double> &state )

// RRC end

// CDR begin
	Manchester(rds_data_RRC,result)
// CDR end

	return 0;
}
