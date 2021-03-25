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
	std::vector<int> RDS_data;


/////// MIXER SHIT - ENDp
	fmPll(tone_ds,PLLfreq,PLLfs,PLLNCOscale,phaseAdjust,normBandwidth,ncoOut,prevstate);
	for(auto i = 0; i < rds_ds.size(),i++){
		RDS_data[i] = ncoOut[i] * rds_ds[i] * 2;
	}
	/////// MIXER SHIT - END


// low pass + resampler - begin
	impulseResponseLPF(double Fs, double Fc, unsigned short int num_taps, std::vector<double> &h, double decim);
	resampler(decimator, stereo_data_ds, stereo_data, audio_coeff,state_stereo_data);
	// low pass + resampler - end

// RRC begin
impulseResponseRootRaisedCosine(Fs, N_taps,impulseResponseRRC)
// RRC end

// CDR begin
impulseResponseRootRaisedCosine(Fs, N_taps,impulseResponseRRC)
// CDR end

	return 0;
}
