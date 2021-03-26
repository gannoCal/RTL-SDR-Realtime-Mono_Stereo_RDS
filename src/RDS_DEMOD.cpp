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
int Fc = 3000;
int k = 24;


	impulseResponseLPF( Fs, Fc, num_taps,RDS_coeff, 1);
if(mode == 0){resampler(80,19, rds_data_ds, rds_data, rds_coeff,state_stereo_data);}
if(mode == 1){resampler(80,19, rds_data_ds, rds_data, rds_coeff,state_stereo_data);}

	// low pass + resampler - end

// RRC begin
impulseResponseRootRaisedCosine(Fs, N_taps,impulseResponseRRC)
convolveFIR_N_dec(1, rds_data_RRC, rds_data_ds, impulseResponseRRC, std::vector<double> &state )

// RRC end

// CDR begin

// CDR end

	return 0;
}
