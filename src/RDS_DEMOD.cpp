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
std::vector<double> rds_carrier((block_size)/20, 0); // get from before
double PLL_freq = 3e3;
double Fs = 0;
if(mode == 1){Fs = 250e3;int k = 25;}
else{Fs = 240e3;int k = 24;}
double PLLNCOscale = 2.0;
double phaseAdjust = 0.0;
double normBandwidth = 0.01;
std::vector<double> ncoOut((block_size)/20, 0);
std::vector<double> prevstate(6,0);

	fmPll(rds_carrier,PLLfreq,PLLfs,PLLNCOscale,phaseAdjust,normBandwidth,ncoOut,prevstate);
	for(auto i = 0; i < rds_ds.size(),i++){
		rds_ds[i] = ncoOut[i] * rds_carrier[i] * 2;
	}
	/////// MIXER SHIT - END


// low pass + resampler - begin
int Fc = 3000;



	impulseResponseLPF( Fs, Fc, num_taps,RDS_coeff, 1);
	resampler(80,19, rds_data_ds, rds_ds, rds_coeff,state_rds_data);

	// low pass + resampler - end

// RRC begin
impulseResponseRootRaisedCosine(Fs, N_taps,impulseResponseRRC);
convolveFIR_N_dec(1, rds_data_RRC, rds_data_ds, impulseResponseRRC, std::vector<double> &state );

// RRC end

// CDR begin

// to find the best point to sample, take the first 10 symbols and find the largest points, then avarage them
CDR(rds_data_RRC,k, sample_point)


// CDR end

// DATA Processing begin
Manchester_and_differntial(rds_data_RRC,sample_point,logicdata);
// DATA Processing end
	return 0;
}
