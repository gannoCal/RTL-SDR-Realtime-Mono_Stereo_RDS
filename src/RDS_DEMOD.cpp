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
#include "RDS.h"
#include "fmPll.h"


int main(int argc,char** argv)
{
	const int rf_decim = 10;
	int Fc = 3000;
	int upsample = 80;
	int downsample = 19;
	double PLLfreq = 3e3;
	double Fs = 0;
	int k = 0;
	int num_taps = 51;
	if(argc == 1){Fs = 250e3; k = 25;num_taps = num_taps*upsample;}
	else{Fs = 240e3; k = 24;}
	double PLLNCOscale = 2.0;
	double phaseAdjust = 0.0;
	double normBandwidth = 0.01;
	const int block_size = 1024 * rf_decim * downsample * 2;
	int block_count = 0;
	int sample_point = 0;
	std::vector<double> rds_ds;
	std::vector<double> rds_data;
	std::vector<double> ncoOut((int)(block_size)/20, 0);
	std::vector<double> prevstate(6,0);

	std::vector<double> rds_carrier((int)(block_size)/20, 0); // get from before
	std::vector<double> logicdata(1, 0);
	std::vector<double> logicdata_ds;
	std::vector<double>rds_data_ds((int)(block_size)/20, 0);
	std::vector<double> rds_data_RRC((int)(block_size)/20, 0);
	std::vector<double> rds_coeff;
	std::vector<double> state(num_taps-1, 0);
	std::vector<double> state_rds_data((int)(block_size)/20, 0);
	std::vector<double> impulseResponseRRC;

/////// MIXER SHIT - ENDp
while ((block_count+1)*block_size < rds_carrier.size() /*&& iii == 0*/)
{
	std::vector<double> rds_carrier_ds(rds_carrier.begin() + (block_count*block_size)/2, rds_carrier.begin() + ((block_count+1)*block_size)/2);


	fmPll(rds_carrier_ds,PLLfreq,Fs,PLLNCOscale,phaseAdjust,normBandwidth,ncoOut,prevstate);
	for(auto i = 0; i < rds_ds.size();i++){
		rds_ds[i] = ncoOut[i] * rds_carrier[i] * 2;
	}
	/////// MIXER SHIT - END


// low pass + resampler - begin



	impulseResponseLPF( Fs, Fc, num_taps,rds_coeff, 1);
	resampler(upsample,downsample, rds_data_ds, rds_ds, rds_coeff,state_rds_data);

	// low pass + resampler - end

// RRC begin
impulseResponseRootRaisedCosine(Fs, num_taps,impulseResponseRRC);
convolveFIR_N_dec(1, rds_data_RRC, rds_data_ds, impulseResponseRRC, state );

// RRC end

// CDR begin

// to find the best point to sample, take the first 10 symbols and find the largest points, then avarage them
CDR(rds_data_RRC,k, sample_point);


// CDR end

// DATA Processing begin
Manchester_and_differntial(rds_data_RRC,sample_point,logicdata_ds);
// DATA Processing end
logicdata.insert(logicdata.end(), logicdata_ds.begin(), logicdata_ds.end());
block_count++;
}
	return 0;
}
