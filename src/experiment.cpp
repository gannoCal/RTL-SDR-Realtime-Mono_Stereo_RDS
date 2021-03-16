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


void upsampleIQSamples(std::vector<float> &,std::vector<float> &);

int main(int argc,char** argv)
{
	int mode = 0;
	if(argc > 1){
		mode = (int)(*argv[1] & 0xb0001);
	}
	std::cout << "Mode is initially : " << mode << " \n";
// 	// binary files can be generated through the
// 	// Python models from the "../model/" sub-folder
// 	const std::string in_fname = "../data/fm_demod_10.bin";
// 	std::vector<float> bin_data;
// 	readBinData(in_fname, bin_data);

// 	// generate an index vector to be used by logVector on the X axis
// 	std::vector<float> vector_index;
// 	genIndexVector(vector_index, bin_data.size());
// 	// log time data in the "../data/" subfolder in a file with the following name
// 	// note: .dat suffix will be added to the log file in the logVector function
// 	logVector("demod_time", vector_index, bin_data);

// 	// take a slice of data with a limited number of samples for the Fourier transform
// 	// note: NFFT constant is actually just the number of points for the
// 	// Fourier transform - there is no FFT implementation ... yet
// 	// unless you wish to wait for a very long time, keep NFFT at 1024 or below
// 	std::vector<float> slice_data = \
// 		std::vector<float>(bin_data.begin(), bin_data.begin() + NFFT);
// 	// note: make sure that binary data vector is big enough to take the slice

// 	// declare a vector of complex values for DFT
//   std::vector<std::complex<float>> Xf;
// 	// ... in-lab ...
// 	// compute the Fourier transform
// 	// the function is already provided in fourier.cpp

// 	// compute the magnitude of each frequency bin
// 	// note: we are concerned only with the magnitude of the frequency bin
// 	// (there is NO logging of the phase response, at least not at this time)
// 	std::vector<float> Xmag;
// 	// ... in-lab ...
// 	// compute the magnitude of each frequency bin
// 	// the function is already provided in fourier.cpp

// 	// log the frequency magnitude vector
// 	vector_index.clear();
// 	genIndexVector(vector_index, Xmag.size());
// 	logVector("demod_freq", vector_index, Xmag); // log only positive freq

// 	// for your take-home exercise - repeat the above after implementing
// 	// your OWN function for PSD based on the Python code that has been provided
// 	// note the estimate PSD function should use the entire block of "bin_data"
// 	//
// 	// ... complete as part of the take-home ...
// 	//

// 	// if you wish to write some binary files, see below example
// 	// const std::string out_fname = "../data/outdata.bin";
// 	// writeBinData(out_fname, bin_data);

	int upsample = 0;
	int rf_Fs = 0;
	int decimator = 0;
	const int rf_Fc = 100e3;
	const int rf_taps = 151;
	const int rf_decim = 10;

	const int audio_Fs = 48e3;
	const int audio_decim = 5;
	const int audio_Fc = 16e3;
	const int audio_taps = 151;

	if(mode == 1)
	{
		rf_Fs = 2.5e6;
		decimator = 125;
	}
	else
	// Mode 0 assumed to be default if Mode 1 not selected
	{
		rf_Fs = 2.4e6;
		decimator = audio_decim;
	}

	const std::string in_fname = "../data/my_samples_u8.raw";
	std::vector<uint8_t> bin_data;

	readRawData(in_fname, bin_data);

	std::vector<float> iq_data(bin_data.size(), 0);
	for(auto ii = 0 ; ii < bin_data.size() ; ii++){
		iq_data[ii] = (float)(bin_data[ii] - 128) / (float)128;
	}
	std::cout << "iq data : " << iq_data[0] << " \n";
	std::vector<float> rf_coeff;
	std::vector<float> audio_coeff;
	std::vector<float> fm_demod;


	impulseResponseLPF(rf_Fs, rf_Fc, rf_taps, rf_coeff);
	impulseResponseLPF(audio_Fs, audio_Fc, audio_taps, audio_coeff);

	std::cout << "rf_coeff 0 : " << rf_coeff[0] << " \n";
	std::cout << "rf_coeff 1 : " << rf_coeff[1] << " \n";

	const int block_size = 1024 * rf_decim * audio_decim * 2;
	int block_count = 0;
	std::vector<float> state_i_lpf_100k(rf_taps-1, 0);
	std::vector<float> state_q_lpf_100k(rf_taps-1, 0);
	//vector<float> i_filt(rf_taps-1, 0);
	//vector<float> q_filt(rf_taps-1, 0);
	std::vector<float> i_ds((block_size-1)/10, 0);
	std::vector<float> q_ds((block_size-1)/10, 0);
	std::vector<float> data_to_keep(audio_taps-1, 0);

	std::vector<float> audio_ds((block_size-1)/50, 0);
	std::vector<float> state_conv(audio_taps-1, 0);
	std::vector<float> state_phase(2, 0);

	std::vector<float> audio_data;

	std::vector<float> i_samples, q_samples;
	i_samples.resize(iq_data.size()/2);
	q_samples.resize(iq_data.size()/2);

	int sample_counter = 0;
	for (auto i = 0; i < iq_data.size() - 1; i = i + 2)
	{
		i_samples[sample_counter] = iq_data[i];
		q_samples[sample_counter] = iq_data[i+1];
		sample_counter++;
	}
	std::cout << "i data 0 : " << i_samples[0] << " \n";
	std::cout << "q data 0: " << q_samples[0] << " \n";
	// Upsampling the IQ samples if mode1 was selected
	if(mode == 1) {std::cout << "Mode is 1\n"; upsampleIQSamples(i_samples, q_samples);}
	int iii = 0;
	while ((block_count+1)*block_size < iq_data.size() /*&& iii == 0*/)
	{
		std::cout << "IN while loop\n";
		// Seperate the necessary I/Q samples for this block
		std::vector<float> i_samples_block(i_samples.begin() + (block_count*block_size), i_samples.begin() + ((block_count+1)*block_size));
		std::vector<float> q_samples_block(q_samples.begin() + (block_count*block_size), q_samples.begin() + ((block_count+1)*block_size));
		
		std::cout << "input size : " <<  i_samples_block.size() << " \n";
		std::cout << "output size : " <<  i_ds.size() << " \n";
		//std::cout << "i block rand# : " <<  i_samples_block[45] << " \n";

		// Next step -- grab every second value for I grab every other value for Q
		convolveFIR_N_dec(10, i_ds, i_samples_block, rf_coeff,state_i_lpf_100k);
		std::cout << "output size post : " <<  i_ds.size() << " \n";
		convolveFIR_N_dec(10, q_ds, q_samples_block, rf_coeff,state_q_lpf_100k);
		std::cout << "ids data 0 : " << i_ds[0] << " \n";
		std::cout << "qds data 0: " << q_ds[0] << " \n";
		std::cout << "ids data 1 : " << i_ds[1] << " \n";
		std::cout << "qds data 1: " << q_ds[1] << " \n";
		fmDemodArctanBlock(fm_demod,i_ds, q_ds, state_phase);
		std::cout << "output size tan: " <<  i_ds.size() << " \n";
		std::cout << "size: " << fm_demod.size() << " \n";
		for(auto jjj = 0 ; jjj < fm_demod.size() ; jjj ++){
			//std::cout << "demod block post-convolution rand# : " <<  fm_demod[jjj] << " \n";
		// if(std::isnan(fm_demod[jjj])){
		// 	std::cout << "New Bad Samples" << fm_demod[jjj] << " at : " << jjj << " \n";
		// }
		}

		if(block_count == (iq_data.size()/block_size)/2){
			std::vector<float> vector_index;
			genIndexVector(vector_index, 256);
			// log time data in the "../data/" subfolder in a file with the following name
			// note: .dat suffix will be added to the log file in the logVector function
			std::vector<float> freq, psd_est;

			float Fs = 240;
			int NFFT_in = (int) NFFT;

			estimatePSD(fm_demod, NFFT_in, Fs, freq, psd_est);
			std::cout << "\nCalculated PSD\n";
			logVector("demod_psd", vector_index, psd_est);
			std::cout << "Generated PSD log" <<   " \n";
		}



		convolveFIR_N_dec(decimator, audio_ds,fm_demod,audio_coeff,state_conv);
		std::cout << "LAst Convolve Good!" <<   " \n";
		


		audio_data.insert(audio_data.end(), audio_ds.begin(), audio_ds.end());
		iii = 1;
		block_count++;
	}

	// naturally, you can comment the line below once you are comfortable to run gnuplot
	std::cout << "Run: gnuplot -e 'set terminal png size 1024,768' example.gnuplot > ../data/example.png\n";


	writeBinData("../data/new_binary_test.bin",audio_data);

	return 0;
}

void upsampleIQSamples(std::vector<float> &i_samples,std::vector<float> &q_samples)
{
	std::vector<float> i_samples_upsampled, q_samples_upsampled;
	auto upsampledLength = i_samples.size() * 24;

	i_samples_upsampled.resize(upsampledLength);
	q_samples_upsampled.resize(upsampledLength);
	auto sample_counter = 0;
	for(auto ii = 0; ii < i_samples_upsampled.size(); ii++){
		if(ii%24)
		{
			i_samples_upsampled[ii] = i_samples[sample_counter];
			q_samples_upsampled[ii] = q_samples[sample_counter];
			sample_counter++;
		}
		else
		{
			i_samples_upsampled[ii] = 0;
			q_samples_upsampled[ii] = 0;
		}
	}

	// Reassigning/resizing i_samples vector and q_samples vector to re-use some mode0 code
	i_samples.resize(upsampledLength);
	q_samples.resize(upsampledLength);

	for(auto i = 0; i < i_samples.size();i++)
	{
		i_samples[i] = i_samples_upsampled[i];
		q_samples[i] = q_samples_upsampled[i];
	}
}