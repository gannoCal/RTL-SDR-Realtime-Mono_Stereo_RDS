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

int main()
{
	// binary files can be generated through the
	// Python models from the "../model/" sub-folder
	const std::string in_fname = "../data/fm_demod_10.bin";
	std::vector<float> bin_data;
	readBinData(in_fname, bin_data);

	// generate an index vector to be used by logVector on the X axis
	std::vector<float> vector_index;
	genIndexVector(vector_index, bin_data.size());
	// log time data in the "../data/" subfolder in a file with the following name
	// note: .dat suffix will be added to the log file in the logVector function
	logVector("demod_time", vector_index, bin_data);

	// take a slice of data with a limited number of samples for the Fourier transform
	// note: NFFT constant is actually just the number of points for the
	// Fourier transform - there is no FFT implementation ... yet
	// unless you wish to wait for a very long time, keep NFFT at 1024 or below
	std::vector<float> slice_data = \
		std::vector<float>(bin_data.begin(), bin_data.begin() + NFFT);
	// note: make sure that binary data vector is big enough to take the slice

	// declare a vector of complex values for DFT
  std::vector<std::complex<float>> Xf;
	// ... in-lab ...
	// compute the Fourier transform
	// the function is already provided in fourier.cpp

	// compute the magnitude of each frequency bin
	// note: we are concerned only with the magnitude of the frequency bin
	// (there is NO logging of the phase response, at least not at this time)
	std::vector<float> Xmag;
	// ... in-lab ...
	// compute the magnitude of each frequency bin
	// the function is already provided in fourier.cpp

	// log the frequency magnitude vector
	vector_index.clear();
	genIndexVector(vector_index, Xmag.size());
	logVector("demod_freq", vector_index, Xmag); // log only positive freq

	// for your take-home exercise - repeat the above after implementing
	// your OWN function for PSD based on the Python code that has been provided
	// note the estimate PSD function should use the entire block of "bin_data"
	//
	// ... complete as part of the take-home ...
	//

	// if you wish to write some binary files, see below example
	// const std::string out_fname = "../data/outdata.bin";
	// writeBinData(out_fname, bin_data);
	const int rf_Fs = 2.4e6;
	const int rf_Fc = 100e3;
	const int rf_taps = 151;
	const int rf_decim = 10;

	const int audio_Fs = 48e3;
	const int audio_decim = 5;
	const int audio_Fc = 16e3;
	const int audio_taps = 151;

	const string in_fname = "../data/iq_samples.raw";
	std::vector<float> bin_data;
	readBinData(in_fname, &bin_data)

	std::vector<float> rf_coeff;
	std::vector<float> fm_demod;

	impulseResponseLPF(audio_Fs, audio_Fc, audio_taps, rf_coeff);


	const int block_size = 1024 * rf_decim * audio_decim * 2;
	int block_count = 0;
	vector<float> state_i_lpf_100k(rf_taps-1, 0);
	vector<float> state_q_lpf_100k(rf_taps-1, 0);
	vector<float> i_filt(rf_taps-1, 0);
	vector<float> q_filt(rf_taps-1, 0);
	vector<float> i_ds;
	vector<float> q_ds;
	vector<float> data_to_keep(audio_taps-1, 0);
	vector<float> state_conv(audio_taps-1, 0);
	vector<float> state_phase(2, 0);

	while ((block_count+1)*block_size < len(iq_data)){
// still need to change this
		int convolveSpacing;//To set the interval on which we convolve

		//convolveFIR(i_filt, iq_data[(block_count)*block_size:(block_count+1)*block_size:2],rf_coeff)
		//convolveFIR(q_filt,iq_data[(block_count)*block_size+1:(block_count+1)*block_size:2],rf_coeff)
		convolveFIR_N_step(convolveSpacing, i_filt, iq_data[(block_count)*block_size:(block_count+1)*block_size:2],rf_coeff)
		convolveFIR_N_step(convolveSpacing, q_filt,iq_data[(block_count)*block_size+1:(block_count+1)*block_size:2],rf_coeff)

		i_ds.resize(int(i_filt/5), 0.0);
		q_ds.resize(int(q_filt/5), 0.0);

		int j = 0;

		for(auto i = 0; i < i_filt.size();i++){
			//if(i%rf_decim == 0){
			if(i%convolveSpacing == 0){	
				i_ds[j] = i_filt[i];
				q_ds[j] = q_filt[i];
				j++;

			}
		}


		fmDemodArctanBlock(fm_demod,i_ds, q_ds, state_phase)


		audio_filt, state_conv = conv(audio_coeff, fm_demod, state_conv)
		convolveFIR(audio_filt,	audio_coeff,fm_demod)

		 j = 0;
		for(auto i = 0; i < i_filt.size();i++){
			if(i%rf_decim == 0){
				i_ds[j] = i_filt[i];
				q_ds[j] = q_filt[i];
				j++;

			}
		}

		audio_data.insert(audio_data.end(), audio_block.begin(), audio_block.end());

}

	// naturally, you can comment the line below once you are comfortable to run gnuplot
	std::cout << "Run: gnuplot -e 'set terminal png size 1024,768' example.gnuplot > ../data/example.png\n";

	return 0;
}
