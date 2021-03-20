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


void upsampleIQSamples(std::vector<double> &,std::vector<double> &);

int main(int argc,char** argv)
{
	int mode = 0;
	if(argc > 1){
		mode = (int)(*argv[1] & 0xb0001);
	}
	std::cout << "Mode is initially : " << mode << " \n";
mode = 1;
	int upsample = 0;
	int rf_Fs = 2.4e6;
	int decimator = 5;
	const int rf_Fc = 100e3;
	const int rf_taps = 151;
	const int rf_decim = 10;

	int audio_Fs = 240e3;
	const int audio_decim = 5;
	const int audio_Fc = 16e3;
	const int audio_taps = 151;

	if(mode == 1)
	{
		rf_Fs = 2.5e6;
		decimator = 125;
        audio_Fs = 250e3;
	}
	else
	// Mode 0 assumed to be default if Mode 1 not selected
	{
		rf_Fs = 2.4e6;
		decimator = audio_decim;
	}

	const std::string in_fname = "../data/my_samples_2.4_march16.raw";
	std::vector<uint8_t> bin_data;

	readRawData(in_fname, bin_data);

	std::vector<double> iq_data(bin_data.size(), (double)0.0);
	for(auto ii = 0 ; ii < bin_data.size() ; ii++){
		iq_data[ii] = ((double)bin_data[ii] - (double)128.0) / (double)128.0;
	//	std::cout << "bin_data " << iq_data[ii] <<  "\n";
	}
	std::cout << "iq data : " << iq_data[0] << " \n";
	std::vector<double> rf_coeff;
	std::vector<double> audio_coeff;
	std::vector<double> fm_demod;


	impulseResponseLPF(rf_Fs, rf_Fc, rf_taps, rf_coeff,1);
	impulseResponseLPF(audio_Fs, audio_Fc, audio_taps, audio_coeff,1);

	std::cout << "rf_coeff 0 : " << rf_coeff[0] << " \n";
	std::cout << "rf_coeff 1 : " << rf_coeff[1] << " \n";

	const int block_size = 1024 * rf_decim * audio_decim * 2;
	int block_count = 0;
	std::vector<double> state_i_lpf_100k(rf_taps-1, 0);
	std::vector<double> state_q_lpf_100k(rf_taps-1, 0);
	//vector<float> i_filt(rf_taps-1, 0);
	//vector<float> q_filt(rf_taps-1, 0);
	std::vector<double> i_ds((block_size)/20, 0);
	std::vector<double> q_ds((block_size)/20, 0);

	std::vector<double> audio_ds(((block_size)/100), 0);
	std::vector<double> state_conv(audio_taps-1, 0);
	std::vector<double> state_phase(2, 0.0);

	std::vector<double> audio_data;

	std::vector<double> i_samples, q_samples;
	i_samples.resize(iq_data.size()/2);
	q_samples.resize(iq_data.size()/2);

	int sample_counter = 0;
	for (auto i = 0; i < iq_data.size() - 1; i = i + 2)
	{
		i_samples[sample_counter] = (double)iq_data[i];
		q_samples[sample_counter] = (double)iq_data[i+1];
		sample_counter++;
	//std::cout << "iq_data " << iq_data[i] << " i_samples " << i_samples[sample_counter] << " q_samples: " << q_samples[sample_counter] << " \n";
	}



	std::cout << "i data 0 : " << i_samples[0] << " \n";
	std::cout << "q data 0: " << q_samples[0] << " \n";
    int outputBlock = 1;
	// Upsampling the IQ samples if mode1 was selected
	while ((block_count+1)*block_size < iq_data.size() /*&& iii == 0*/)
	{
		std::cout << "IN while loop\n";
		// Seperate the necessary I/Q samples for this block
		std::vector<double> i_samples_block(i_samples.begin() + (block_count*block_size)/2, i_samples.begin() + ((block_count+1)*block_size)/2);
		std::vector<double> q_samples_block(q_samples.begin() + (block_count*block_size)/2, q_samples.begin() + ((block_count+1)*block_size)/2);	///for speed - declare these outside loop


		if(block_count == (outputBlock)){
			std::vector<float> vector_index;
            // genIndexVector(vector_index, iq_data.size());

			// logVector("1_iq_data", vector_index, iq_data);


			genIndexVector(vector_index, i_samples_block.size());

			logVector("1_i_samples_block0", vector_index, i_samples_block);
			logVector("1_q_samples_block0", vector_index, q_samples_block);

			genIndexVector(vector_index, state_i_lpf_100k.size());
			logVector("1_BEFORE_state_i_lpf_100k_block0", vector_index, state_i_lpf_100k);

			genIndexVector(vector_index, state_q_lpf_100k.size());
			logVector("1_BEFORE_state_q_lpf_100k_block0", vector_index, state_q_lpf_100k);

			genIndexVector(vector_index, state_phase.size());

			logVector("1_BEFORE_state_phase_block0", vector_index, state_phase);

			genIndexVector(vector_index, state_conv.size());

			logVector("1_BEFORE_state_conv_block0", vector_index, state_conv);

		}


		// std::cout << "input size : " <<  i_samples_block.size() << " \n";
		// std::cout << "output size : " <<  i_ds.size() << " \n";
		//std::cout << "i block rand# : " <<  i_samples_block[45] << " \n";

		// Next step -- grab every second value for I grab every other value for Q
		convolveFIR_N_dec(10, i_ds, i_samples_block, rf_coeff,state_i_lpf_100k);
		// std::cout << "output size post : " <<  i_ds.size() << " \n";
		convolveFIR_N_dec(10, q_ds, q_samples_block, rf_coeff,state_q_lpf_100k);
		// std::cout << "ids data 0 : " << i_ds[0] << " \n";
		// std::cout << "qds data 0: " << q_ds[0] << " \n";
		// std::cout << "ids data 1 : " << i_ds[1] << " \n";
		// std::cout << "qds data 1: " << q_ds[1] << " \n";

		if(block_count == (outputBlock /*iq_data.size()/block_size)/2*/)){
			std::vector<float> vector_index;
			genIndexVector(vector_index, i_ds.size());

			logVector("1_i_ds_block0", vector_index, i_ds);
			logVector("1_q_ds_block0", vector_index, q_ds);

			genIndexVector(vector_index, rf_coeff.size());
			logVector("1_rf_coeff_block0", vector_index, rf_coeff);

			genIndexVector(vector_index, state_i_lpf_100k.size());
			logVector("1_AFTER_state_i_lpf_100k_block0", vector_index, state_i_lpf_100k);

			genIndexVector(vector_index, state_q_lpf_100k.size());
			logVector("1_AFTER_state_q_lpf_100k_block0", vector_index, state_q_lpf_100k);


		}

		fmDemodArctanBlock(fm_demod,i_ds, q_ds, state_phase);
		std::cout << "output size tan: " <<  i_ds.size() << " \n";
		std::cout << "size: " << fm_demod.size() << " \n";
if(mode == 1){std::cout << "MODE 1 \n";convolve_UPSAMPLE_N_dec(5,24, audio_ds,fm_demod,audio_coeff,state_conv);}
else{convolveFIR_N_dec(decimator, audio_ds,fm_demod,audio_coeff,state_conv);}
		// for(auto jjj = 0 ; jjj < fm_demod.size() ; jjj ++){
		// 	//std::cout << "demod block post-convolution rand# : " <<  fm_demod[jjj] << " \n";
		// // if(std::isnan(fm_demod[jjj])){
		// // 	std::cout << "New Bad Samples" << fm_demod[jjj] << " at : " << jjj << " \n";
		// // }
		// }
		if(block_count == (outputBlock /*iq_data.size()/block_size)/2*/)){
			std::vector<float> vector_index;
			genIndexVector(vector_index, fm_demod.size());

			logVector("1_fm_demod_block0", vector_index, fm_demod);

			genIndexVector(vector_index, state_phase.size());

			logVector("1_AFTER_state_phase_block0", vector_index, state_phase);
		}




		// if(block_count == (outputBlock)){
		// 	std::vector<float> vector_index;
		// 	genIndexVector(vector_index, 256);
		// 	// log time data in the "../data/" subfolder in a file with the following name
		// 	// note: .dat suffix will be added to the log file in the logVector function
		// 	std::vector<float> freq, psd_est;

		// 	float Fs = 240;
		// 	int NFFT_in = (int) NFFT;

		// 	estimatePSD(fm_demod, NFFT_in, Fs, freq, psd_est);
		// 	std::cout << "\nCalculated PSD\n";
		// 	logVector("demod_psd", vector_index, psd_est);
		// 	std::cout << "Generated PSD log" <<   " \n";
		// }



	//	std::cout << "LAst Convolve Good!" <<   " \n";


		if(block_count == (outputBlock /*iq_data.size()/block_size)/2*/)){
			std::vector<float> vector_index;
			genIndexVector(vector_index, audio_ds.size());

			logVector("1_audio_ds_block0", vector_index, audio_ds);

			genIndexVector(vector_index, audio_coeff.size());

			logVector("1_audio_coeff_block0", vector_index, audio_coeff);

			genIndexVector(vector_index, state_conv.size());

			logVector("1_AFTER_state_conv_block0", vector_index, state_conv);
		}



        if(block_count == (outputBlock)){
			std::vector<float> vector_index;
			genIndexVector(vector_index, 256);
			// log time data in the "../data/" subfolder in a file with the following name
			// note: .dat suffix will be added to the log file in the logVector function
			std::vector<double> freq, psd_est;

			double Fs = 240;
			int NFFT_in = (int) NFFT;

			estimatePSD(audio_ds, NFFT_in, Fs, freq, psd_est);
			std::cout << "\nCalculated PSD\n";
			logVector("demod_psd", vector_index, psd_est);
			std::cout << "Generated PSD log" <<   " \n";
		}


        // for(auto i = 0 ; i < audio_ds.size() ; i++){
        //     audio_ds[i] = audio_ds[i] * 500;
        // }


		audio_data.insert(audio_data.end(), audio_ds.begin(), audio_ds.end());
		if(block_count == (outputBlock /*iq_data.size()/block_size)/2*/)){
			std::vector<float> vector_index;
			genIndexVector(vector_index, audio_data.size());

			logVector("1_audio_data_block0", vector_index, audio_data);
		}
		block_count++;
	}

	// naturally, you can comment the line below once you are comfortable to run gnuplot
	std::cout << "Run: gnuplot -e 'set terminal png size 1024,768' example.gnuplot > ../data/example.png\n";

	for(auto i = 0; i<audio_data.size();i++){std::cout << audio_data[i] << " \n";}
	writeBinData("../data/new_binary69_test6.bin",audio_data);

	return 0;
}
