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

int main(int argc,char* argv[])
{
	int mode = 0;
    if(argc < 2){
        std::cerr << "Default - Mode 0";
    
	}else if(argc == 2){
		mode = atoi(argv[1]);
        if(mode != 1 && mode != 0){
            std::cerr << "Wrong Mode! Exiting...";
            exit(1);
        }
        std::cerr << "Mode is : " << mode << " \n";
	}else{
        std::cerr << "Error. Please fix your input. " << " \n";
    }


    


	

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

	std::vector<double> rf_coeff;
	std::vector<double> audio_coeff;
	std::vector<double> fm_demod;


	impulseResponseLPF(rf_Fs, rf_Fc, rf_taps, rf_coeff,1);
	impulseResponseLPF(audio_Fs, audio_Fc, audio_taps, audio_coeff,1);

	const int block_size = 1024 * rf_decim * audio_decim * 2;
	int block_count = 0;
	std::vector<double> state_i_lpf_100k(rf_taps-1, 0);
	std::vector<double> state_q_lpf_100k(rf_taps-1, 0);
	//vector<float> i_filt(rf_taps-1, 0);
	//vector<float> q_filt(rf_taps-1, 0);
	std::vector<double> i_ds((block_size)/20, 0);
	std::vector<double> q_ds((block_size)/20, 0);

	std::vector<double> audio_ds((block_size)/100, 0);
	std::vector<double> state_conv(audio_taps-1, 0);
	std::vector<double> state_phase(2, 0.0);

	

    int outputBlock = 1;
	// Upsampling the IQ samples if mode1 was selected
    std::vector<short int> audio_data(block_size/100);
	std::cerr << "Enter Loop"  << " \n";
        for (unsigned int block_id = 0; ; block_id++){

        std::vector<double> iq_data(block_size);

        readStdinBlockData(block_size,block_id,iq_data);
        if((std::cin.rdstate())!=0){
            std::cerr << "End of Input Stream";
            exit(1);
        }

        std::vector<double> i_samples_block, q_samples_block;
        i_samples_block.resize(iq_data.size()/2,(double)0.0);
        q_samples_block.resize(iq_data.size()/2,(double)0.0);

        int sample_counter = 0;
        for (auto i = 0; i < iq_data.size() - 1; i = i + 2)
        {
            i_samples_block[sample_counter] = iq_data[i];
            q_samples_block[sample_counter] = iq_data[i+1];
            sample_counter++;
        }
        if(mode == 1) {
        std::cout << "Mode is 1\n"; upsampleIQSamples(i_samples_block, q_samples_block);
        }
        
		
		// if(block_id == (outputBlock)){
		// 	std::vector<float> vector_index;
            


		// 	genIndexVector(vector_index, i_samples_block.size());
			
		// 	logVector("1_i_samples_block0", vector_index, i_samples_block);
		// 	logVector("1_q_samples_block0", vector_index, q_samples_block);

		// 	genIndexVector(vector_index, state_i_lpf_100k.size());
		// 	logVector("1_BEFORE_state_i_lpf_100k_block0", vector_index, state_i_lpf_100k);

		// 	genIndexVector(vector_index, state_q_lpf_100k.size());
		// 	logVector("1_BEFORE_state_q_lpf_100k_block0", vector_index, state_q_lpf_100k);

		// 	genIndexVector(vector_index, state_phase.size());
			
		// 	logVector("1_BEFORE_state_phase_block0", vector_index, state_phase);

		// 	genIndexVector(vector_index, state_conv.size());
			
		// 	logVector("1_BEFORE_state_conv_block0", vector_index, state_conv);
			
		// }
        
		convolveFIR_N_dec(10, i_ds, i_samples_block, rf_coeff,state_i_lpf_100k);
		convolveFIR_N_dec(10, q_ds, q_samples_block, rf_coeff,state_q_lpf_100k);

		// if(block_id == (outputBlock)){
		// 	std::vector<float> vector_index;
		// 	genIndexVector(vector_index, i_ds.size());
			
		// 	logVector("1_i_ds_block0", vector_index, i_ds);
		// 	logVector("1_q_ds_block0", vector_index, q_ds);

		// 	genIndexVector(vector_index, rf_coeff.size());
		// 	logVector("1_rf_coeff_block0", vector_index, rf_coeff);

		// 	genIndexVector(vector_index, state_i_lpf_100k.size());
		// 	logVector("1_AFTER_state_i_lpf_100k_block0", vector_index, state_i_lpf_100k);

		// 	genIndexVector(vector_index, state_q_lpf_100k.size());
		// 	logVector("1_AFTER_state_q_lpf_100k_block0", vector_index, state_q_lpf_100k);

			
		// }

		fmDemodArctanBlock(fm_demod,i_ds, q_ds, state_phase);
		// if(block_id == (outputBlock)){
		// 	std::vector<float> vector_index;
		// 	genIndexVector(vector_index, fm_demod.size());
			
		// 	logVector("1_fm_demod_block0", vector_index, fm_demod);

		// 	genIndexVector(vector_index, state_phase.size());
			
		// 	logVector("1_AFTER_state_phase_block0", vector_index, state_phase);
		// }




		// if(block_id == (outputBlock)){
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



		convolveFIR_N_dec(decimator, audio_ds,fm_demod,audio_coeff,state_conv);

		// if(block_id == (outputBlock /*iq_data.size()/block_size)/2*/)){
		// 	std::vector<float> vector_index;
		// 	genIndexVector(vector_index, audio_ds.size());
			
		// 	logVector("1_audio_ds_block0", vector_index, audio_ds);

		// 	genIndexVector(vector_index, audio_coeff.size());
			
		// 	logVector("1_audio_coeff_block0", vector_index, audio_coeff);

		// 	genIndexVector(vector_index, state_conv.size());
			
		// 	logVector("1_AFTER_state_conv_block0", vector_index, state_conv);
		// }



        // if(block_id == (outputBlock)){
		// 	std::vector<float> vector_index;
		// 	genIndexVector(vector_index, 256);
		// 	// log time data in the "../data/" subfolder in a file with the following name
		// 	// note: .dat suffix will be added to the log file in the logVector function
		// 	std::vector<double> freq, psd_est;

		// 	double Fs = 240;
		// 	int NFFT_in = (int) NFFT;

		// 	estimatePSD(audio_ds, NFFT_in, Fs, freq, psd_est);
		// 	std::cout << "\nCalculated PSD\n";
		// 	logVector("demod_psd", vector_index, psd_est);
		// 	std::cout << "Generated PSD log" <<   " \n";
		// }



		// //audio_data.insert(audio_data.end(), audio_ds.begin(), audio_ds.end());
        for(unsigned int k=0 ; k < audio_ds.size() ; k++){
            if(std::isnan(audio_ds[k])) audio_data[k] = 0;
            else audio_data[k] = audio_ds[k] * 16384;
        }
        fwrite(&audio_data[0], sizeof(short int),audio_data.size(),stdout);
		// if(block_id == (outputBlock)){
		// 	std::vector<float> vector_index;
		// 	genIndexVector(vector_index, audio_data.size());
			
		// 	logVector("1_audio_data_block0", vector_index, audio_data);
		// }
}

	// naturally, you can comment the line below once you are comfortable to run gnuplot
	std::cout << "Run: gnuplot -e 'set terminal png size 1024,768' example.gnuplot > ../data/example.png\n";


	//writeBinData("../data/new_binary69_test6.bin",audio_data);

	return 0;
}

void upsampleIQSamples(std::vector<double> &i_samples,std::vector<double> &q_samples)
{
	std::vector<double> i_samples_clone(i_samples.size() , 0.0); 
    std::vector<double> q_samples_clone(i_samples.size(), 0.0);
	for(auto ii = 0; ii < i_samples.size(); ii++){
			i_samples_clone[ii] = i_samples[ii];
			q_samples_clone[ii] = q_samples[ii];
	}

	// Reassigning/resizing i_samples vector and q_samples vector to re-use some mode0 code
	i_samples.resize(i_samples.size() * 24 , 0.0);
	q_samples.resize(i_samples.size() * 24 , 0.0);

    auto sample_counter = 0;
	for(auto i = 0; i < i_samples.size()*24;i = i + 24)
	{
		i_samples[i] = i_samples_clone[sample_counter];
		q_samples[i] = q_samples_clone[sample_counter];
        sample_counter++;
	}
}