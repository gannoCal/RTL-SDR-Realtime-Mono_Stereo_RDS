/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Copyright by Nicola Nicolici
Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"
#include "fmPll.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"



void upsampleIQSamples(std::vector<double> &,std::vector<double> &);

int main(int argc,char* argv[])
{
    bool floatMode = false;
    bool stereoMode = false;
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
	}else if(argc == 3){
        mode = atoi(argv[1]);
        char* ModeIn = argv[2];
        if(mode != 1 && mode != 0){
            std::cerr << "Wrong Mode! Exiting...";
            exit(1);
        }

        int len = sizeof(ModeIn)/sizeof(*ModeIn);

        if(len < 2){
            std::cerr << "Bad Param! Exiting...";
            exit(1);
        }

        if(ModeIn[0] == '-' && ModeIn[1]== 'f'){
            floatMode = true;
        }else if(ModeIn[0] == '-' && ModeIn[1]== 's'){
            stereoMode = true;
        }else{
            std::cerr << "Bad Param! Exiting...";
            exit(1);
        }
        std::cerr << "Mode is : " << mode << " \n";
        std::cerr << "Param is : " << ModeIn[0] << ModeIn[1] << " \n";

    }else if(argc == 4){
        mode = atoi(argv[1]);
        char* ModeIn1 = argv[2];
        char* ModeIn2 = argv[3];
        if(mode != 1 && mode != 0){
            std::cerr << "Wrong Mode! Exiting...";
            exit(1);
        }
        int len1 = sizeof(ModeIn1)/sizeof(*ModeIn1);
        int len2 = sizeof(ModeIn2)/sizeof(*ModeIn2);

        if(len1 < 2 || len2 < 2){
            std::cerr << "Bad Param! Exiting...";
            exit(1);
        }

        if( ModeIn1[0] == ModeIn2[0] && ModeIn1[1] == ModeIn2[1]){
            std::cerr << "Duplicate param! Exiting...";
            exit(1);
        }

        if(ModeIn1[0] == '-' && ModeIn1[1]== 'f'){
            floatMode = true;
        }else if(ModeIn1[0] == '-' && ModeIn1[1]== 's'){
            stereoMode = true;
        }else{
            std::cerr << "Bad Param! Exiting...";
            exit(1);
        }

        if(ModeIn2[0] == '-' && ModeIn2[1]== 'f'){
            floatMode = true;
        }else if(ModeIn2[0] == '-' && ModeIn2[1]== 's'){
            stereoMode = true;
        }else{
            std::cerr << "Bad Param! Exiting...";
            exit(1);
        }


        

        std::cerr << "Mode is : " << mode << " \n";
        std::cerr << "Param1 is : " << ModeIn1[0] << ModeIn1[1] << " \n";
        std::cerr << "Param2 is : " << ModeIn2[0] << ModeIn2[1] << " \n";

    }else{
        std::cerr << "Error. Please fix your input. " << " \n";
    }


    


	

	int upsample = 0;
	int rf_Fs = 2.4e6;
	int decimator = 5;
	const int rf_Fc = 100e3;
	const int rf_taps = 51;
	const int rf_decim = 10;

	int audio_Fs = 240e3;
	const int audio_decim = 5;
	const int audio_Fc = 16e3;
	const int audio_taps = 51;

    int if_Fs = 240e3;
	const int sbanddecim = 1;
	const int sband_Fb = 22e3;
    const int sband_Fe = 54e3;
	const int sband_taps = 51;
    
	const int tbanddecim = 1;
	const int tband_Fb = 18.5e3;
    const int tband_Fe = 19.5e3;
	const int tband_taps = 51;

	if(mode == 1)
	{
		rf_Fs = 2.5e6;
		decimator = 125;
        audio_Fs = 250e3;
        if_Fs = 250e3;
	}
	else
	// Mode 0 assumed to be default if Mode 1 not selected
	{
		rf_Fs = 2.4e6;
		decimator = audio_decim;
	}

	std::vector<double> rf_coeff;
	std::vector<double> audio_coeff;
    std::vector<double> st_coeff;
    std::vector<double> tone_coeff;
	std::vector<double> fm_demod;


	impulseResponseLPF(rf_Fs, rf_Fc, rf_taps, rf_coeff,1);
	impulseResponseLPF(audio_Fs, audio_Fc, audio_taps, audio_coeff,1);
    impulseResponseBPF(if_Fs, sband_Fb, sband_Fe, audio_taps, st_coeff,1);
    impulseResponseBPF(if_Fs, tband_Fb, tband_Fe, audio_taps, tone_coeff,1);

    std::vector<float> vector_index;
    genIndexVector(vector_index, st_coeff.size());
    logVector("2_st_coeff", vector_index, st_coeff);

    genIndexVector(vector_index, tone_coeff.size());
    logVector("2_tone_coeff", vector_index, tone_coeff);

	const int block_size = 1024 * rf_decim * audio_decim * 2;
	int block_count = 0;
	std::vector<double> state_i_lpf_100k(rf_taps-1, 0);
	std::vector<double> state_q_lpf_100k(rf_taps-1, 0);
	//vector<float> i_filt(rf_taps-1, 0);
	//vector<float> q_filt(rf_taps-1, 0);
	std::vector<double> i_ds((block_size)/20, 0);
	std::vector<double> q_ds((block_size)/20, 0);

    std::vector<double> st_ds((block_size)/20, 0);
    std::vector<double> tone_ds((block_size)/20, 0);

	std::vector<double> audio_ds((block_size)/100, 0);
	std::vector<double> state_conv(audio_taps-1, 0);
	std::vector<double> state_phase(2, 0.0);
    
    std::vector<double> state_st_240k(audio_taps-1, 0);
	std::vector<double> state_tone_240k(audio_taps-1, 0);

    std::vector<double> ncoOut((block_size)/20, 0);
    std::vector<double> stereo_data((block_size)/20, 0);
    std::vector<double> prevstate(6,0);
    prevstate[0] = 0.0;            //state_integrator = 0.0;
    prevstate[1] = 0.0;            //state_phaseEst = 0.0;
    prevstate[2] = 1.0;            //state_feedbackI = 1.0;
    prevstate[3] = 0.0;            //state_feedbackQ = 0.0;
    prevstate[4] = 1.0;            //state_ncoOut0 = 1.0;
    prevstate[5] = 0.0;            //state_trigOffset = 0.0;

    double PLLfreq = 19e3;
    double PLLfs = 240e3;
    double PLLNCOscale = 2.0;
    double phaseAdjust = 0.0;
    double normBandwidth = 0.01;

    std::vector<double> stereo_data_ds((block_size)/100, 0);
    std::vector<short int> stereo_block((block_size)/50, 0);
    std::vector<double> state_stereo_data(audio_taps-1, 0);

	
    //std::vector<short int> stereo_out_data(block_size/50);
    std::vector<short int> audio_data(block_size/100);
    std::vector<double> iq_data(block_size);
    std::vector<double> i_samples_block, q_samples_block;
    i_samples_block.resize(iq_data.size()/2,(double)0.0);
    q_samples_block.resize(iq_data.size()/2,(double)0.0);
    int sample_counter = 0;


	std::cerr << "Enter Loop"  << " \n";
        for (unsigned int block_id = 0; ; block_id++){          //Loop Begin

        
        if(floatMode == false)
            readStdinBlockData(block_size,block_id,iq_data);
        else
            readStdinBlockDataFloat(block_size,block_id,iq_data);

        if((std::cin.rdstate())!=0){
            std::cerr << "End of Input Stream";
            exit(1);
        }

        

        sample_counter = 0;
        for (auto i = 0; i < iq_data.size() - 1; i = i + 2)
        {
            i_samples_block[sample_counter] = iq_data[i];
            q_samples_block[sample_counter] = iq_data[i+1];
            sample_counter++;
        }
        if(mode == 1) {
        std::cout << "Mode is 1\n"; upsampleIQSamples(i_samples_block, q_samples_block);
        }
        
        
		convolveFIR_N_dec(10, i_ds, i_samples_block, rf_coeff,state_i_lpf_100k);
		convolveFIR_N_dec(10, q_ds, q_samples_block, rf_coeff,state_q_lpf_100k);

		
		fmDemodArctanBlock(fm_demod,i_ds, q_ds, state_phase);
		
		convolveFIR_N_dec(decimator, audio_ds,fm_demod,audio_coeff,state_conv);


        if (stereoMode){
            convolveFIR_N_dec(1, st_ds, fm_demod, st_coeff,state_st_240k);
            convolveFIR_N_dec(1, tone_ds, fm_demod, tone_coeff,state_tone_240k);
            
            fmPll(tone_ds,PLLfreq,PLLfs,PLLNCOscale,phaseAdjust,normBandwidth,ncoOut,prevstate);
            
            for(auto i = 0 ; i < st_ds.size();i++){
                stereo_data[i] = ncoOut[i] * st_ds[i] * 2;
            }

            convolveFIR_N_dec(decimator, stereo_data_ds, stereo_data, audio_coeff,state_stereo_data);

            for(unsigned int i=0 ; i < audio_ds.size() ; i++){
            
            stereo_block[2*i+1] = std::isnan(audio_ds[i]) || std::isnan(stereo_data_ds[i]) ? (short int)0 :  (short int)(16384*(audio_ds[i] - stereo_data_ds[i]  ) / 2)   ;     //r
            stereo_block[2*i] = std::isnan(audio_ds[i]) || std::isnan(stereo_data_ds[i]) ? (short int)0 : (short int)(16384*(audio_ds[i]  + stereo_data_ds[i]  ) / 2)    ;   //l
            
            }
            // for(unsigned int k=0 ; k < stereo_block.size() ; k++){
            //     if(std::isnan(stereo_block[k])) stereo_out_data[k] = 0;
            //     else stereo_out_data[k] = (short int)(stereo_block[k] * 16384);
            // }
            fwrite(&stereo_block[0], sizeof(short int),stereo_block.size(),stdout);
            //fwrite(&stereo_data[0], sizeof(short int),stereo_data.size(),stdout);
            // if(block_id == 1){
            //     int NFFT_in = (int)NFFT;
            //     std::vector<double> freq, psd_est;
            //     double fsIN = (double)48e3;
            //     estimatePSD(audio_ds, NFFT_in , fsIN, freq, psd_est);

	
	        //     std::vector<float> vector_index;
	        //     genIndexVector(vector_index, psd_est.size());
	        //     logVector("plot1", vector_index, psd_est);

            //     fsIN = (double)48e3;
            //     estimatePSD(audio_ds, NFFT_in , fsIN, freq, psd_est);
            //     genIndexVector(vector_index, psd_est.size());
	        //     logVector("plot2", vector_index, psd_est);

            //     fsIN = (double)48e3;
            //     estimatePSD(stereo_data_ds, NFFT_in , fsIN, freq, psd_est);
            //     genIndexVector(vector_index, psd_est.size());
	        //     logVector("plot3", vector_index, psd_est);
            // }
        
        
        }else{
            for(unsigned int k=0 ; k < audio_ds.size() ; k++){
                if(std::isnan(audio_ds[k])) audio_data[k] = 0;
                else audio_data[k] = (short int)(audio_ds[k] * 16384);
            }
            fwrite(&audio_data[0], sizeof(short int),audio_data.size(),stdout);
        }

        
		
}//Loop end

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