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


///////////////////Semaphore-eqsue construction
// class semaphore {
//     std::mutex mutex1;
//     std::condition_variable cvar;
//     long count = 0;

// public:
//     void signal() {
//         std::lock_guard<decltype(mutex1)> lock(mutex1);
//         ++count;
//         cvar.notify_one();
//     }

//     void wait() {
//         std::unique_lock<decltype(mutex1)> lock(mutex1);
//         while(!count)
//             cvar.wait(lock);
//         --count;
//     }
// };

void upsampleIQSamples(std::vector<double> &,std::vector<double> &);

void floatFEthreadMethod_00(bool &floatMode, 
    const int &block_size, 
    std::vector<double> &iq_data,
    int &sample_counter,
    std::vector<double> &i_samples_block,
    std::vector<double> &q_samples_block,
    std::vector<double> &rf_coeff,
    int &mode,
    std::vector<double> &i_ds,
    std::vector<double> &q_ds,
    std::vector<double> &state_i_lpf_100k,
    std::vector<double> &state_q_lpf_100k,
    std::vector<double> &state_phase,
    std::vector<double> &fm_demod,
    std::queue<std::vector<double>> &my_queue,
    std::mutex &my_mutex,
    std::condition_variable &my_cvar,
    std::queue<std::vector<double>> &RDS_queue,
    std::mutex &RDS_mutex,
    std::condition_variable &RDS_cvar

);

void floatFEthreadMethod_01(bool &floatMode, 
    const int &block_size, 
    std::vector<double> &iq_data,
    int &sample_counter,
    std::vector<double> &i_samples_block,
    std::vector<double> &q_samples_block,
    std::vector<double> &rf_coeff,
    int &mode,
    std::vector<double> &i_ds,
    std::vector<double> &q_ds,
    std::vector<double> &state_i_lpf_100k,
    std::vector<double> &state_q_lpf_100k,
    std::vector<double> &state_phase,
    std::vector<double> &fm_demod,
    std::queue<std::vector<double>> &RDS_queue,
    std::condition_variable &RDS_cvar,
    std::mutex &RDS_mutex

);

void binFEthreadMethod_0(bool &floatMode, 
    const int &block_size, 
    std::vector<double> &iq_data,
    int &sample_counter,
    std::vector<double> &i_samples_block,
    std::vector<double> &q_samples_block,
    std::vector<double> &rf_coeff,
    int &mode,
    std::vector<double> &i_ds,
    std::vector<double> &q_ds,
    std::vector<double> &state_i_lpf_100k,
    std::vector<double> &state_q_lpf_100k,
    std::vector<double> &state_phase,
    std::vector<double> &fm_demod,
    std::queue<std::vector<double>> &my_queue,
    std::mutex &my_mutex,
    std::condition_variable &my_cvar

);

void monoAuDiOtHrEaDmEtHoD_0(
    int &decimator,
    std::vector<double> &audio_ds,
    std::vector<double> &fm_demod,
    std::vector<double> &audio_coeff,
    std::vector<double> &state_conv,
    std::vector<double> &st_ds,
    std::vector<double> &tone_ds,
    std::vector<double> &st_coeff,
    std::vector<double> &tone_coeff,
    std::vector<double> &state_st_240k,
    std::vector<double> &state_tone_240k,
    double &PLLfreq,
    double &PLLfs,
    double &PLLNCOscale,
    double &phaseAdjust,
    double &normBandwidth,
    std::vector<double> &ncoOut,
    std::vector<double> &prevstate,
    std::vector<double> &stereo_data_ds,
    std::vector<double> &stereo_data,
    std::vector<double> &state_stereo_data,
    std::vector<short int> &stereo_block,
    std::vector<short int> &audio_data,
    std::queue<std::vector<double>> &my_queue,
    std::mutex &my_mutex,
    std::condition_variable &my_cvar

);

void stereoAuDiOtHrEaDmEtHoD_0(
    int &decimator,
    std::vector<double> &audio_ds,
    std::vector<double> &fm_demod,
    std::vector<double> &audio_coeff,
    std::vector<double> &state_conv,
    std::vector<double> &st_ds,
    std::vector<double> &tone_ds,
    std::vector<double> &st_coeff,
    std::vector<double> &tone_coeff,
    std::vector<double> &state_st_240k,
    std::vector<double> &state_tone_240k,
    double &PLLfreq,
    double &PLLfs,
    double &PLLNCOscale,
    double &phaseAdjust,
    double &normBandwidth,
    std::vector<double> &ncoOut,
    std::vector<double> &prevstate,
    std::vector<double> &stereo_data_ds,
    std::vector<double> &stereo_data,
    std::vector<double> &state_stereo_data,
    std::vector<short int> &stereo_block,
    std::vector<short int> &audio_data,
    std::queue<std::vector<double>> &my_queue,
    std::mutex &my_mutex,
    std::condition_variable &my_cvar

);


void floatFEthreadMethod_1(bool &floatMode, 
    const int &block_size, 
    std::vector<double> &iq_data,
    int &sample_counter,
    std::vector<double> &i_samples_block,
    std::vector<double> &q_samples_block,
    std::vector<double> &rf_coeff,
    int &mode,
    std::vector<double> &i_ds,
    std::vector<double> &q_ds,
    std::vector<double> &state_i_lpf_100k,
    std::vector<double> &state_q_lpf_100k,
    std::vector<double> &state_phase,
    std::vector<double> &fm_demod,
    std::queue<std::vector<double>> &my_queue,
    std::mutex &my_mutex,
    std::condition_variable &my_cvar

);

void binFEthreadMethod_1(bool &floatMode, 
    const int &block_size, 
    std::vector<double> &iq_data,
    int &sample_counter,
    std::vector<double> &i_samples_block,
    std::vector<double> &q_samples_block,
    std::vector<double> &rf_coeff,
    int &mode,
    std::vector<double> &i_ds,
    std::vector<double> &q_ds,
    std::vector<double> &state_i_lpf_100k,
    std::vector<double> &state_q_lpf_100k,
    std::vector<double> &state_phase,
    std::vector<double> &fm_demod,
    std::queue<std::vector<double>> &my_queue,
    std::mutex &my_mutex,
    std::condition_variable &my_cvar

);

void monoAuDiOtHrEaDmEtHoD_1(
    int &decimator,
    std::vector<double> &audio_ds,
    std::vector<double> &fm_demod,
    std::vector<double> &audio_coeff,
    std::vector<double> &state_conv,
    std::vector<double> &st_ds,
    std::vector<double> &tone_ds,
    std::vector<double> &st_coeff,
    std::vector<double> &tone_coeff,
    std::vector<double> &state_st_240k,
    std::vector<double> &state_tone_240k,
    double &PLLfreq,
    double &PLLfs,
    double &PLLNCOscale,
    double &phaseAdjust,
    double &normBandwidth,
    std::vector<double> &ncoOut,
    std::vector<double> &prevstate,
    std::vector<double> &stereo_data_ds,
    std::vector<double> &stereo_data,
    std::vector<double> &state_stereo_data,
    std::vector<short int> &stereo_block,
    std::vector<short int> &audio_data,
    std::queue<std::vector<double>> &my_queue,
    std::mutex &my_mutex,
    std::condition_variable &my_cvar

);

void stereoAuDiOtHrEaDmEtHoD_1(
    int &decimator,
    std::vector<double> &audio_ds,
    std::vector<double> &fm_demod,
    std::vector<double> &audio_coeff,
    std::vector<double> &state_conv,
    std::vector<double> &st_ds,
    std::vector<double> &tone_ds,
    std::vector<double> &st_coeff,
    std::vector<double> &tone_coeff,
    std::vector<double> &state_st_240k,
    std::vector<double> &state_tone_240k,
    double &PLLfreq,
    double &PLLfs,
    double &PLLNCOscale,
    double &phaseAdjust,
    double &normBandwidth,
    std::vector<double> &ncoOut,
    std::vector<double> &prevstate,
    std::vector<double> &stereo_data_ds,
    std::vector<double> &stereo_data,
    std::vector<double> &state_stereo_data,
    std::vector<short int> &stereo_block,
    std::vector<short int> &audio_data,
    std::queue<std::vector<double>> &my_queue,
    std::mutex &my_mutex,
    std::condition_variable &my_cvar

);



void RDS_0(
    const int &block_size,
    int &if_Fs,
    const int &rds_exreco_taps,
    std::queue<std::vector<double>> &my_queue,
    std::mutex &my_mutex,
    std::condition_variable &my_cvar
);




int main(int argc,char* argv[])
{
    bool floatMode = true;
    bool stereoMode = true;
	int mode = 0;
    if(argc < 2){
        //std::cerr << "Default - Mode 0";
    
	}else if(argc == 2){
		mode = atoi(argv[1]);
        if(mode != 1 && mode != 0){
            //std::cerr << "Wrong Mode! Exiting...";
            exit(1);
        }
        //std::cerr << "Mode is : " << mode << " \n";
	}else if(argc == 3){
        mode = atoi(argv[1]);
        char* ModeIn = argv[2];
        if(mode != 1 && mode != 0){
            //std::cerr << "Wrong Mode! Exiting...";
            exit(1);
        }

        int len = sizeof(ModeIn)/sizeof(*ModeIn);

        if(len < 2){
            //std::cerr << "Bad Param! Exiting...";
            exit(1);
        }

        if(ModeIn[0] == '-' && ModeIn[1]== 'f'){
            floatMode = true;
        }else if(ModeIn[0] == '-' && ModeIn[1]== 's'){
            stereoMode = true;
        }else{
            //std::cerr << "Bad Param! Exiting...";
            exit(1);
        }
        //std::cerr << "Mode is : " << mode << " \n";
        //std::cerr << "Param is : " << ModeIn[0] << ModeIn[1] << " \n";

    }else if(argc == 4){
        mode = atoi(argv[1]);
        char* ModeIn1 = argv[2];
        char* ModeIn2 = argv[3];
        if(mode != 1 && mode != 0){
            //std::cerr << "Wrong Mode! Exiting...";
            exit(1);
        }
        int len1 = sizeof(ModeIn1)/sizeof(*ModeIn1);
        int len2 = sizeof(ModeIn2)/sizeof(*ModeIn2);

        if(len1 < 2 || len2 < 2){
            //std::cerr << "Bad Param! Exiting...";
            exit(1);
        }

        if( ModeIn1[0] == ModeIn2[0] && ModeIn1[1] == ModeIn2[1]){
            //std::cerr << "Duplicate param! Exiting...";
            exit(1);
        }

        if(ModeIn1[0] == '-' && ModeIn1[1]== 'f'){
            floatMode = true;
        }else if(ModeIn1[0] == '-' && ModeIn1[1]== 's'){
            stereoMode = true;
        }else{
            //std::cerr << "Bad Param! Exiting...";
            exit(1);
        }

        if(ModeIn2[0] == '-' && ModeIn2[1]== 'f'){
            floatMode = true;
        }else if(ModeIn2[0] == '-' && ModeIn2[1]== 's'){
            stereoMode = true;
        }else{
            //std::cerr << "Bad Param! Exiting...";
            exit(1);
        }


        

        //std::cerr << "Mode is : " << mode << " \n";
        //std::cerr << "Param1 is : " << ModeIn1[0] << ModeIn1[1] << " \n";
        //std::cerr << "Param2 is : " << ModeIn2[0] << ModeIn2[1] << " \n";

    }else{
        //std::cerr << "Error. Please fix your input. " << " \n";
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
	const int audio_taps = Taps;

    int if_Fs = 240e3;
	const int sbanddecim = 1;
	const int sband_Fb = 22e3;
    const int sband_Fe = 54e3;
	const int sband_taps = Taps;
    
	const int tbanddecim = 1;
	const int tband_Fb = 18.5e3;
    const int tband_Fe = 19.5e3;
	const int tband_taps = Taps;

    const int rds_exreco_taps = Taps;

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
	std::vector<double> fm_demod_fe;
    std::vector<double> fm_demod_fe_rds;
    std::vector<double> fm_demod_audio;
    
    


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

    std::vector<double> state_i_lpf_100k_rds(rf_taps-1, 0);
	std::vector<double> state_q_lpf_100k_rds(rf_taps-1, 0);
	//vector<float> i_filt(rf_taps-1, 0);
	//vector<float> q_filt(rf_taps-1, 0);
	std::vector<double> i_ds((block_size)/20, 0);
	std::vector<double> q_ds((block_size)/20, 0);

    std::vector<double> i_ds_rds((block_size)/20, 0);
	std::vector<double> q_ds_rds((block_size)/20, 0);

    std::vector<double> st_ds((block_size)/20, 0);
    std::vector<double> tone_ds((block_size)/20, 0);

	std::vector<double> audio_ds((block_size)/100, 0);
	std::vector<double> state_conv(audio_taps-1, 0);
	std::vector<double> state_phase(2, 0.0);

    std::vector<double> state_phase_rds(2, 0.0);
    
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
    double PLLfs = rf_Fs/rf_decim;
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
    int sample_counter_rds = 0;


    std::vector<short int> audio_data_rds(block_size/100);
    std::vector<double> iq_data_rds(block_size);
    std::vector<double> i_samples_block_rds, q_samples_block_rds;
    i_samples_block_rds.resize(iq_data.size()/2,(double)0.0);
    q_samples_block_rds.resize(iq_data.size()/2,(double)0.0);

    std::queue<std::vector<double>> my_queue;
    std::queue<std::vector<double>> RDS_queue;
    std::mutex my_mutex;
    std::condition_variable my_cvar;
    std::mutex RDS_mutex;
    std::condition_variable RDS_cvar;


    
    



    



// if(mode == 0){
//     if(floatMode){
//         if(stereoMode){
//         /////////////////////////////////////////////////////////Stereo_float_0/////////////////////////////////////////////////////////
        
        std::thread fe = std::thread(
            floatFEthreadMethod_00,
            std::ref(floatMode), 
            std::cref(block_size), 
            std::ref(iq_data),
            std::ref(sample_counter),
            std::ref(i_samples_block),
            std::ref(q_samples_block),
            std::ref(rf_coeff),
            std::ref(mode),
            std::ref(i_ds),
            std::ref(q_ds),
            std::ref(state_i_lpf_100k),
            std::ref(state_q_lpf_100k),
            std::ref(state_phase),
            std::ref(fm_demod_fe),
            std::ref(my_queue),
            std::ref(my_mutex),
            std::ref(my_cvar),
            std::ref(RDS_queue),
            std::ref(RDS_mutex),
            std::ref(RDS_cvar)

        );

        // std::thread rdsfe = std::thread(
        //     floatFEthreadMethod_00,
        //     std::ref(floatMode), 
        //     std::cref(block_size), 
        //     std::ref(iq_data_rds),
        //     std::ref(sample_counter_rds),
        //     std::ref(i_samples_block_rds),
        //     std::ref(q_samples_block_rds),
        //     std::ref(rf_coeff),
        //     std::ref(mode),
        //     std::ref(i_ds_rds),
        //     std::ref(q_ds_rds),
        //     std::ref(state_i_lpf_100k_rds),
        //     std::ref(state_q_lpf_100k_rds),
        //     std::ref(state_phase_rds),
        //     std::ref(fm_demod_fe_rds),
        //     std::ref(RDS_queue),
        //     std::ref(RDS_mutex),
        //     std::ref(RDS_cvar)

        // );

        std::thread rds0 = std::thread(    
            RDS_0,
            std::cref(block_size),
            std::ref(if_Fs),
            std::cref(rds_exreco_taps),
            std::ref(RDS_queue),
            std::ref(RDS_mutex),
            std::ref(RDS_cvar)
        );


        std::thread audio = std::thread(
            stereoAuDiOtHrEaDmEtHoD_0,
            std::ref(decimator),
            std::ref(audio_ds),
            std::ref(fm_demod_audio),
            std::ref(audio_coeff),
            std::ref(state_conv),
            std::ref(st_ds),
            std::ref(tone_ds),
            std::ref(st_coeff),
            std::ref(tone_coeff),
            std::ref(state_st_240k),
            std::ref(state_tone_240k),
            std::ref(PLLfreq),
            std::ref(PLLfs),
            std::ref(PLLNCOscale),
            std::ref(phaseAdjust),
            std::ref(normBandwidth),
            std::ref(ncoOut),
            std::ref(prevstate),
            std::ref(stereo_data_ds),
            std::ref(stereo_data),
            std::ref(state_stereo_data),
            std::ref(stereo_block),
            std::ref(audio_data),
            std::ref(my_queue),
            std::ref(my_mutex),
            std::ref(my_cvar)
        );

        fe.join();
        //rdsfe.join();
        rds0.join();
        audio.join();
//         /////////////////////////////////////////////////////////Stereo_float_0/////////////////////////////////////////////////////////
//         }else{
//         /////////////////////////////////////////////////////////Mono_float_0/////////////////////////////////////////////////////////
        
            
            
            
//             std::thread fe = std::thread(
//             floatFEthreadMethod_00,
//             std::ref(floatMode), 
//             std::cref(block_size), 
//             std::ref(iq_data),
//             std::ref(sample_counter),
//             std::ref(i_samples_block),
//             std::ref(q_samples_block),
//             std::ref(rf_coeff),
//             std::ref(mode),
//             std::ref(i_ds),
//             std::ref(q_ds),
//             std::ref(state_i_lpf_100k),
//             std::ref(state_q_lpf_100k),
//             std::ref(state_phase),
//             std::ref(fm_demod_fe),
//             std::ref(my_queue),
//             std::ref(my_mutex),
//             std::ref(my_cvar)
//         );


//             std::thread rds0 = std::thread(    
//             RDS_0,
//             std::cref(block_size),
//             std::ref(if_Fs),
//             std::cref(rds_exreco_taps),
//             std::ref(RDS_queue),
//             std::ref(RDS_mutex),
//             std::ref(RDS_cvar)
//         );


//         std::thread audio = std::thread(
//             monoAuDiOtHrEaDmEtHoD_0,
//             std::ref(decimator),
//             std::ref(audio_ds),
//             std::ref(fm_demod_audio),
//             std::ref(audio_coeff),
//             std::ref(state_conv),
//             std::ref(st_ds),
//             std::ref(tone_ds),
//             std::ref(st_coeff),
//             std::ref(tone_coeff),
//             std::ref(state_st_240k),
//             std::ref(state_tone_240k),
//             std::ref(PLLfreq),
//             std::ref(PLLfs),
//             std::ref(PLLNCOscale),
//             std::ref(phaseAdjust),
//             std::ref(normBandwidth),
//             std::ref(ncoOut),
//             std::ref(prevstate),
//             std::ref(stereo_data_ds),
//             std::ref(stereo_data),
//             std::ref(state_stereo_data),
//             std::ref(stereo_block),
//             std::ref(audio_data),
//             std::ref(my_queue),
//             std::ref(my_mutex),
//             std::ref(my_cvar)
//         );

//         fe.join();
//         rds0.join();
//         audio.join();

//         /////////////////////////////////////////////////////////Mono_float_0/////////////////////////////////////////////////////////
//         }
//     }else{
//         if(stereoMode){
//         /////////////////////////////////////////////////////////Stereo_bin_0/////////////////////////////////////////////////////////
        
        
        
//         std::thread fe = std::thread(
//             binFEthreadMethod_0,
//             std::ref(floatMode), 
//             std::cref(block_size), 
//             std::ref(iq_data),
//             std::ref(sample_counter),
//             std::ref(i_samples_block),
//             std::ref(q_samples_block),
//             std::ref(rf_coeff),
//             std::ref(mode),
//             std::ref(i_ds),
//             std::ref(q_ds),
//             std::ref(state_i_lpf_100k),
//             std::ref(state_q_lpf_100k),
//             std::ref(state_phase),
//             std::ref(fm_demod_fe),
//             std::ref(my_queue),
//             std::ref(my_mutex),
//             std::ref(my_cvar)
//         );


//         std::thread rds0 = std::thread(    
//             RDS_0,
//             std::cref(block_size),
//             std::ref(if_Fs),
//             std::cref(rds_exreco_taps),
//             std::ref(my_queue),
//             std::ref(my_mutex),
//             std::ref(my_cvar)
//         );


//         std::thread audio = std::thread(
//             stereoAuDiOtHrEaDmEtHoD_0,
//             std::ref(decimator),
//             std::ref(audio_ds),
//             std::ref(fm_demod_audio),
//             std::ref(audio_coeff),
//             std::ref(state_conv),
//             std::ref(st_ds),
//             std::ref(tone_ds),
//             std::ref(st_coeff),
//             std::ref(tone_coeff),
//             std::ref(state_st_240k),
//             std::ref(state_tone_240k),
//             std::ref(PLLfreq),
//             std::ref(PLLfs),
//             std::ref(PLLNCOscale),
//             std::ref(phaseAdjust),
//             std::ref(normBandwidth),
//             std::ref(ncoOut),
//             std::ref(prevstate),
//             std::ref(stereo_data_ds),
//             std::ref(stereo_data),
//             std::ref(state_stereo_data),
//             std::ref(stereo_block),
//             std::ref(audio_data),
//             std::ref(my_queue),
//             std::ref(my_mutex),
//             std::ref(my_cvar)
//         );

//         fe.join();
//         rds0.join();
//         audio.join();
//         /////////////////////////////////////////////////////////Stereo_bin_0/////////////////////////////////////////////////////////
//         }else{
//         /////////////////////////////////////////////////////////Mono_bin_0/////////////////////////////////////////////////////////
           
            
            
//             std::thread fe = std::thread(
//             binFEthreadMethod_0,
//             std::ref(floatMode), 
//             std::cref(block_size), 
//             std::ref(iq_data),
//             std::ref(sample_counter),
//             std::ref(i_samples_block),
//             std::ref(q_samples_block),
//             std::ref(rf_coeff),
//             std::ref(mode),
//             std::ref(i_ds),
//             std::ref(q_ds),
//             std::ref(state_i_lpf_100k),
//             std::ref(state_q_lpf_100k),
//             std::ref(state_phase),
//             std::ref(fm_demod_fe),
//             std::ref(my_queue),
//             std::ref(my_mutex),
//             std::ref(my_cvar)
//         );


//         std::thread rds0 = std::thread(    
//             RDS_0,
//             std::cref(block_size),
//             std::ref(if_Fs),
//             std::cref(rds_exreco_taps),
//             std::ref(my_queue),
//             std::ref(my_mutex),
//             std::ref(my_cvar)
//         ); 


//         std::thread audio = std::thread(
//             monoAuDiOtHrEaDmEtHoD_0,
//             std::ref(decimator),
//             std::ref(audio_ds),
//             std::ref(fm_demod_audio),
//             std::ref(audio_coeff),
//             std::ref(state_conv),
//             std::ref(st_ds),
//             std::ref(tone_ds),
//             std::ref(st_coeff),
//             std::ref(tone_coeff),
//             std::ref(state_st_240k),
//             std::ref(state_tone_240k),
//             std::ref(PLLfreq),
//             std::ref(PLLfs),
//             std::ref(PLLNCOscale),
//             std::ref(phaseAdjust),
//             std::ref(normBandwidth),
//             std::ref(ncoOut),
//             std::ref(prevstate),
//             std::ref(stereo_data_ds),
//             std::ref(stereo_data),
//             std::ref(state_stereo_data),
//             std::ref(stereo_block),
//             std::ref(audio_data),
//             std::ref(my_queue),
//             std::ref(my_mutex),
//             std::ref(my_cvar)
//         );

//         fe.join();
//         rds0.join();
//         audio.join();

//         /////////////////////////////////////////////////////////Mono_bin_0/////////////////////////////////////////////////////////
//         }


//     }
// }else{
// if(floatMode){
//         if(stereoMode){
//         /////////////////////////////////////////////////////////Stereo_float_1/////////////////////////////////////////////////////////
//         std::thread fe = std::thread(
//             floatFEthreadMethod_1,
//             std::ref(floatMode), 
//             std::cref(block_size), 
//             std::ref(iq_data),
//             std::ref(sample_counter),
//             std::ref(i_samples_block),
//             std::ref(q_samples_block),
//             std::ref(rf_coeff),
//             std::ref(mode),
//             std::ref(i_ds),
//             std::ref(q_ds),
//             std::ref(state_i_lpf_100k),
//             std::ref(state_q_lpf_100k),
//             std::ref(state_phase),
//             std::ref(fm_demod_fe),
//             std::ref(my_queue),
//             std::ref(my_mutex),
//             std::ref(my_cvar)
//         );


//         std::thread audio = std::thread(
//             stereoAuDiOtHrEaDmEtHoD_1,
//             std::ref(decimator),
//             std::ref(audio_ds),
//             std::ref(fm_demod_audio),
//             std::ref(audio_coeff),
//             std::ref(state_conv),
//             std::ref(st_ds),
//             std::ref(tone_ds),
//             std::ref(st_coeff),
//             std::ref(tone_coeff),
//             std::ref(state_st_240k),
//             std::ref(state_tone_240k),
//             std::ref(PLLfreq),
//             std::ref(PLLfs),
//             std::ref(PLLNCOscale),
//             std::ref(phaseAdjust),
//             std::ref(normBandwidth),
//             std::ref(ncoOut),
//             std::ref(prevstate),
//             std::ref(stereo_data_ds),
//             std::ref(stereo_data),
//             std::ref(state_stereo_data),
//             std::ref(stereo_block),
//             std::ref(audio_data),
//             std::ref(my_queue),
//             std::ref(my_mutex),
//             std::ref(my_cvar)
//         );

//         fe.join();
//         audio.join();
//         /////////////////////////////////////////////////////////Stereo_float_1/////////////////////////////////////////////////////////
//         }else{
//         /////////////////////////////////////////////////////////Mono_float_1/////////////////////////////////////////////////////////
//             std::thread fe = std::thread(
//             floatFEthreadMethod_1,
//             std::ref(floatMode), 
//             std::cref(block_size), 
//             std::ref(iq_data),
//             std::ref(sample_counter),
//             std::ref(i_samples_block),
//             std::ref(q_samples_block),
//             std::ref(rf_coeff),
//             std::ref(mode),
//             std::ref(i_ds),
//             std::ref(q_ds),
//             std::ref(state_i_lpf_100k),
//             std::ref(state_q_lpf_100k),
//             std::ref(state_phase),
//             std::ref(fm_demod_fe),
//             std::ref(my_queue),
//             std::ref(my_mutex),
//             std::ref(my_cvar)
//         );


//         std::thread audio = std::thread(
//             monoAuDiOtHrEaDmEtHoD_1,
//             std::ref(decimator),
//             std::ref(audio_ds),
//             std::ref(fm_demod_audio),
//             std::ref(audio_coeff),
//             std::ref(state_conv),
//             std::ref(st_ds),
//             std::ref(tone_ds),
//             std::ref(st_coeff),
//             std::ref(tone_coeff),
//             std::ref(state_st_240k),
//             std::ref(state_tone_240k),
//             std::ref(PLLfreq),
//             std::ref(PLLfs),
//             std::ref(PLLNCOscale),
//             std::ref(phaseAdjust),
//             std::ref(normBandwidth),
//             std::ref(ncoOut),
//             std::ref(prevstate),
//             std::ref(stereo_data_ds),
//             std::ref(stereo_data),
//             std::ref(state_stereo_data),
//             std::ref(stereo_block),
//             std::ref(audio_data),
//             std::ref(my_queue),
//             std::ref(my_mutex),
//             std::ref(my_cvar)
//         );

//         fe.join();
//         audio.join();

//         /////////////////////////////////////////////////////////Mono_float_1/////////////////////////////////////////////////////////
//         }
//     }else{
//         if(stereoMode){
//         /////////////////////////////////////////////////////////Stereo_bin_1/////////////////////////////////////////////////////////
//         std::thread fe = std::thread(
//             binFEthreadMethod_1,
//             std::ref(floatMode), 
//             std::cref(block_size), 
//             std::ref(iq_data),
//             std::ref(sample_counter),
//             std::ref(i_samples_block),
//             std::ref(q_samples_block),
//             std::ref(rf_coeff),
//             std::ref(mode),
//             std::ref(i_ds),
//             std::ref(q_ds),
//             std::ref(state_i_lpf_100k),
//             std::ref(state_q_lpf_100k),
//             std::ref(state_phase),
//             std::ref(fm_demod_fe),
//             std::ref(my_queue),
//             std::ref(my_mutex),
//             std::ref(my_cvar)
//         );


//         std::thread audio = std::thread(
//             stereoAuDiOtHrEaDmEtHoD_1,
//             std::ref(decimator),
//             std::ref(audio_ds),
//             std::ref(fm_demod_audio),
//             std::ref(audio_coeff),
//             std::ref(state_conv),
//             std::ref(st_ds),
//             std::ref(tone_ds),
//             std::ref(st_coeff),
//             std::ref(tone_coeff),
//             std::ref(state_st_240k),
//             std::ref(state_tone_240k),
//             std::ref(PLLfreq),
//             std::ref(PLLfs),
//             std::ref(PLLNCOscale),
//             std::ref(phaseAdjust),
//             std::ref(normBandwidth),
//             std::ref(ncoOut),
//             std::ref(prevstate),
//             std::ref(stereo_data_ds),
//             std::ref(stereo_data),
//             std::ref(state_stereo_data),
//             std::ref(stereo_block),
//             std::ref(audio_data),
//             std::ref(my_queue),
//             std::ref(my_mutex),
//             std::ref(my_cvar)
//         );

//         fe.join();
//         audio.join();
//         /////////////////////////////////////////////////////////Stereo_bin_1/////////////////////////////////////////////////////////
//         }else{
//         /////////////////////////////////////////////////////////Mono_bin_1/////////////////////////////////////////////////////////
//             std::thread fe = std::thread(
//             binFEthreadMethod_1,
//             std::ref(floatMode), 
//             std::cref(block_size), 
//             std::ref(iq_data),
//             std::ref(sample_counter),
//             std::ref(i_samples_block),
//             std::ref(q_samples_block),
//             std::ref(rf_coeff),
//             std::ref(mode),
//             std::ref(i_ds),
//             std::ref(q_ds),
//             std::ref(state_i_lpf_100k),
//             std::ref(state_q_lpf_100k),
//             std::ref(state_phase),
//             std::ref(fm_demod_fe),
//             std::ref(my_queue),
//             std::ref(my_mutex),
//             std::ref(my_cvar)
//         );


//         std::thread audio = std::thread(
//             monoAuDiOtHrEaDmEtHoD_1,
//             std::ref(decimator),
//             std::ref(audio_ds),
//             std::ref(fm_demod_audio),
//             std::ref(audio_coeff),
//             std::ref(state_conv),
//             std::ref(st_ds),
//             std::ref(tone_ds),
//             std::ref(st_coeff),
//             std::ref(tone_coeff),
//             std::ref(state_st_240k),
//             std::ref(state_tone_240k),
//             std::ref(PLLfreq),
//             std::ref(PLLfs),
//             std::ref(PLLNCOscale),
//             std::ref(phaseAdjust),
//             std::ref(normBandwidth),
//             std::ref(ncoOut),
//             std::ref(prevstate),
//             std::ref(stereo_data_ds),
//             std::ref(stereo_data),
//             std::ref(state_stereo_data),
//             std::ref(stereo_block),
//             std::ref(audio_data),
//             std::ref(my_queue),
//             std::ref(my_mutex),
//             std::ref(my_cvar)
//         );

//         fe.join();
//         audio.join();

//         /////////////////////////////////////////////////////////Mono_bin_1/////////////////////////////////////////////////////////
//         }


//     }






// }

        
		
		

        

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


//////////////////////////////////////////////////////////////////////////Mode 0 Methods///////////////////////////////////////////////////////////////////

void RDS_0(
    const int &block_size,
    int &if_Fs,
    const int &rds_exreco_taps,
    std::queue<std::vector<double>> &my_queue,
    std::mutex &my_mutex,
    std::condition_variable &my_cvar
){

    

    std::vector<double> fm_demod_RDS;

    std::vector<double> extr_coeff;
    std::vector<double> reco_coeff;

    impulseResponseBPF(if_Fs, 54e3, 60e3, rds_exreco_taps, extr_coeff,1);
    impulseResponseBPF(if_Fs, 113.5e3, 114.5e3, rds_exreco_taps, reco_coeff,1);

    std::vector<double> extr_filt;
    std::vector<double> state_extr(rds_exreco_taps-1,0.0);


    std::vector<double> sq_extr_filt;
    std::vector<double> reco_filt;
    std::vector<double> state_reco(rds_exreco_taps-1,0.0);

    // double I_adj = PI/6 - PI/2 - PI/18 - 27*PI/180;
    // double Q_adj = -PI/2 + PI/6 - PI/2 - PI/18 - 27*PI/180;

    double I_adj = PI/6 - PI/2 - PI/18;
    double Q_adj = -PI/2 + PI/6 - PI/2 - PI/18;
    

    std::vector<double> prevstate(6,0);
    prevstate[0] = 0.0;            //state_integrator = 0.0;
    prevstate[1] = 0.0;            //state_phaseEst = 0.0;
    prevstate[2] = 1.0;            //state_feedbackI = 1.0;
    prevstate[3] = 0.0;            //state_feedbackQ = 0.0;
    prevstate[4] = 1.0;            //state_ncoOut0 = 1.0;
    prevstate[5] = 0.0;  

    double PLLfreq = 114e3;
    double PLLfs = 240e3;
    double PLLNCOscale = 0.5;
    double phaseAdjust = Q_adj;
    double normBandwidth = 0.01;

    std::vector<double> prevstate1(6,0);
    prevstate1[0] = 0.0;            //state_integrator = 0.0;
    prevstate1[1] = 0.0;            //state_phaseEst = 0.0;
    prevstate1[2] = 1.0;            //state_feedbackI = 1.0;
    prevstate1[3] = 0.0;            //state_feedbackQ = 0.0;
    prevstate1[4] = 1.0;            //state_ncoOut0 = 1.0;
    prevstate1[5] = 0.0;            //state_trigOffset = 0.0;

    double PLLfreq1 = 114e3;
    double PLLfs1 = 240e3;
    double PLLNCOscale1 = 0.5;
    double phaseAdjust1 = I_adj;
    double normBandwidth1 = 0.01;

    std::vector<double> ncoOut_I((block_size)/20, 0);
    std::vector<double> ncoOut_Q((block_size)/20, 0);

    std::vector<double> RDS_I;
    std::vector<double> RDS_Q;

    std::vector<double> RDS_I_filt;
    std::vector<double> state_RDS_I(rds_exreco_taps-1,0.0);

    std::vector<double> RDS_Q_filt;
    std::vector<double> state_RDS_Q(rds_exreco_taps-1,0.0);

    std::vector<double> rds_3k_coeff;

    impulseResponseLPF(240e3, 3e3, rds_exreco_taps, rds_3k_coeff,1);

    std::vector<double> rrc_coeff;
    impulseResponseRootRaisedCosine(57e3,rds_exreco_taps,rrc_coeff);

    std::vector<double> RRC_I;
    std::vector<double> state_RRC_I(rds_exreco_taps - 1,0.0);

    std::vector<double> RRC_I_samples;

    std::vector<double> RRC_Q;
    std::vector<double> state_RRC_Q(rds_exreco_taps - 1,0.0);

    std::vector<double> RRC_Q_samples;

    int last_sign = 0;
    int zero_at = 0;

    int t_shift_state = 0;
    int sample_spacing = 24;

    int bin_50 = 0;

    int shift_state = 0;
    int saved_i = 0;

    int e_count_limit = 2;
    int HH_LL_detected;
    int error_counter;

    int bin_50_state = 0;
    double next_symbol = 0.0;

    std::vector<double> binary_data;


    

    int block_count = 0;


    ///////Other Thread

    std::vector< std::vector<int> > parityArray = { {1,0,0,0,0,0,0,0,0,0}, 
                        {0,1,0,0,0,0,0,0,0,0}, 
                        {0,0,1,0,0,0,0,0,0,0}, 
                        {0,0,0,1,0,0,0,0,0,0}, 
                        {0,0,0,0,1,0,0,0,0,0}, 
                        {0,0,0,0,0,1,0,0,0,0}, 
                        {0,0,0,0,0,0,1,0,0,0}, 
                        {0,0,0,0,0,0,0,1,0,0}, 
                        {0,0,0,0,0,0,0,0,1,0}, 
                        {0,0,0,0,0,0,0,0,0,1}, 
                        {1,0,1,1,0,1,1,1,0,0}, 
                        {0,1,0,1,1,0,1,1,1,0}, 
                        {0,0,1,0,1,1,0,1,1,1}, 
                        {1,0,1,0,0,0,0,1,1,1}, 
                        {1,1,1,0,0,1,1,1,1,1}, 
                        {1,1,0,0,0,1,0,0,1,1}, 
                        {1,1,0,1,0,1,0,1,0,1}, 
                        {1,1,0,1,1,1,0,1,1,0}, 
                        {0,1,1,0,1,1,1,0,1,1}, 
                        {1,0,0,0,0,0,0,0,0,1}, 
                        {1,1,1,1,0,1,1,1,0,0}, 
                        {0,1,1,1,1,0,1,1,1,0}, 
                        {0,0,1,1,1,1,0,1,1,1}, 
                        {1,0,1,0,1,0,0,1,1,1}, 
                        {1,1,1,0,0,0,1,1,1,1}, 
                        {1,1,0,0,0,1,1,0,1,1} }; 

    std::vector<int> manchester_binary;
    int manchester_binary_state;

    std::vector<int> manchester_binary_processing_array(26);

    int previous_match = 0; 
    int is_nSync = 0; 
    int nSync_Hit_counter = 0; 
    int nSync_Flop_counter = 0; 
    int allowed_nSync_Flops = 6; 
    std::vector<std::vector<int>> found_array;
    std::vector<int> syndrome(10,0);

//printRealVectormdint(parityArray);

while(true){
    //std::cerr << "RDS" << block_count  << "\n";
        
    //////////////////Queue/////////////////////////
    std::unique_lock<std::mutex> my_lock(my_mutex);
        if(my_queue.empty()){
            //std::cerr << "RDS Audio Waiting...";
            my_cvar.wait(my_lock);
        }

        fm_demod_RDS = my_queue.front();
        my_queue.pop();
    



    //////////////////Queue/////////////////////////
    



    //////////////////Compute/////////////////////////
        //std::this_thread::sleep_for (std::chrono::milliseconds(10));
        convolveFIR_N_dec(1,extr_filt,fm_demod_RDS,extr_coeff,state_extr);
        sq_extr_filt.resize(extr_filt.size());
        for(auto i = 0 ; i < extr_filt.size() ; i++){
            sq_extr_filt[i] = extr_filt[i] * extr_filt[i]; 
        }
        convolveFIR_N_dec(1,reco_filt,sq_extr_filt,reco_coeff,state_reco);

        fmPll_RDS(reco_filt,PLLfreq1,PLLfs1,
        PLLNCOscale1,phaseAdjust1,normBandwidth1,ncoOut_I,prevstate1);

        fmPll_RDS(reco_filt,PLLfreq,PLLfs,
        PLLNCOscale,phaseAdjust,normBandwidth,ncoOut_Q,prevstate);
        
        RDS_Q.resize(ncoOut_I.size());
        RDS_I.resize(ncoOut_I.size());
        for(auto i = 0 ; i< extr_filt.size() ; i++){
            RDS_I[i] = ncoOut_I[i] * extr_filt[i] * 2;
            RDS_Q[i] = ncoOut_Q[i] * extr_filt[i] * 2;
        }

        // convolve_UPSAMPLE_N_dec(decimator, 24, stereo_data_ds, stereo_data, audio_coeff,state_stereo_data);
        convolve_UPSAMPLE_N_dec(80, 19, RDS_I_filt, RDS_I, rds_3k_coeff,state_RDS_I);
        

        //convolveFIR_N_dec(80,RDS_I_filt,RDS_I,rds_3k_coeff,state_RDS_I);
        //convolveFIR_N_dec(80,RDS_Q_filt,RDS_Q,rds_3k_coeff,state_RDS_Q);


        convolveFIR_N_dec(1,RRC_I,RDS_I_filt,rrc_coeff,state_RRC_I);

        convolve_UPSAMPLE_N_dec(80, 19, RDS_Q_filt, RDS_Q, rds_3k_coeff,state_RDS_Q);
        
        convolveFIR_N_dec(1,RRC_Q,RDS_Q_filt,rrc_coeff,state_RRC_Q);
        
        RRC_I_samples.resize(trunc(RRC_I.size()/24.0)+1);
        RRC_Q_samples.resize(trunc(RRC_Q.size()/24.0)+1);

        //////Zero Finder
        if(block_count == 0 /*|| error_counter == e_count_limit*/){
            shift_state = 0;
            last_sign = 0;
            zero_at = 0;
            for( auto i = 0; i < (RRC_I.size()) ; i++){
                if( (RRC_I[i] > 0 && last_sign == -1) || (RRC_I[i] < 0 && last_sign == 1)){
                    if(abs(RRC_I[i]) < abs(RRC_I[i-1]))
                        zero_at = i;
                    else
                        zero_at = i-1;
                    break;
                }
                if(RRC_I[i] < 0)
                    last_sign = -1;
                else
                    last_sign = 1;
            }
            t_shift_state = (zero_at + 12)%24;
            std::cerr <<t_shift_state <<" ss \n";
        }else{
            t_shift_state = 0;
        }

        //t_shift_state = 14;    
        ///////Sampler
        bin_50 = 0;
        saved_i = 0;
        for(auto i = 0 ; i < RRC_I_samples.size() ; i++){
            if(shift_state+sample_spacing*i + t_shift_state > RRC_I.size()-1){
                //saved_i++;
                bin_50 = 1;
                break;
            }else{
                RRC_I_samples[i] = RRC_I[shift_state+sample_spacing*i + t_shift_state]*20;
                RRC_Q_samples[i] = RRC_Q[shift_state+sample_spacing*i + t_shift_state]*20;
            }
            saved_i++;
        }
        saved_i = 50;
        //std::cerr << "i: " <<saved_i <<" \n";
        shift_state = (shift_state+sample_spacing*(saved_i-bin_50) + t_shift_state + 24) - (RRC_I.size()-1);
        shift_state--;
        //std::cerr << "Saved State: " <<shift_state <<" \n";

        if(block_count == 1){
            //std::vector<float> vector_index;
            //genIndexVector(vector_index,RRC_I.size());
            logVector2("RRC_I_c",RRC_I_samples,RRC_Q_samples );

            std::vector<float> vector_index;
            genIndexVector(vector_index,RRC_I.size());
            logVector("RRC_I",vector_index,RRC_I );

            //std::vector<float> vector_index;
            genIndexVector(vector_index,RRC_Q.size());
            logVector("RRC_Q",vector_index,RRC_Q );

            //std::vector<float> vector_index;
            genIndexVector(vector_index,rrc_coeff.size());
            logVector("RRC_coeff",vector_index,rrc_coeff );
        }


        ////////Symbol Analysis
        HH_LL_detected = 0;
        error_counter = 0;
        binary_data.clear();
        
        if(bin_50_state == 0){
            if(bin_50 == 0){
                binary_data.resize(25,0.0);

                for(auto i = 0; i <binary_data.size(); i++){
                    if(RRC_I_samples[2*i+1] > 0 && RRC_I_samples[2*i] < 0)       //LH
                        binary_data[i] = 0;
                    else if(RRC_I_samples[2*i+1] < 0 && RRC_I_samples[2*i] > 0)        //HL
                        binary_data[i] = 1;
                    else{
                        error_counter = error_counter + 1;
                        if(error_counter == e_count_limit){
                            HH_LL_detected = 1;
                            bin_50_state = 0;
                            next_symbol = RRC_I_samples[50];
                            break;
                        }else{
                            if(abs(RRC_I_samples[2*i+1]) > abs(RRC_I_samples[2*i])){
                                RRC_I_samples[2*i] = -1 * RRC_I_samples[2*i];
                                if(RRC_I_samples[2*i+1] > 0 && RRC_I_samples[2*i] < 0)        //LH
                                    binary_data[i] = 0;
                                else if(RRC_I_samples[2*i+1] < 0 && RRC_I_samples[2*i] > 0)        //HL
                                    binary_data[i] = 1;
                            }else{
                                RRC_I_samples[2*i+1] = -1 * RRC_I_samples[2*i+1];
                                if(RRC_I_samples[2*i+1] > 0 && RRC_I_samples[2*i] < 0)        //LH
                                    binary_data[i] = 0;
                                else if(RRC_I_samples[2*i+1] < 0 && RRC_I_samples[2*i] > 0)       //HL
                                    binary_data[i] = 1;
                            }
                        }
                    }
                }
                if(HH_LL_detected == 0){
                    bin_50_state = 1;
                    next_symbol = RRC_I_samples[50];
                }
            }else if(bin_50 == 1){
                binary_data.resize(25,0.0);

                for(auto i = 0; i <binary_data.size(); i++){
                    if(RRC_I_samples[2*i+1] > 0 && RRC_I_samples[2*i] < 0)        //LH
                        binary_data[i] = 0;

                    else if(RRC_I_samples[2*i+1] < 0 && RRC_I_samples[2*i] > 0)        //HL
                        binary_data[i] = 1;
                    else{
                        error_counter = error_counter + 1;
                        if(error_counter == e_count_limit){
                            HH_LL_detected = 1;
                            bin_50_state = 1;
                            next_symbol = RRC_I_samples[49];
                            break;
                        }else{
                            if(abs(RRC_I_samples[2*i+1]) > abs(RRC_I_samples[2*i])){
                                RRC_I_samples[2*i] = -1 * RRC_I_samples[2*i];
                                if(RRC_I_samples[2*i+1] > 0 && RRC_I_samples[2*i] < 0)        //LH
                                    binary_data[i] = 0;
                                else if(RRC_I_samples[2*i+1] < 0 && RRC_I_samples[2*i] > 0)        //HL
                                    binary_data[i] = 1;
                            }else{
                                RRC_I_samples[2*i+1] = -1 * RRC_I_samples[2*i+1];
                                if(RRC_I_samples[2*i+1] > 0 && RRC_I_samples[2*i] < 0)        //LH
                                    binary_data[i] = 0;
                                else if(RRC_I_samples[2*i+1] < 0 && RRC_I_samples[2*i] > 0)        //HL
                                    binary_data[i] = 1;
                            }
                        }
                    }
                }
                if(HH_LL_detected == 0){
                    bin_50_state = 0;
                    next_symbol = RRC_I_samples[50];
                }
            }
        }else if(bin_50_state == 1){
            if(bin_50 == 0){
                binary_data.resize(26,0.0);

                if(RRC_I_samples[0] > 0 && next_symbol < 0)       //LH
                    binary_data[0] = 0;

                else if(RRC_I_samples[0] < 0 && next_symbol > 0)       //HL
                    binary_data[0] = 1;
                else{
                    error_counter = error_counter + 1;
                    if(error_counter == e_count_limit){
                        HH_LL_detected = 1;
                        bin_50_state = 1;
                        next_symbol = RRC_I_samples[50];
                    }else{
                        if(abs(RRC_I_samples[0]) > abs(next_symbol)){
                            next_symbol = -1 * next_symbol;
                            if(RRC_I_samples[0] > 0 && next_symbol < 0)       //LH
                                binary_data[0] = 0;
                            else if(RRC_I_samples[0] < 0 && next_symbol > 0)       //HL
                                binary_data[0] = 1;
                        }else{
                            RRC_I_samples[0] = -1 * RRC_I_samples[0];
                            if(RRC_I_samples[0] > 0 && next_symbol < 0)        //LH
                                binary_data[0] = 0;
                            else if(RRC_I_samples[0] < 0 && next_symbol > 0)        //HL
                                binary_data[0] = 1;
                        }
                    }
                }

                if(HH_LL_detected == 0){
                    for(auto i = 1; i <binary_data.size(); i++){
                        if(RRC_I_samples[2*i] > 0 && RRC_I_samples[2*i-1] < 0)        //LH
                            binary_data[i] = 0;

                        else if(RRC_I_samples[2*i] < 0 && RRC_I_samples[2*i-1] > 0)        //HL
                            binary_data[i] = 1;
                        else{
                            error_counter = error_counter + 1;
                            if(error_counter == e_count_limit){
                                HH_LL_detected = 1;
                                bin_50_state = 1;
                                next_symbol = RRC_I_samples[50];
                                break;
                            }else{
                                if(abs(RRC_I_samples[2*i]) > abs(RRC_I_samples[2*i-1])){
                                    RRC_I_samples[2*i-1] = -1 * RRC_I_samples[2*i-1];
                                    if(RRC_I_samples[2*i] > 0 && RRC_I_samples[2*i-1] < 0)        //LH
                                        binary_data[i] = 0;
                                    else if(RRC_I_samples[2*i] < 0 && RRC_I_samples[2*i-1] > 0)        //HL
                                        binary_data[i] = 1;
                                }else{
                                    RRC_I_samples[2*i] = -1 * RRC_I_samples[2*i];
                                    if(RRC_I_samples[2*i] > 0 && RRC_I_samples[2*i-1] < 0)        //LH
                                        binary_data[i] = 0;
                                    else if(RRC_I_samples[2*i] < 0 && RRC_I_samples[2*i-1] > 0)        //HL
                                        binary_data[i] = 1;
                                }
                            }
                        }
                    }
                }
                if(HH_LL_detected == 0){
                    bin_50_state = 0;
                    next_symbol = RRC_I_samples[50];
                }
            }else if(bin_50 == 1){
                binary_data.resize(25,0.0);

                if(RRC_I_samples[0] > 0 && next_symbol < 0)        //LH
                    binary_data[0] = 0;

                else if(RRC_I_samples[0] < 0 && next_symbol > 0)        //HL
                    binary_data[0] = 1;
                else{
                    error_counter = error_counter + 1;
                    if(error_counter == e_count_limit){
                        HH_LL_detected = 1;
                        bin_50_state = 0;
                        next_symbol = RRC_I_samples[49];
                    }else{
                        if(abs(RRC_I_samples[0]) > abs(next_symbol)){
                            next_symbol = -1 * next_symbol;
                            if(RRC_I_samples[0] > 0 && next_symbol < 0)        //LH
                                binary_data[0] = 0;
                            else if(RRC_I_samples[0] < 0 && next_symbol > 0)        //HL
                                binary_data[0] = 1;
                        }else{
                            RRC_I_samples[0] = -1 * RRC_I_samples[0];
                            if(RRC_I_samples[0] > 0 && next_symbol < 0)        //LH
                                binary_data[0] = 0;
                            else if(RRC_I_samples[0] < 0 && next_symbol > 0)        //HL
                                binary_data[0] = 1;
                        }
                    }
                }
                if(HH_LL_detected == 0){
                    for(auto i = 1; i <binary_data.size(); i++){
                        if(RRC_I_samples[2*i] > 0 && RRC_I_samples[2*i-1] < 0)        //LH
                            binary_data[i] = 0;

                        else if(RRC_I_samples[2*i] < 0 && RRC_I_samples[2*i-1] > 0)        //HL
                            binary_data[i] = 1;
                        else{
                            error_counter = error_counter + 1;
                            if(error_counter == e_count_limit){
                                HH_LL_detected = 1;
                                bin_50_state = 0;
                                next_symbol = RRC_I_samples[49];
                                break;
                            }else{
                                if(abs(RRC_I_samples[2*i]) > abs(RRC_I_samples[2*i-1])){
                                    RRC_I_samples[2*i-1] = -1 * RRC_I_samples[2*i-1];
                                    if(RRC_I_samples[2*i] > 0 && RRC_I_samples[2*i-1] < 0)        //LH
                                        binary_data[i] = 0;
                                    else if(RRC_I_samples[2*i] < 0 && RRC_I_samples[2*i-1] > 0)        //HL
                                        binary_data[i] = 1;
                                }else{
                                    RRC_I_samples[2*i] = -1 * RRC_I_samples[2*i];
                                    if(RRC_I_samples[2*i] > 0 && RRC_I_samples[2*i-1] < 0)        //LH
                                        binary_data[i] = 0;
                                    else if(RRC_I_samples[2*i] < 0 && RRC_I_samples[2*i-1] > 0)        //HL
                                        binary_data[i] = 1;
                                }
                            }
                        }
                    }
                }
                if(HH_LL_detected == 0){
                    bin_50_state = 1;
                    next_symbol = RRC_I_samples[49];
                }
            }
        }

        //std::cerr << "HHLL " << ((error_counter == e_count_limit) ? "Detected" : "No!") <<" - Block: " << block_count << " \n";
        
        




        // /////////////New Thread
        manchester_binary.resize( binary_data.size(), 0);
        manchester_binary[0] = ((int)round(binary_data[0]) != manchester_binary_state) ? 1 : 0;

        process_MBA(manchester_binary[0],manchester_binary_processing_array,previous_match, is_nSync, 
        nSync_Hit_counter, nSync_Flop_counter, allowed_nSync_Flops, 
        found_array,
        parityArray, syndrome);
        // if(found_array.size() > 0)
        //     printRealVectormdint(found_array);
        for(auto i = 1 ; i < manchester_binary.size() ; i++){
            manchester_binary[i] = ((int)round(binary_data[i]) != (int)round(binary_data[i-1])) ? 1 : 0;
            process_MBA(manchester_binary[i],manchester_binary_processing_array,previous_match, is_nSync, 
            nSync_Hit_counter, nSync_Flop_counter, allowed_nSync_Flops, 
            found_array,
            parityArray, syndrome);
            // if(found_array.size() > 0)
            //     printRealVectormdint(found_array);
        }
        manchester_binary_state = binary_data[binary_data.size() - 1];
        








        block_count++;
        //////////////////Compute/////////////////////////


        //////////////////Queue/////////////////////////
        my_lock.unlock();
        my_cvar.notify_one();

        //////////////////Queue/////////////////////////


    }
}




void binFEthreadMethod_0(bool &floatMode, 
    const int &block_size, 
    std::vector<double> &iq_data,
    int &sample_counter,
    std::vector<double> &i_samples_block,
    std::vector<double> &q_samples_block,
    std::vector<double> &rf_coeff,
    int &mode,
    std::vector<double> &i_ds,
    std::vector<double> &q_ds,
    std::vector<double> &state_i_lpf_100k,
    std::vector<double> &state_q_lpf_100k,
    std::vector<double> &state_phase,
    std::vector<double> &fm_demod,
    std::queue<std::vector<double>> &my_queue,
    std::mutex &my_mutex,
    std::condition_variable &my_cvar


){
    while(true){


    //////////////////Generate/////////////////////////    
        
                readStdinBlockData(block_size,iq_data);
                

                if((std::cin.rdstate())!=0){
                    //std::cerr << "End of Input Stream";
                    exit(1);
                }

                

                sample_counter = 0;
                for (auto i = 0; i < iq_data.size() - 1; i = i + 2)
                {
                    i_samples_block[sample_counter] = iq_data[i];
                    q_samples_block[sample_counter] = iq_data[i+1];
                    sample_counter++;
                }
                // if(mode == 1) {
                // std::cout << "Mode is 1\n"; upsampleIQSamples(i_samples_block, q_samples_block);
                // }
                
                
                convolveFIR_N_dec(10, i_ds, i_samples_block, rf_coeff,state_i_lpf_100k);
                convolveFIR_N_dec(10, q_ds, q_samples_block, rf_coeff,state_q_lpf_100k);

                
                fmDemodArctanBlock(fm_demod,i_ds, q_ds, state_phase);
    //////////////////Generate/////////////////////////


    //////////////////Queue/////////////////////////
        std::unique_lock<std::mutex> my_lock(my_mutex);
        if(my_queue.size() == QUEUE_bLoCks){
            //std::cerr << "QUEUE Waiting...";
            my_cvar.wait(my_lock);
        }
        my_queue.push(fm_demod);
        

        my_lock.unlock();
        my_cvar.notify_one();







    //////////////////Queue/////////////////////////



    }
}

void floatFEthreadMethod_00(bool &floatMode, 
    const int &block_size, 
    std::vector<double> &iq_data,
    int &sample_counter,
    std::vector<double> &i_samples_block,
    std::vector<double> &q_samples_block,
    std::vector<double> &rf_coeff,
    int &mode,
    std::vector<double> &i_ds,
    std::vector<double> &q_ds,
    std::vector<double> &state_i_lpf_100k,
    std::vector<double> &state_q_lpf_100k,
    std::vector<double> &state_phase,
    std::vector<double> &fm_demod,
    std::queue<std::vector<double>> &my_queue,
    std::mutex &my_mutex,
    std::condition_variable &my_cvar,
    std::queue<std::vector<double>> &RDS_queue,
    std::mutex &RDS_mutex,
    std::condition_variable &RDS_cvar


){
    while(true){
        
        //std::cerr << "FE\n";
        //std::cerr << "Size AUDIO : "<<my_queue.size()<<" Size RDS : "<<RDS_queue.size()<< " -FE\n";

    //////////////////Generate/////////////////////////    
        
                readStdinBlockDataFloat(block_size,iq_data);

                if((std::cin.rdstate())!=0){
                    //std::cerr << "End of Input Stream";
                    exit(1);
                }

                

                sample_counter = 0;
                for (auto i = 0; i < iq_data.size() - 1; i = i + 2)
                {
                    i_samples_block[sample_counter] = iq_data[i];
                    q_samples_block[sample_counter] = iq_data[i+1];
                    sample_counter++;
                }
                // if(mode == 1) {
                // std::cout << "Mode is 1\n"; upsampleIQSamples(i_samples_block, q_samples_block);
                // }
                
                
                convolveFIR_N_dec(10, i_ds, i_samples_block, rf_coeff,state_i_lpf_100k);
                convolveFIR_N_dec(10, q_ds, q_samples_block, rf_coeff,state_q_lpf_100k);

                
                fmDemodArctanBlock(fm_demod,i_ds, q_ds, state_phase);
    //////////////////Generate/////////////////////////

        //std::cerr << "FE about to lock\n";
    //////////////////Queue/////////////////////////
        
        //std::cerr << "FE locked\n";
        std::unique_lock<std::mutex> my_lock(my_mutex);
        if(my_queue.size() == QUEUE_bLoCks){
                
                //std::cerr << "QUEUE Waiting...\n";
                my_cvar.wait(my_lock);
                
        }
        my_queue.push(fm_demod);

        std::unique_lock<std::mutex> RDS_lock(RDS_mutex);
        if(RDS_queue.size() == QUEUE_bLoCks){
                
                //std::cerr << "QUEUE Waiting...\n";
                RDS_cvar.wait(RDS_lock);
                
            
        }



        RDS_queue.push(fm_demod);
        //std::cerr << "FE Pushed\n";

        RDS_lock.unlock();

        
        RDS_cvar.notify_all();

        
        

        //std::cerr << "FE Pushed\n";

        my_lock.unlock();
        my_cvar.notify_all();

        
        







    //////////////////Queue/////////////////////////



    }
}

void floatFEthreadMethod_01(bool &floatMode, 
    const int &block_size, 
    std::vector<double> &iq_data,
    int &sample_counter,
    std::vector<double> &i_samples_block,
    std::vector<double> &q_samples_block,
    std::vector<double> &rf_coeff,
    int &mode,
    std::vector<double> &i_ds,
    std::vector<double> &q_ds,
    std::vector<double> &state_i_lpf_100k,
    std::vector<double> &state_q_lpf_100k,
    std::vector<double> &state_phase,
    std::vector<double> &fm_demod,
    std::queue<std::vector<double>> &RDS_queue,
    std::condition_variable &RDS_cvar,
    std::mutex &RDS_mutex


){
    while(true){
        
        //std::cerr << "FE\n";
        //std::cerr << "Size AUDIO : "<<my_queue.size()<<" Size RDS : "<< RDS_queue.size()<< " -FE\n";

    //////////////////Generate/////////////////////////    
        
                readStdinBlockDataFloat(block_size,iq_data);

                if((std::cin.rdstate())!=0){
                    //std::cerr << "End of Input Stream";
                    exit(1);
                }

                

                sample_counter = 0;
                for (auto i = 0; i < iq_data.size() - 1; i = i + 2)
                {
                    i_samples_block[sample_counter] = iq_data[i];
                    q_samples_block[sample_counter] = iq_data[i+1];
                    sample_counter++;
                }
                // if(mode == 1) {
                // std::cout << "Mode is 1\n"; upsampleIQSamples(i_samples_block, q_samples_block);
                // }
                
                
                convolveFIR_N_dec(10, i_ds, i_samples_block, rf_coeff,state_i_lpf_100k);
                convolveFIR_N_dec(10, q_ds, q_samples_block, rf_coeff,state_q_lpf_100k);

                
                fmDemodArctanBlock(fm_demod,i_ds, q_ds, state_phase);
    //////////////////Generate/////////////////////////

        //std::cerr << "FE about to lock\n";
    //////////////////Queue/////////////////////////
        
        //std::cerr << "FE locked\n";
        

        
        //std::cerr << "FE locked\n";
        std::unique_lock<std::mutex> RDS_lock(RDS_mutex);
        if(RDS_queue.size() == QUEUE_bLoCks){
                
                //std::cerr << "QUEUE Waiting...\n";
                RDS_cvar.wait(RDS_lock);
                
            
        }



        RDS_queue.push(fm_demod);
        //std::cerr << "FE Pushed\n";

        RDS_lock.unlock();

        
        RDS_cvar.notify_all();







    //////////////////Queue/////////////////////////



    }
}


void stereoAuDiOtHrEaDmEtHoD_0(
    int &decimator,
    std::vector<double> &audio_ds,
    std::vector<double> &fm_demod,
    std::vector<double> &audio_coeff,
    std::vector<double> &state_conv,
    std::vector<double> &st_ds,
    std::vector<double> &tone_ds,
    std::vector<double> &st_coeff,
    std::vector<double> &tone_coeff,
    std::vector<double> &state_st_240k,
    std::vector<double> &state_tone_240k,
    double &PLLfreq,
    double &PLLfs,
    double &PLLNCOscale,
    double &phaseAdjust,
    double &normBandwidth,
    std::vector<double> &ncoOut,
    std::vector<double> &prevstate,
    std::vector<double> &stereo_data_ds,
    std::vector<double> &stereo_data,
    std::vector<double> &state_stereo_data,
    std::vector<short int> &stereo_block,
    std::vector<short int> &audio_data,
    std::queue<std::vector<double>> &my_queue,
    std::mutex &my_mutex,
    std::condition_variable &my_cvar

){
    //int loopcnt = 0;

    

    while(true){

        //std::cerr << "Stereo" << loopcnt  << "\n";
        
    //////////////////Queue/////////////////////////
        std::unique_lock<std::mutex> my_lock(my_mutex);
        if(my_queue.empty()){
            //std::cerr << "St Audio Waiting...";
            my_cvar.wait(my_lock);
        }

        fm_demod = my_queue.front();
        my_queue.pop();




    //////////////////Queue/////////////////////////




    //////////////////Compute/////////////////////////
        convolveFIR_N_dec(decimator, audio_ds,fm_demod,audio_coeff,state_conv);


        
        convolveFIR_N_dec(1, st_ds, fm_demod, st_coeff,state_st_240k);
        convolveFIR_N_dec(1, tone_ds, fm_demod, tone_coeff,state_tone_240k);
        
        fmPll(tone_ds,PLLfreq,PLLfs,PLLNCOscale,phaseAdjust,normBandwidth,ncoOut,prevstate);

        // std::vector<float> index_vector;
        // genIndexVector(index_vector, ncoOut.size());

        // if(loopcnt == 0){
        //     logVector("pll_block_0",index_vector,ncoOut);
        //     logVector("st_ds_0",index_vector,st_ds);
        // }
        // if(loopcnt == 1){
        //     logVector("pll_block_1",index_vector,ncoOut);
        //     logVector("st_ds_1",index_vector,st_ds);
        // }
        // if(loopcnt == 2){
        //     logVector("pll_block_2",index_vector,ncoOut);
        //     logVector("st_ds_2",index_vector,st_ds);
        // }
        // loopcnt++;

        for(auto i = 0 ; i < st_ds.size();i++){
            stereo_data[i] = ncoOut[i] * st_ds[i] * 2;
        }

        convolveFIR_N_dec(decimator, stereo_data_ds, stereo_data, audio_coeff,state_stereo_data);

        for(unsigned int i=0 ; i < audio_ds.size() ; i++){
        
        stereo_block[2*i+1] = std::isnan(audio_ds[i]) || std::isnan(stereo_data_ds[i]) ? (short int)0 :  (short int)(16384*(audio_ds[i] - stereo_data_ds[i]  ) / 2)   ;     //r
        stereo_block[2*i] = std::isnan(audio_ds[i]) || std::isnan(stereo_data_ds[i]) ? (short int)0 : (short int)(16384*(audio_ds[i]  + stereo_data_ds[i]  ) / 2)    ;   //l
        //stereo_block[2*i+1] =  (short int)(16384*(audio_ds[i]));     //r
        //stereo_block[2*i] = (short int)(16384*(audio_ds[i]));   //l
        
        }
        
        fwrite(&stereo_block[0], sizeof(short int),stereo_block.size(),stdout);
        
    
    
    
        //////////////////Compute/////////////////////////


        //////////////////Queue/////////////////////////
        my_lock.unlock();
            my_cvar.notify_one();

        //////////////////Queue/////////////////////////

        //loopcnt++;
    }
}

void monoAuDiOtHrEaDmEtHoD_0(
    int &decimator,
    std::vector<double> &audio_ds,
    std::vector<double> &fm_demod,
    std::vector<double> &audio_coeff,
    std::vector<double> &state_conv,
    std::vector<double> &st_ds,
    std::vector<double> &tone_ds,
    std::vector<double> &st_coeff,
    std::vector<double> &tone_coeff,
    std::vector<double> &state_st_240k,
    std::vector<double> &state_tone_240k,
    double &PLLfreq,
    double &PLLfs,
    double &PLLNCOscale,
    double &phaseAdjust,
    double &normBandwidth,
    std::vector<double> &ncoOut,
    std::vector<double> &prevstate,
    std::vector<double> &stereo_data_ds,
    std::vector<double> &stereo_data,
    std::vector<double> &state_stereo_data,
    std::vector<short int> &stereo_block,
    std::vector<short int> &audio_data,
    std::queue<std::vector<double>> &my_queue,
    std::mutex &my_mutex,
    std::condition_variable &my_cvar

){
    while(true){
    //////////////////Queue/////////////////////////
        







    //////////////////Queue/////////////////////////




    //////////////////Compute/////////////////////////
        convolveFIR_N_dec(decimator, audio_ds,fm_demod,audio_coeff,state_conv);


            
        for(unsigned int k=0 ; k < audio_ds.size() ; k++){
            if(std::isnan(audio_ds[k])) audio_data[k] = 0;
            else audio_data[k] = (short int)(audio_ds[k] * 16384);
        }
        fwrite(&audio_data[0], sizeof(short int),audio_data.size(),stdout);
    
        //////////////////Compute/////////////////////////


        //////////////////Queue/////////////////////////
             
        

            






        //////////////////Queue/////////////////////////


    }
}



//////////////////////////////////////////////////////////////////////////////////Mode 1 Methods///////////////////////////////////////////////////////////////////////








void binFEthreadMethod_1(bool &floatMode, 
    const int &block_size, 
    std::vector<double> &iq_data,
    int &sample_counter,
    std::vector<double> &i_samples_block,
    std::vector<double> &q_samples_block,
    std::vector<double> &rf_coeff,
    int &mode,
    std::vector<double> &i_ds,
    std::vector<double> &q_ds,
    std::vector<double> &state_i_lpf_100k,
    std::vector<double> &state_q_lpf_100k,
    std::vector<double> &state_phase,
    std::vector<double> &fm_demod,
    std::queue<std::vector<double>> &my_queue,
    std::mutex &my_mutex,
    std::condition_variable &my_cvar


){
    while(true){


    //////////////////Generate/////////////////////////    
        
                readStdinBlockData(block_size,iq_data);
                

                if((std::cin.rdstate())!=0){
                    //std::cerr << "End of Input Stream";
                    exit(1);
                }

                

                sample_counter = 0;
                for (auto i = 0; i < iq_data.size() - 1; i = i + 2)
                {
                    i_samples_block[sample_counter] = iq_data[i];
                    q_samples_block[sample_counter] = iq_data[i+1];
                    sample_counter++;
                }
                // if(mode == 1) {
                // std::cout << "Mode is 1\n"; upsampleIQSamples(i_samples_block, q_samples_block);
                // }
                
                
                convolveFIR_N_dec(10, i_ds, i_samples_block, rf_coeff,state_i_lpf_100k);
                convolveFIR_N_dec(10, q_ds, q_samples_block, rf_coeff,state_q_lpf_100k);

                
                fmDemodArctanBlock(fm_demod,i_ds, q_ds, state_phase);
    //////////////////Generate/////////////////////////


    //////////////////Queue/////////////////////////
        std::unique_lock<std::mutex> my_lock(my_mutex);
        if(my_queue.size() == QUEUE_bLoCks){
            //std::cerr << "QUEUE Waiting...";
            my_cvar.wait(my_lock);
        }
        my_queue.push(fm_demod);


        my_lock.unlock();
        my_cvar.notify_one();







    //////////////////Queue/////////////////////////



    }
}

void floatFEthreadMethod_1(bool &floatMode, 
    const int &block_size, 
    std::vector<double> &iq_data,
    int &sample_counter,
    std::vector<double> &i_samples_block,
    std::vector<double> &q_samples_block,
    std::vector<double> &rf_coeff,
    int &mode,
    std::vector<double> &i_ds,
    std::vector<double> &q_ds,
    std::vector<double> &state_i_lpf_100k,
    std::vector<double> &state_q_lpf_100k,
    std::vector<double> &state_phase,
    std::vector<double> &fm_demod,
    std::queue<std::vector<double>> &my_queue,
    std::mutex &my_mutex,
    std::condition_variable &my_cvar


){
    while(true){


    //////////////////Generate/////////////////////////    
        
                readStdinBlockDataFloat(block_size,iq_data);

                if((std::cin.rdstate())!=0){
                    //std::cerr << "End of Input Stream";
                    exit(1);
                }

                

                sample_counter = 0;
                for (auto i = 0; i < iq_data.size() - 1; i = i + 2)
                {
                    i_samples_block[sample_counter] = iq_data[i];
                    q_samples_block[sample_counter] = iq_data[i+1];
                    sample_counter++;
                }
                // if(mode == 1) {
                // std::cout << "Mode is 1\n"; upsampleIQSamples(i_samples_block, q_samples_block);
                // }
                
                
                convolveFIR_N_dec(10, i_ds, i_samples_block, rf_coeff,state_i_lpf_100k);
                convolveFIR_N_dec(10, q_ds, q_samples_block, rf_coeff,state_q_lpf_100k);

                
                fmDemodArctanBlock(fm_demod,i_ds, q_ds, state_phase);
    //////////////////Generate/////////////////////////


    //////////////////Queue/////////////////////////
        std::unique_lock<std::mutex> my_lock(my_mutex);
        if(my_queue.size() == QUEUE_bLoCks){
            //std::cerr << "QUEUE Waiting...";
            my_cvar.wait(my_lock);
        }
        my_queue.push(fm_demod);


        my_lock.unlock();
        my_cvar.notify_one();







    //////////////////Queue/////////////////////////



    }
}


void stereoAuDiOtHrEaDmEtHoD_1(
    int &decimator,
    std::vector<double> &audio_ds,
    std::vector<double> &fm_demod,
    std::vector<double> &audio_coeff,
    std::vector<double> &state_conv,
    std::vector<double> &st_ds,
    std::vector<double> &tone_ds,
    std::vector<double> &st_coeff,
    std::vector<double> &tone_coeff,
    std::vector<double> &state_st_240k,
    std::vector<double> &state_tone_240k,
    double &PLLfreq,
    double &PLLfs,
    double &PLLNCOscale,
    double &phaseAdjust,
    double &normBandwidth,
    std::vector<double> &ncoOut,
    std::vector<double> &prevstate,
    std::vector<double> &stereo_data_ds,
    std::vector<double> &stereo_data,
    std::vector<double> &state_stereo_data,
    std::vector<short int> &stereo_block,
    std::vector<short int> &audio_data,
    std::queue<std::vector<double>> &my_queue,
    std::mutex &my_mutex,
    std::condition_variable &my_cvar

){
    while(true){
    //////////////////Queue/////////////////////////
        std::unique_lock<std::mutex> my_lock(my_mutex);
        if(my_queue.empty()){
            //std::cerr << "St 1 Audio Waiting...";
            my_cvar.wait(my_lock);
        }

        fm_demod = my_queue.front();
        my_queue.pop();







    //////////////////Queue/////////////////////////




    //////////////////Compute/////////////////////////
        convolve_UPSAMPLE_N_dec(decimator, 24, audio_ds,fm_demod,audio_coeff,state_conv);


        
        convolveFIR_N_dec(1, st_ds, fm_demod, st_coeff,state_st_240k);
        convolveFIR_N_dec(1, tone_ds, fm_demod, tone_coeff,state_tone_240k);
        
        fmPll(tone_ds,PLLfreq,PLLfs,PLLNCOscale,phaseAdjust,normBandwidth,ncoOut,prevstate);
        
        for(auto i = 0 ; i < st_ds.size();i++){
            stereo_data[i] = ncoOut[i] * st_ds[i] * 2;
        }

        convolve_UPSAMPLE_N_dec(decimator, 24, stereo_data_ds, stereo_data, audio_coeff,state_stereo_data);

        for(unsigned int i=0 ; i < audio_ds.size() ; i++){
        
        stereo_block[2*i+1] = std::isnan(audio_ds[i]) || std::isnan(stereo_data_ds[i]) ? (short int)0 :  (short int)(16384*(audio_ds[i] - stereo_data_ds[i]  ) / 2)   ;     //r
        stereo_block[2*i] = std::isnan(audio_ds[i]) || std::isnan(stereo_data_ds[i]) ? (short int)0 : (short int)(16384*(audio_ds[i]  + stereo_data_ds[i]  ) / 2)    ;   //l
        //stereo_block[2*i+1] =  (short int)(16384*(audio_ds[i]));     //r
        //stereo_block[2*i] = (short int)(16384*(audio_ds[i]));   //l
        
        }
        
        fwrite(&stereo_block[0], sizeof(short int),stereo_block.size(),stdout);
        
    
    
    
        //////////////////Compute/////////////////////////


        //////////////////Queue/////////////////////////
            my_lock.unlock();
            my_cvar.notify_one();







        //////////////////Queue/////////////////////////


    }
}

void monoAuDiOtHrEaDmEtHoD_1(
    int &decimator,
    std::vector<double> &audio_ds,
    std::vector<double> &fm_demod,
    std::vector<double> &audio_coeff,
    std::vector<double> &state_conv,
    std::vector<double> &st_ds,
    std::vector<double> &tone_ds,
    std::vector<double> &st_coeff,
    std::vector<double> &tone_coeff,
    std::vector<double> &state_st_240k,
    std::vector<double> &state_tone_240k,
    double &PLLfreq,
    double &PLLfs,
    double &PLLNCOscale,
    double &phaseAdjust,
    double &normBandwidth,
    std::vector<double> &ncoOut,
    std::vector<double> &prevstate,
    std::vector<double> &stereo_data_ds,
    std::vector<double> &stereo_data,
    std::vector<double> &state_stereo_data,
    std::vector<short int> &stereo_block,
    std::vector<short int> &audio_data,
    std::queue<std::vector<double>> &my_queue,
    std::mutex &my_mutex,
    std::condition_variable &my_cvar

){
    while(true){
    //////////////////Queue/////////////////////////
        std::unique_lock<std::mutex> my_lock(my_mutex);
        if(my_queue.empty()){
            //std::cerr << "Audio Waiting...";
            my_cvar.wait(my_lock);
        }

        fm_demod = my_queue.front();
        my_queue.pop();







    //////////////////Queue/////////////////////////




    //////////////////Compute/////////////////////////
        convolve_UPSAMPLE_N_dec(decimator,24, audio_ds,fm_demod,audio_coeff,state_conv);


            
        for(unsigned int k=0 ; k < audio_ds.size() ; k++){
            if(std::isnan(audio_ds[k])) audio_data[k] = 0;
            else audio_data[k] = (short int)(audio_ds[k] * 16384);
        }
        fwrite(&audio_data[0], sizeof(short int),audio_data.size(),stdout);
    
        //////////////////Compute/////////////////////////


        //////////////////Queue/////////////////////////
            my_lock.unlock();
            my_cvar.notify_one();







        //////////////////Queue/////////////////////////


    }
}




