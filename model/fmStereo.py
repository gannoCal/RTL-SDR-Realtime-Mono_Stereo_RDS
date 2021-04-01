#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

"""
After you have installed the drivers to work with the RF dongle,
the 8-bit unsigned values for the I/Q pairs can be recorded as follows:

rtl_sdr -f 99.9M -s 2.4M - > my_samples_u8.raw

The above assumes that we are tuned to the FM station at 99.9 MHz,
we use an RF sample rate of 2.4 Msamples/sec and our file is called
my_samples_u8.raw (change as you see fit).

For the above use case, the data acquisition runs indefinitely,
hence the recording needs to be stopped by pressing Ctrl+C.
If we wish to stop it after a pre-defined number of samples,
e.g., 12 million I/Q pairs, we can use an extra argument:

rtl_sdr -f 99.9M -s 2.4M -n 12000000 - > my_samples_u8.raw

To check if the raw I/Q data has been recorded properly, place the file
in the "data" sub-folder from your project repo and run this Python file
from the "model" sub-folder. It should produce both the .png image files
(of the PSD estimates) for a couple of blocks, as well as the .wav file.

In the source code below (line 66) you can observe where the
normalization of the 8-bit unsigned raw I/Q samples is done;
while the range (-1 to +1) is an optional choice done by many,
it is at the discretion of each group how to handle the 8-bit
unsigned I/Q samples after they are read from the stdin in C++.

"""

import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy import signal
import numpy as np
import math

# use "custom" fmDemodArctan
from fmSupportLib import fmDemodArctan,fmPll

rf_Fs = 2.4e6
rf_Fc = 100e3
rf_taps = 151
rf_decim = 10

audio_Fs = 48e3
audio_Fc = 16e3
audio_taps = 151
audio_decim = 5


def generateSin(Fs, interval, frequency = 7.0, amplitude = 5.0, phase = 0.0,start = 0):

    dt = 1.0/Fs                          # sampling period (increment in time)
    time = np.arange(start, interval, dt)    # time vector over interval

    # generate the sin signal
    x = amplitude*np.sin(2*math.pi*frequency*time+phase)

    return time, x

if __name__ == "__main__":

    # read the raw IQ data from the recorded file
    # IQ data is assumed to be in 8-bits unsigned (and interleaved)
    in_fname = "../data/test1.raw"
    raw_data = np.fromfile(in_fname, dtype='float32')
    # IQ data is normalized between -1 and +1
    #iq_data = (raw_data - 128.0)/128.0
    iq_data = raw_data
    print("Read raw RF data from \"" + in_fname + "\" in unsigned 969-bit format")

    # coefficients for the front-end low-pass filter
    rf_coeff = signal.firwin(rf_taps, \
        rf_Fc/(rf_Fs/2), window=('hann'))

    # coefficients for the filter to extract the mono audio
    audio_coeff = signal.firwin(audio_taps, \
        audio_Fc/((rf_Fs/rf_decim)/2), window=('hann'))


    indx = np.arange(0,len(rf_coeff),1,int)
    indxd_data = np.column_stack((indx,rf_coeff))
    
    np.savetxt('fmfilter.out', indxd_data, delimiter=',', fmt='%d  %1.5e')

    indx = np.arange(0,len(audio_coeff),1,int)
    indxd_data = np.column_stack((indx,audio_coeff))
    
    np.savetxt('audiofilter.out', indxd_data, delimiter=',', fmt='%d  %1.5e')

    # set up the better drawing
    fig1, (ax96) = plt.subplots(nrows=1)
    fig1.subplots_adjust(hspace = 1.0)

    # set up the drawing
    fig, (ax0, ax1, ax2) = plt.subplots(nrows=3)
    fig.subplots_adjust(hspace = 1.0)

    # select a block_size that is a multiple of KB
    # and a multiple of decimation factors
    block_size = 1024 * rf_decim * audio_decim * 2
    block_count = 0

    # states needed for continuity in block processing
    state_i_lpf_100k = np.zeros(rf_taps-1)
    state_q_lpf_100k = np.zeros(rf_taps-1)
    state_lpf_16k = np.zeros(audio_taps-1)
    state_lpf_st = np.zeros(audio_taps-1)
    state_tone = np.zeros(150)
    state_stereo = np.zeros(150)
    state_phase = 0

    # audio buffer that stores all the audio blocks
    audio_data = np.array([])
    audio_data_l = np.array([])
    audio_data_r = np.array([])
    audio_data_s = np.array([])

    audio_block_r = np.zeros(1024)
    audio_block_l = np.zeros(1024)
    audio_block_s = np.zeros(2048)

    state_integrator = 0.0
    state_phaseEst = 0.0
    state_feedbackI = 1.0
    state_feedbackQ = 0.0
    state_ncoOut0 = 1.0
    state_trigOffset = 0

    tone_pass_coeff = signal.firwin(151, [18.5e3/((rf_Fs/rf_decim)/2) , 19.5e3/((rf_Fs/rf_decim)/2)],window=('hann'),pass_zero=False)
    stereo_pass_coeff = signal.firwin(151, [22e3/((rf_Fs/rf_decim)/2) , 54e3/((rf_Fs/rf_decim)/2)],window=('hann'),pass_zero=False)

    # if the number of samples in the last block is less than the block size
    # it is fine to ignore the last few samples from the raw I/Q file
    while (block_count+1)*block_size < len(iq_data):
        # if you wish to have shorter runtimes while troubleshooting
        # you can control the above loop exit condition as you see fit

        print('Processing block ' + str(block_count))

        # filter to extract the FM channel (I samples are even, Q samples are odd)
        i_filt, state_i_lpf_100k = signal.lfilter(rf_coeff, 1.0, \
                iq_data[(block_count)*block_size:(block_count+1)*block_size:2],
                zi=state_i_lpf_100k)
        q_filt, state_q_lpf_100k = signal.lfilter(rf_coeff, 1.0, \
                iq_data[(block_count)*block_size+1:(block_count+1)*block_size:2],
                zi=state_q_lpf_100k)

        # downsample the I/Q data from the FM channel
        i_ds = i_filt[::rf_decim]
        q_ds = q_filt[::rf_decim]

        

        # FM demodulator
        fm_demod, state_phase = fmDemodArctan(i_ds, q_ds, state_phase)


        

        tone_filt, state_tone = signal.lfilter(tone_pass_coeff, 1.0, \
                fm_demod,
                zi=state_tone)

        stereo_filt, state_stereo = signal.lfilter(stereo_pass_coeff, 1.0, \
                fm_demod,
                zi=state_stereo)        


        ncoOut, state_integrator, state_phaseEst, state_feedbackI, state_feedbackQ, state_ncoOut0, state_trigOffset \
             = fmPll(tone_filt,19e3,rf_Fs/rf_decim,ncoScale = 2, \
        state_integrator = state_integrator, \
        state_phaseEst = state_phaseEst, \
        state_feedbackI = state_feedbackI, \
        state_feedbackQ = state_feedbackQ, \
        state_ncoOut0 = state_ncoOut0, \
        state_trigOffset = state_trigOffset)

        st_data = np.zeros(len(stereo_filt))
        
        for i in range(0 , len(st_data)):
            st_data[i] = ncoOut[i] * stereo_filt[i] * 2

        st_filt, state_lpf_st = signal.lfilter(audio_coeff, 1.0, \
                st_data, zi=state_lpf_st)


        # time1,x1 = generateSin(240e3,6*6e-5,19e3,1.0,0)
        # time2,x2 = generateSin(240e3,2*6*6e-5,19e3,1.0,np.pi/2,start=6*6e-5)
        # x = np.append(x1,x2)
        # newTime = np.array([(time2[3]-time2[2])+time2[len(time2)-1]])
        
        # time2 = np.append(time2,newTime)
        # time = np.append(time1,time2)

        # x, state_integrator, state_phaseEst, state_feedbackI, state_feedbackQ, state_ncoOut0, state_trigOffset \
        #      = fmPll(x,19e3,rf_Fs/rf_decim, \
        # state_integrator = state_integrator, \
        # state_phaseEst = state_phaseEst, \
        # state_feedbackI = state_feedbackI, \
        # state_feedbackQ = state_feedbackQ, \
        # state_ncoOut0 = state_ncoOut0, \
        # state_trigOffset = state_trigOffset)

        # extract the mono audio data through filtering
        audio_filt, state_lpf_16k = signal.lfilter(audio_coeff, 1.0, \
                fm_demod, zi=state_lpf_16k)

        # downsample audio data
        audio_block = audio_filt[::audio_decim]
        st_block = st_filt[::audio_decim]

        for i in range(0,len(audio_block)):
            audio_block_s[2*i+1] = (audio_block[i] - st_block[i]  )   /2     #r
            audio_block_s[2*i] = (audio_block[i] + st_block[i]  )    /2    #l


        # concatenate the most recently processed audio_block
        # to the previous blocks stored already in audio_data
        #audio_data = np.concatenate((audio_data, audio_block))
        #audio_data_r = np.concatenate((audio_data_r, audio_block_r))
        #audio_data_l = np.concatenate((audio_data_l, audio_block_l))
        audio_data_s = np.concatenate((audio_data_s, audio_block_s))

        if(block_count == 1):
            #time = np.arange(0,round(len(step1)),1,int)
            # ax96.clear()
            # ax96.plot(time, x)
            # ax96.set_ylabel('ax96.set_ylabel')
            # ax96.set_xlabel('ax96.set_xlabel')
            # ax96.set_title('ax96.set_title')

            # fig1.savefig("../data/mrunobars5" + ".png")

            indx = np.arange(0,tone_pass_coeff.size,1,int)
            indxd_data = np.column_stack((indx,tone_pass_coeff))
    
            np.savetxt('tone_pass_coeff.out', indxd_data, delimiter=',', fmt='%d  %1.5e')

            indx = np.arange(0,stereo_pass_coeff.size,1,int)
            indxd_data = np.column_stack((indx,stereo_pass_coeff))
    
            np.savetxt('stereo_pass_coeff.out', indxd_data, delimiter=',', fmt='%d  %1.5e')

            indx = np.arange(0,state_phase.size,1,int)
            indxd_data = np.column_stack((indx,state_phase))
    
            np.savetxt('state_phase.out', indxd_data, delimiter=',', fmt='%d  %1.5e')

            indx = np.arange(0,state_lpf_16k.size,1,int)
            indxd_data = np.column_stack((indx,state_lpf_16k))
    
            np.savetxt('state_lpf_16k.out', indxd_data, delimiter=',', fmt='%d  %1.5e')

            indx = np.arange(0,state_q_lpf_100k.size,1,int)
            indxd_data = np.column_stack((indx,state_q_lpf_100k))
    
            np.savetxt('state_q_lpf_100k.out', indxd_data, delimiter=',', fmt='%d  %1.5e')

            indx = np.arange(0,state_i_lpf_100k.size,1,int)
            indxd_data = np.column_stack((indx,state_i_lpf_100k))
    
            np.savetxt('state_i_lpf_100k.out', indxd_data, delimiter=',', fmt='%d  %1.5e')


            indx = np.arange(0,len(audio_data),1,int)
            indxd_data = np.column_stack((indx,audio_data))
    
            np.savetxt('audio_data.out', indxd_data, delimiter=',', fmt='%d  %1.5e')


            indx = np.arange(0,len(audio_block),1,int)
            indxd_data = np.column_stack((indx,audio_block))
    
            np.savetxt('audio_block.out', indxd_data, delimiter=',', fmt='%d  %1.5e')

            firstIQ = iq_data[0:block_size]

            indx = np.arange(0,len(firstIQ),1,int)
            indxd_data = np.column_stack((indx,firstIQ))
    
            np.savetxt('iq.out', indxd_data, delimiter=',', fmt='%d  %1.5e')

            isamples_c = iq_data[(block_count)*block_size:(block_count+1)*block_size:2]
            qsamples_c = iq_data[(block_count)*block_size+1:(block_count+1)*block_size:2]

            indx = np.arange(0,len(isamples_c),1,int)
            indxd_data = np.column_stack((indx,isamples_c))
    
            np.savetxt('isamples.out', indxd_data, delimiter=',', fmt='%d  %1.5e')

            indx = np.arange(0,len(qsamples_c),1,int)
            indxd_data = np.column_stack((indx,qsamples_c))
    
            np.savetxt('qsamples.out', indxd_data, delimiter=',', fmt='%d  %1.5e')

            indx = np.arange(0,len(i_ds),1,int)
            indxd_data = np.column_stack((indx,i_ds))
    
            np.savetxt('ids.out', indxd_data, delimiter=',', fmt='%d  %1.5e')

            indx = np.arange(0,len(q_ds),1,int)
            indxd_data = np.column_stack((indx,q_ds))
    
            np.savetxt('qds.out', indxd_data, delimiter=',', fmt='%d  %1.5e')

            indx = np.arange(0,len(fm_demod),1,int)
            indxd_data = np.column_stack((indx,fm_demod))
    
            np.savetxt('fm_demod.out', indxd_data, delimiter=',', fmt='%d  %1.5e')
            #break

        # to save runtime select the range of blocks to log data
        # this includes both saving binary files as well plotting PSD
        # below we assume we want to plot for graphs for blocks 10 and 11
        if block_count >= 10 and block_count < 12:
            # PSD after FM demodulation
            ax0.clear()
            ax0.psd(fm_demod, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
            ax0.set_ylabel('PSD (dB/Hz)')
            ax0.set_xlabel('Freq (kHz)')
            ax0.set_title('Demodulated FM (block ' + str(block_count) + ')')
            # output binary file name (where samples are written from Python)
            fm_demod_fname = "../data/fm_demod_" + str(block_count) + ".bin"
            # create binary file where each sample is a 32-bit float
            fm_demod.astype('float32').tofile(fm_demod_fname)

            # PSD after extracting mono audio
            ax1.clear()
            ax1.psd(tone_filt, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
            ax1.set_ylabel('PSD (dB/Hz)')
            ax1.set_xlabel('Freq (kHz)')
            ax1.set_title('Extracted tone')

            # PSD after decimating mono audio
            ax2.clear()
            ax2.psd(stereo_filt, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
            ax2.set_ylabel('PSD (dB/Hz)')
            ax2.set_xlabel('Freq (kHz)')
            ax2.set_title('Extracted stereo')

            # save figure to file
            fig.savefig("../data/fmAudioStereo" + str(block_count) + ".png")

        block_count += 1

    print('Finished processing the raw I/Q samples')

    reshaped_data = np.reshape(audio_data_s, (-1, 2))
    # write audio data to a .wav file (assumes audio_data samples are -1 to +1)
    wavfile.write("../data/fmAudio2eroster.wav", int(audio_Fs), np.int16((reshaped_data/2)*32767))

    # uncomment assuming you wish to show some plots
    # plt.show()
