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
import queue



# use "custom" fmDemodArctan
from fmSupportLib import fmDemodArctan,fmPll,impulseResponseRootRaisedCosine

rf_Fs = 2.4e6
rf_Fc = 100e3
rf_taps = 151
rf_decim = 10

audio_Fs = 48e3
audio_Fc = 16e3
audio_taps = 151
audio_decim = 5

found_array = []

is_nSync = False
nSync_Hit_counter = 0
nSync_Flop_counter = 0
allowed_nSync_Flops = 6
previous_match = "None"




parityArray = np.array([ [1,0,0,0,0,0,0,0,0,0], \
                        [0,1,0,0,0,0,0,0,0,0], \
                        [0,0,1,0,0,0,0,0,0,0], \
                        [0,0,0,1,0,0,0,0,0,0], \
                        [0,0,0,0,1,0,0,0,0,0], \
                        [0,0,0,0,0,1,0,0,0,0], \
                        [0,0,0,0,0,0,1,0,0,0], \
                        [0,0,0,0,0,0,0,1,0,0], \
                        [0,0,0,0,0,0,0,0,1,0], \
                        [0,0,0,0,0,0,0,0,0,1], \
                        [1,0,1,1,0,1,1,1,0,0], \
                        [0,1,0,1,1,0,1,1,1,0], \
                        [0,0,1,0,1,1,0,1,1,1], \
                        [1,0,1,0,0,0,0,1,1,1], \
                        [1,1,1,0,0,1,1,1,1,1], \
                        [1,1,0,0,0,1,0,0,1,1], \
                        [1,1,0,1,0,1,0,1,0,1], \
                        [1,1,0,1,1,1,0,1,1,0], \
                        [0,1,1,0,1,1,1,0,1,1], \
                        [1,0,0,0,0,0,0,0,0,1], \
                        [1,1,1,1,0,1,1,1,0,0], \
                        [0,1,1,1,1,0,1,1,1,0], \
                        [0,0,1,1,1,1,0,1,1,1], \
                        [1,0,1,0,1,0,0,1,1,1], \
                        [1,1,1,0,0,0,1,1,1,1], \
                        [1,1,0,0,0,1,1,0,1,1] ] \
                                                )


def generateSin(Fs, interval, frequency = 7.0, amplitude = 5.0, phase = 0.0,start = 0):

    dt = 1.0/Fs                          # sampling period (increment in time)
    time = np.arange(start, interval, dt)    # time vector over interval

    # generate the sin signal
    x = amplitude*np.sin(2*math.pi*frequency*time+phase)

    return time, x


def Gfield_mult(MBA, parity_matrix):
    syndrome = np.zeros(10)
    for i in range(10):
        working_column = bool(MBA[0]) and bool(parityArray[0][i])
        for j in range(1,26):
            working_column = (bool(MBA[j]) and bool(parityArray[j][i])) ^ working_column
        syndrome[i] = bool_to_int(working_column)
    return syndrome

def Gfield_mult_reverse(MBA, parity_matrix):
    syndrome = np.zeros(10)
    for i in range(10):
        working_column = bool(MBA[25]) and bool(parityArray[0][i])
        for j in range(1,26):
            working_column = (bool(MBA[25-j]) and bool(parityArray[j][i])) ^ working_column
        syndrome[i] = bool_to_int(working_column)
    return syndrome

def process_MBA(new_bit,MBA):
    global found_array
    global is_nSync
    global nSync_Hit_counter
    global nSync_Flop_counter
    global allowed_nSync_Flops
    global previous_match
    UpdatedMBA = np.insert(MBA,0,new_bit)
    UpdatedMBA = np.delete(UpdatedMBA, -1)
    #print(UpdatedMBA)

    if(not is_nSync):
        syndrome = Gfield_mult_reverse(UpdatedMBA, parityArray)
        #print(syndrome)
        previous_match = "None"
    
        i = 0
        if(\
            syndrome[i+0] == 1 and\
            syndrome[i+1] == 1 and\
            syndrome[i+2] == 1 and\
            syndrome[i+3] == 1 and\
            syndrome[i+4] == 0 and\
            syndrome[i+5] == 1 and\
            syndrome[i+6] == 1 and\
            syndrome[i+7] == 0 and\
            syndrome[i+8] == 0 and\
            syndrome[i+9] == 0):
            #if found, add 3-block to array
            found_array.append(["A",0])
            previous_match = "D"
            print("Found Code A")

        elif(\
            syndrome[i+0] == 1 and\
            syndrome[i+1] == 1 and\
            syndrome[i+2] == 1 and\
            syndrome[i+3] == 1 and\
            syndrome[i+4] == 0 and\
            syndrome[i+5] == 1 and\
            syndrome[i+6] == 0 and\
            syndrome[i+7] == 1 and\
            syndrome[i+8] == 0 and\
            syndrome[i+9] == 0):
            found_array.append(["B",0])
            previous_match = "A"
            print("Found Code B")

        elif(\
            syndrome[i+0] == 1 and\
            syndrome[i+1] == 0 and\
            syndrome[i+2] == 0 and\
            syndrome[i+3] == 1 and\
            syndrome[i+4] == 0 and\
            syndrome[i+5] == 1 and\
            syndrome[i+6] == 1 and\
            syndrome[i+7] == 1 and\
            syndrome[i+8] == 0 and\
            syndrome[i+9] == 0):
            found_array.append(["C",0])
            previous_match = "B"
            print("Found Code C")

        elif(\
            syndrome[i+0] == 1 and\
            syndrome[i+1] == 1 and\
            syndrome[i+2] == 1 and\
            syndrome[i+3] == 1 and\
            syndrome[i+4] == 0 and\
            syndrome[i+5] == 0 and\
            syndrome[i+6] == 1 and\
            syndrome[i+7] == 1 and\
            syndrome[i+8] == 0 and\
            syndrome[i+9] == 0):
            found_array.append(["C",0])
            previous_match = "B"
            print("Found Code C'")
        
        elif(\
            syndrome[i+0] == 1 and\
            syndrome[i+1] == 0 and\
            syndrome[i+2] == 0 and\
            syndrome[i+3] == 1 and\
            syndrome[i+4] == 0 and\
            syndrome[i+5] == 1 and\
            syndrome[i+6] == 1 and\
            syndrome[i+7] == 0 and\
            syndrome[i+8] == 0 and\
            syndrome[i+9] == 0):
            found_array.append(["D",0])
            previous_match = "C"
            print("Found Code D")

        # parse array, look for count of 26 that matches current found with (next from) intital
        del_array = []
        for i in range(len(found_array)):
            if(found_array[i][0] == previous_match and found_array[i][1] == 26):
                is_nSync = True
                #previous_match = found_array[len(found_array)-1][0]
                if(previous_match == "A"):
                    previous_match = "B"
                elif(previous_match == "B"):
                    previous_match = "C"
                elif(previous_match == "C"):
                    previous_match = "D"
                elif(previous_match == "D"):
                    previous_match = "A"
                found_array = []
                print("Synchronized...")
                break
            elif(found_array[i][0] != previous_match and found_array[i][1] >= 26):
                del_array.append(i)
            else:
                found_array[i][1] = found_array[i][1] + 1

        if(not is_nSync):    
            for i in range(len(del_array)):
                del found_array[del_array[i]]
    else:
        if(nSync_Hit_counter == 26):
            codeFound = False
            syndrome = Gfield_mult_reverse(UpdatedMBA, parityArray)
            #print(syndrome)
        
            i = 0
            if(\
                syndrome[i+0] == 1 and\
                syndrome[i+1] == 1 and\
                syndrome[i+2] == 1 and\
                syndrome[i+3] == 1 and\
                syndrome[i+4] == 0 and\
                syndrome[i+5] == 1 and\
                syndrome[i+6] == 1 and\
                syndrome[i+7] == 0 and\
                syndrome[i+8] == 0 and\
                syndrome[i+9] == 0):
                codeFound = True
                if(previous_match != "D"):
                    if(nSync_Flop_counter >= allowed_nSync_Flops):
                        is_nSync = False
                        nSync_Hit_counter = 0
                        nSync_Flop_counter = 0
                        print("De-synchronized...")
                    else:
                        nSync_Flop_counter = nSync_Flop_counter + 1
                        print("Expected Code D")
                else:
                    nSync_Flop_counter = 0
                    print("Found Code A")

                if(previous_match == "A"):
                    previous_match = "B"
                elif(previous_match == "B"):
                    previous_match = "C"
                elif(previous_match == "C"):
                    previous_match = "D"
                elif(previous_match == "D"):
                    previous_match = "A"


                    
                #previous_match = "A"
                

            elif(\
                syndrome[i+0] == 1 and\
                syndrome[i+1] == 1 and\
                syndrome[i+2] == 1 and\
                syndrome[i+3] == 1 and\
                syndrome[i+4] == 0 and\
                syndrome[i+5] == 1 and\
                syndrome[i+6] == 0 and\
                syndrome[i+7] == 1 and\
                syndrome[i+8] == 0 and\
                syndrome[i+9] == 0):
                codeFound = True
                if(previous_match != "A"):
                    if(nSync_Flop_counter >= allowed_nSync_Flops):
                        is_nSync = False
                        nSync_Hit_counter = 0
                        nSync_Flop_counter = 0
                        print("De-synchronized...")
                    else:
                        nSync_Flop_counter = nSync_Flop_counter + 1
                        print("Expected Code A")
                else:
                    nSync_Flop_counter = 0
                    print("Found Code B")

                if(previous_match == "A"):
                    previous_match = "B"
                elif(previous_match == "B"):
                    previous_match = "C"
                elif(previous_match == "C"):
                    previous_match = "D"
                elif(previous_match == "D"):
                    previous_match = "A"



                #previous_match = "B"
                

            elif(\
                syndrome[i+0] == 1 and\
                syndrome[i+1] == 0 and\
                syndrome[i+2] == 0 and\
                syndrome[i+3] == 1 and\
                syndrome[i+4] == 0 and\
                syndrome[i+5] == 1 and\
                syndrome[i+6] == 1 and\
                syndrome[i+7] == 1 and\
                syndrome[i+8] == 0 and\
                syndrome[i+9] == 0):
                codeFound = True
                if(previous_match != "B"):
                    if(nSync_Flop_counter >= allowed_nSync_Flops):
                        is_nSync = False
                        nSync_Hit_counter = 0
                        nSync_Flop_counter = 0
                        print("De-synchronized...")
                    else:
                        nSync_Flop_counter = nSync_Flop_counter + 1
                        print("Expected Code B")
                else:
                    nSync_Flop_counter = 0
                    print("Found Code C")

                if(previous_match == "A"):
                    previous_match = "B"
                elif(previous_match == "B"):
                    previous_match = "C"
                elif(previous_match == "C"):
                    previous_match = "D"
                elif(previous_match == "D"):
                    previous_match = "A"


                #previous_match = "C"
                

            elif(\
                syndrome[i+0] == 1 and\
                syndrome[i+1] == 1 and\
                syndrome[i+2] == 1 and\
                syndrome[i+3] == 1 and\
                syndrome[i+4] == 0 and\
                syndrome[i+5] == 0 and\
                syndrome[i+6] == 1 and\
                syndrome[i+7] == 1 and\
                syndrome[i+8] == 0 and\
                syndrome[i+9] == 0):
                codeFound = True
                if(previous_match != "B"):
                    if(nSync_Flop_counter >= allowed_nSync_Flops):
                        is_nSync = False
                        nSync_Hit_counter = 0
                        nSync_Flop_counter = 0
                        print("De-synchronized...")
                    else:
                        nSync_Flop_counter = nSync_Flop_counter + 1
                        print("Expected Code B")
                else:
                    nSync_Flop_counter = 0
                    print("Found Code C'")

                if(previous_match == "A"):
                    previous_match = "B"
                elif(previous_match == "B"):
                    previous_match = "C"
                elif(previous_match == "C"):
                    previous_match = "D"
                elif(previous_match == "D"):
                    previous_match = "A"


                #previous_match = "C"
                
            
            elif(\
                syndrome[i+0] == 1 and\
                syndrome[i+1] == 0 and\
                syndrome[i+2] == 0 and\
                syndrome[i+3] == 1 and\
                syndrome[i+4] == 0 and\
                syndrome[i+5] == 1 and\
                syndrome[i+6] == 1 and\
                syndrome[i+7] == 0 and\
                syndrome[i+8] == 0 and\
                syndrome[i+9] == 0):
                codeFound = True
                if(previous_match != "C"):
                    if(nSync_Flop_counter >= allowed_nSync_Flops):
                        nSync_Hit_counter = 0
                        nSync_Flop_counter = 0
                        is_nSync = False
                        print("De-synchronized...")
                    else:
                        nSync_Hit_counter = 0
                        nSync_Flop_counter = nSync_Flop_counter + 1
                        print("Expected Code C")
                else:
                    nSync_Flop_counter = 0
                    print("Found Code D")

                if(previous_match == "A"):
                    previous_match = "B"
                elif(previous_match == "B"):
                    previous_match = "C"
                elif(previous_match == "C"):
                    previous_match = "D"
                elif(previous_match == "D"):
                    previous_match = "A"


                #previous_match = "D"
                


            if(not codeFound):
                if(nSync_Flop_counter >= allowed_nSync_Flops):
                    is_nSync = False
                    nSync_Hit_counter = 0
                    nSync_Flop_counter = 0
                    print("De-synchronized...")
                else:
                    nSync_Hit_counter = 0
                    nSync_Flop_counter = nSync_Flop_counter + 1
                    
                    if(previous_match == "A"):
                        previous_match = "B"
                    elif(previous_match == "B"):
                        previous_match = "C"
                    elif(previous_match == "C"):
                        previous_match = "D"
                    elif(previous_match == "D"):
                        previous_match = "A"
                    print("Expected Code "+previous_match)
                    # print("Found: ")
                    # print(syndrome)
            nSync_Hit_counter = 0        




        #print(nSync_Hit_counter)
        nSync_Hit_counter = nSync_Hit_counter + 1
    #print(previous_match)
    return UpdatedMBA
    

            


def bool_to_int(b):
    return 1 if b else 0

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
    fig1, ((ax96, ax961), (ax962, ax963)) = plt.subplots(nrows=2,ncols=2)
    fig1.subplots_adjust(hspace = 1.0)

    # set up the drawing
    fig, (ax0, ax1, ax2) = plt.subplots(nrows=3)
    fig.subplots_adjust(hspace = 1.0)

    # set up the constellation drawing
    fig2, (ax123) = plt.subplots(nrows=1)
    fig.subplots_adjust(hspace = 1.0)

    # set up the constellation drawing
    iwave, (axi) = plt.subplots(nrows=1)
    fig.subplots_adjust(hspace = 1.0)

    # set up the constellation drawing
    qwave, (axq) = plt.subplots(nrows=1)
    fig.subplots_adjust(hspace = 1.0)

    # set up the constellation drawing
    spiwave, (spaxi) = plt.subplots(nrows=1)
    fig.subplots_adjust(hspace = 1.0)

    # set up the constellation drawing
    spqwave, (spaxq) = plt.subplots(nrows=1)
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

    state_extr = np.zeros(150)
    state_reco = np.zeros(150)

    state_RDS_I = np.zeros(150)
    state_RDS_Q = np.zeros(150)

    state_RRC_I = np.zeros(150)
    state_RRC_Q = np.zeros(150)

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


    state_integrator1 = 0.0
    state_phaseEst1 = 0.0
    state_feedbackI1 = 1.0
    state_feedbackQ1 = 0.0
    state_ncoOut01 = 1.0
    state_trigOffset1 = 0

    shift_state = 0

    tone_pass_coeff = signal.firwin(151, [18.5e3/((rf_Fs/rf_decim)/2) , 19.5e3/((rf_Fs/rf_decim)/2)],window=('hann'),pass_zero='bandpass')
    stereo_pass_coeff = signal.firwin(151, [22e3/((rf_Fs/rf_decim)/2) , 54e3/((rf_Fs/rf_decim)/2)],window=('hann'),pass_zero='bandpass')


    extr_coeff = signal.firwin(151, [54e3 , 60e3],window=('hann'),pass_zero='bandpass', fs=240e3)
    reco_coeff = signal.firwin(151, [113.5e3 , 114.5e3],window=('hann'),pass_zero='bandpass',fs=240e3)

    rds_3k_coeff = signal.firwin(151, 3e3,window=('hann'),fs=240e3*19,pass_zero='lowpass')

    rrc_coeff = impulseResponseRootRaisedCosine(57e3,151)


    next_symbol = 0
    bin_50_state = 0
    manchester_binary_state = 1

    manchester_binary_processing_array = np.zeros(26)

    HH_LL_detected = 0

    # if the number of samples in the last block is less than the block size
    # it is fine to ignore the last few samples from the raw I/Q file
    while (block_count+1)*block_size < len(iq_data):
        # if you wish to have shorter runtimes while troubleshooting
        # you can control the above loop exit condition as you see fit

        #print('Processing block ' + str(block_count))

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


        

        # tone_filt, state_tone = signal.lfilter(tone_pass_coeff, 1.0, \
        #         fm_demod,
        #         zi=state_tone)

        # stereo_filt, state_stereo = signal.lfilter(stereo_pass_coeff, 1.0, \
        #         fm_demod,
        #         zi=state_stereo)        


        # ncoOut, state_integrator, state_phaseEst, state_feedbackI, state_feedbackQ, state_ncoOut0, state_trigOffset \
        #      = fmPll(tone_filt,19e3,rf_Fs/rf_decim,ncoScale = 2, \
        # state_integrator = state_integrator, \
        # state_phaseEst = state_phaseEst, \
        # state_feedbackI = state_feedbackI, \
        # state_feedbackQ = state_feedbackQ, \
        # state_ncoOut0 = state_ncoOut0, \
        # state_trigOffset = state_trigOffset)

        # st_data = np.zeros(len(stereo_filt))
        
        # for i in range(0 , len(st_data)):
        #     st_data[i] = ncoOut[i] * stereo_filt[i] * 2

        # st_filt, state_lpf_st = signal.lfilter(audio_coeff, 1.0, \
        #         st_data, zi=state_lpf_st)


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
        # st_block = st_filt[::audio_decim]

        # for i in range(0,len(audio_block)):
        #     audio_block_s[2*i+1] = (audio_block[i] - st_block[i]  )   /2     #r
        #     audio_block_s[2*i] = (audio_block[i] + st_block[i]  )    /2    #l


        # concatenate the most recently processed audio_block
        # to the previous blocks stored already in audio_data
        #audio_data = np.concatenate((audio_data, audio_block))
        #audio_data_r = np.concatenate((audio_data_r, audio_block_r))
        #audio_data_l = np.concatenate((audio_data_l, audio_block_l))
        audio_data_s = np.concatenate((audio_data_s, audio_block_s))










        #RDS stuff
        extr_filt, state_extr = signal.lfilter(extr_coeff, 1.0, \
                fm_demod, zi=state_extr)

        sq_extr_filt = np.zeros(len(extr_filt))

        for i in range(len(sq_extr_filt)):
            sq_extr_filt[i] = extr_filt[i]*extr_filt[i]

        reco_filt, state_reco = signal.lfilter(reco_coeff, 1.0, \
                sq_extr_filt, zi=state_reco)


        

        I_adj = np.pi/6 - np.pi/2 - np.pi/18
        Q_adj = np.pi/6 - np.pi/2 - np.pi/18

        #I samples
        ncoOut_I, state_integrator1, state_phaseEst1, state_feedbackI1, state_feedbackQ1, state_ncoOut01, state_trigOffset1 \
             = fmPll(reco_filt,114e3,rf_Fs/rf_decim,ncoScale = 0.5, phaseAdjust=0.0 + I_adj, \
        state_integrator = state_integrator1, \
        state_phaseEst = state_phaseEst1, \
        state_feedbackI = state_feedbackI1, \
        state_feedbackQ = state_feedbackQ1, \
        state_ncoOut0 = state_ncoOut01, \
        state_trigOffset = state_trigOffset1)


        #Q samples
        ncoOut_Q, state_integrator, state_phaseEst, state_feedbackI, state_feedbackQ, state_ncoOut0, state_trigOffset \
             = fmPll(reco_filt,114e3,rf_Fs/rf_decim,ncoScale = 0.5, phaseAdjust=-np.pi/2.0 + Q_adj, \
        state_integrator = state_integrator, \
        state_phaseEst = state_phaseEst, \
        state_feedbackI = state_feedbackI, \
        state_feedbackQ = state_feedbackQ, \
        state_ncoOut0 = state_ncoOut0, \
        state_trigOffset = state_trigOffset)

        RDS_I = np.zeros(len(ncoOut_I)*19)
        RDS_Q = np.zeros(len(ncoOut_I)*19)

        for i in range(len(extr_filt)):
            RDS_I[i*19] = ncoOut_I[i] * extr_filt[i]
            RDS_Q[i*19] = ncoOut_Q[i] * extr_filt[i]


        
        

        RDS_I_filt, state_RDS_I = signal.lfilter(rds_3k_coeff, 1.0, \
                RDS_I, zi=state_RDS_I)

        RDS_Q_filt, state_RDS_Q = signal.lfilter(rds_3k_coeff, 1.0, \
                RDS_Q, zi=state_RDS_Q)

        RDS_I_filt_ds = RDS_I_filt[::80]
        RDS_Q_filt_ds = RDS_Q_filt[::80]


        RRC_I, state_RRC_I = signal.lfilter(rrc_coeff, 1.0, \
                RDS_I_filt_ds, zi=state_RRC_I)

        RRC_Q, state_RRC_Q = signal.lfilter(rrc_coeff, 1.0, \
                RDS_Q_filt_ds, zi=state_RRC_Q)

        RRC_I_samples = np.zeros(math.trunc(len(RRC_I)/24)+1)
        RRC_Q_samples = np.zeros(math.trunc(len(RRC_I)/24)+1)

        
        if(block_count == 0):
            last_sign = 0
            zero_at = 0
            for i in range(len(RRC_I)):
                if( (RRC_I[i] > 0 and last_sign == -1) or (RRC_I[i] < 0 and last_sign == 1)):
                    if(abs(RRC_I[i]) < abs(RRC_I[i-1])):
                        zero_at = i
                    else:
                        zero_at = i-1
                    break
                if(RRC_I[i] < 0):
                    last_sign = -1
                else:
                    last_sign = 1
        
        # zero_array = np.array([[zero_at],[RRC_I[zero_at]]])
        # print(zero_array)

            t_shift_state = (zero_at + 12)%24
        else:
            t_shift_state = 0
        # if(t_shift_state > 15):
        #     t_shift_state = 15
        
        sample_shift = 0
        sample_spacing = 24
        
        # t_shift_state = shift_state
        # if(shift_state % 8 != 0 and block_count != 0):
        #     t_shift_state = t_shift_state+4

        bin_50 = 0
        RRC_I_sparse = np.zeros(len(RRC_I))
        RRC_Q_sparse = np.zeros(len(RRC_I))
        for i in range(len(RRC_Q_samples)):
            if(shift_state+sample_spacing*i + t_shift_state > len(RRC_I)-1):
                bin_50 = 1
                break
            else:
                RRC_I_samples[i] = RRC_I[shift_state+sample_spacing*i + t_shift_state]*20
                RRC_Q_samples[i] = RRC_Q[shift_state+sample_spacing*i + t_shift_state]*20
                RRC_I_sparse[shift_state+sample_spacing*i + t_shift_state] = RRC_I[shift_state+sample_spacing*i + t_shift_state]*20
                RRC_Q_sparse[shift_state+sample_spacing*i + t_shift_state] = RRC_Q[shift_state+sample_spacing*i + t_shift_state]*20
        


        # print((len(RRC_I)-1))
        # print(sample_shift+sample_spacing*i - 1 + t_shift_state)
        #print(t_shift_state)
        #print("T_state^")
        # print(shift_state)
        # print('Next')

        shift_state = (shift_state+sample_spacing*(i-bin_50) + t_shift_state + 24) - (len(RRC_I)-1)

        # if(t_shift_state + shift_state < 0):
        #     shift_state
        
        
        #print(shift_state)
        
        
        
        
        

        # if(bin_50_state == 0 and bin_50 == 0):
        #     binary_data = np.zeros(26)

        #     if(RRC_I_samples[0] > 0 and next_symbol < 0):        #HL
        #         binary_data[0] = 1

        #     elif(RRC_I_samples[0] < 0 and next_symbol > 0):        #LH
        #         binary_data[0] = 0
        #     else:
        #         HH_LL_detected = 1

        #     for i in range(1,len(binary_data)):
        #         if(RRC_I_samples[2*i] > 0 and RRC_I_samples[2*i-1] < 0):        #HL
        #             binary_data[i] = 1

        #         elif(RRC_I_samples[2*i] < 0 and RRC_I_samples[2*i-1] > 0):        #LH
        #             binary_data[i] = 0
        #         else:
        #             HH_LL_detected = 1

        # elif(bin_50_state == 1 and bin_50 == 0):
        #     binary_data = np.zeros(25)

        #     for i in range(0,len(binary_data)):
        #         if(RRC_I_samples[2*i+1] > 0 and RRC_I_samples[2*i] < 0):        #HL
        #             binary_data[i] = 1

        #         elif(RRC_I_samples[2*i] < 0 and RRC_I_samples[2*i] > 0):        #LH
        #             binary_data[i] = 0
        #         else:
        #             HH_LL_detected = 1
        
        # elif(bin_50 == 1 and  bin_50_state == 1):
        #     binary_data = np.zeros(25)

        #     for i in range(0,len(binary_data)):
        #         if(RRC_I_samples[2*i+1] > 0 and RRC_I_samples[2*i] < 0):        #HL
        #             binary_data[i] = 1

        #         elif(RRC_I_samples[2*i] < 0 and RRC_I_samples[2*i] > 0):        #LH
        #             binary_data[i] = 0
        #         else:
        #             HH_LL_detected = 1

        # elif(bin_50 == 1  bin_50_state == 1):
        #     print("yup")

        HH_LL_detected = 0
        error_counter = 0
        e_count_limit = 2
        if(bin_50_state == 0):
            if(bin_50 == 0):
                binary_data = np.zeros(25)

                for i in range(0,len(binary_data)):
                    if(RRC_I_samples[2*i+1] > 0 and RRC_I_samples[2*i] < 0):        #LH
                        binary_data[i] = 0
                    elif(RRC_I_samples[2*i+1] < 0 and RRC_I_samples[2*i] > 0):        #HL
                        binary_data[i] = 1
                    else:
                        error_counter = error_counter + 1
                        if(error_counter == e_count_limit):
                            #print("HHLL - 1: "+str(RRC_I_samples[2*i])+", 2: "+str(RRC_I_samples[2*i+1]))
                            HH_LL_detected = 1
                            bin_50_state = 0
                            next_symbol = RRC_I_samples[50]
                            break
                        else:
                            if(abs(RRC_I_samples[2*i+1]) > abs(RRC_I_samples[2*i])):
                                RRC_I_samples[2*i] = -1 * RRC_I_samples[2*i]
                                if(RRC_I_samples[2*i+1] > 0 and RRC_I_samples[2*i] < 0):        #LH
                                    binary_data[i] = 0
                                elif(RRC_I_samples[2*i+1] < 0 and RRC_I_samples[2*i] > 0):        #HL
                                    binary_data[i] = 1
                            else:
                                RRC_I_samples[2*i+1] = -1 * RRC_I_samples[2*i+1]
                                if(RRC_I_samples[2*i+1] > 0 and RRC_I_samples[2*i] < 0):        #LH
                                    binary_data[i] = 0
                                elif(RRC_I_samples[2*i+1] < 0 and RRC_I_samples[2*i] > 0):        #HL
                                    binary_data[i] = 1
                if(HH_LL_detected == 0):
                    bin_50_state = 1
                    next_symbol = RRC_I_samples[50]
            elif(bin_50 == 1):
                binary_data = np.zeros(25)

                for i in range(0,len(binary_data)):
                    if(RRC_I_samples[2*i+1] > 0 and RRC_I_samples[2*i] < 0):        #LH
                        binary_data[i] = 0

                    elif(RRC_I_samples[2*i+1] < 0 and RRC_I_samples[2*i] > 0):        #HL
                        binary_data[i] = 1
                    else:
                        error_counter = error_counter + 1
                        if(error_counter == e_count_limit):
                            #print("HHLL - 1: "+str(RRC_I_samples[2*i])+", 2: "+str(RRC_I_samples[2*i+1]))
                            HH_LL_detected = 1
                            bin_50_state = 1
                            next_symbol = RRC_I_samples[49]
                            break
                        else:
                            if(abs(RRC_I_samples[2*i+1]) > abs(RRC_I_samples[2*i])):
                                RRC_I_samples[2*i] = -1 * RRC_I_samples[2*i]
                                if(RRC_I_samples[2*i+1] > 0 and RRC_I_samples[2*i] < 0):        #LH
                                    binary_data[i] = 0
                                elif(RRC_I_samples[2*i+1] < 0 and RRC_I_samples[2*i] > 0):        #HL
                                    binary_data[i] = 1
                            else:
                                RRC_I_samples[2*i+1] = -1 * RRC_I_samples[2*i+1]
                                if(RRC_I_samples[2*i+1] > 0 and RRC_I_samples[2*i] < 0):        #LH
                                    binary_data[i] = 0
                                elif(RRC_I_samples[2*i+1] < 0 and RRC_I_samples[2*i] > 0):        #HL
                                    binary_data[i] = 1
                if(HH_LL_detected == 0):
                    bin_50_state = 0
                    next_symbol = RRC_I_samples[50]
        elif(bin_50_state == 1):
            if(bin_50 == 0):
                binary_data = np.zeros(26)

                if(RRC_I_samples[0] > 0 and next_symbol < 0):        #LH
                    binary_data[0] = 0

                elif(RRC_I_samples[0] < 0 and next_symbol > 0):        #HL
                    binary_data[0] = 1
                else:
                    error_counter = error_counter + 1
                    if(error_counter == e_count_limit):
                        #print("HHLL - 1: "+str(next_symbol)+", 2: "+str(RRC_I_samples[0]))
                        HH_LL_detected = 1
                        bin_50_state = 1
                        next_symbol = RRC_I_samples[50]
                    else:
                        if(abs(RRC_I_samples[0]) > abs(next_symbol)):
                            next_symbol = -1 * next_symbol
                            if(RRC_I_samples[0] > 0 and next_symbol < 0):        #LH
                                binary_data[0] = 0
                            elif(RRC_I_samples[0] < 0 and next_symbol > 0):        #HL
                                binary_data[0] = 1
                        else:
                            RRC_I_samples[0] = -1 * RRC_I_samples[0]
                            if(RRC_I_samples[0] > 0 and next_symbol < 0):        #LH
                                binary_data[0] = 0
                            elif(RRC_I_samples[0] < 0 and next_symbol > 0):        #HL
                                binary_data[0] = 1

                if(HH_LL_detected == 0):
                    for i in range(1,len(binary_data)):
                        if(RRC_I_samples[2*i] > 0 and RRC_I_samples[2*i-1] < 0):        #LH
                            binary_data[i] = 0

                        elif(RRC_I_samples[2*i] < 0 and RRC_I_samples[2*i-1] > 0):        #HL
                            binary_data[i] = 1
                        else:
                            error_counter = error_counter + 1
                            if(error_counter == e_count_limit):
                                #print("HHLL - 1: "+str(RRC_I_samples[2*i-1])+", 2: "+str(RRC_I_samples[2*i]))
                                HH_LL_detected = 1
                                bin_50_state = 1
                                next_symbol = RRC_I_samples[50]
                                break
                            else:
                                if(abs(RRC_I_samples[2*i]) > abs(RRC_I_samples[2*i-1])):
                                    RRC_I_samples[2*i-1] = -1 * RRC_I_samples[2*i-1]
                                    if(RRC_I_samples[2*i] > 0 and RRC_I_samples[2*i-1] < 0):        #LH
                                        binary_data[i] = 0
                                    elif(RRC_I_samples[2*i] < 0 and RRC_I_samples[2*i-1] > 0):        #HL
                                        binary_data[i] = 1
                                else:
                                    RRC_I_samples[2*i] = -1 * RRC_I_samples[2*i]
                                    if(RRC_I_samples[2*i] > 0 and RRC_I_samples[2*i-1] < 0):        #LH
                                        binary_data[i] = 0
                                    elif(RRC_I_samples[2*i] < 0 and RRC_I_samples[2*i-1] > 0):        #HL
                                        binary_data[i] = 1
                if(HH_LL_detected == 0):
                    bin_50_state = 0
                    next_symbol = RRC_I_samples[50]
            elif(bin_50 == 1):
                binary_data = np.zeros(25)

                if(RRC_I_samples[0] > 0 and next_symbol < 0):        #LH
                    binary_data[0] = 0

                elif(RRC_I_samples[0] < 0 and next_symbol > 0):        #HL
                    binary_data[0] = 1
                else:
                    error_counter = error_counter + 1
                    if(error_counter == e_count_limit):
                        #print("HHLL - 1: "+str(next_symbol)+", 2: "+str(RRC_I_samples[0]))
                        HH_LL_detected = 1
                        bin_50_state = 0
                        next_symbol = RRC_I_samples[49]
                    else:
                        if(abs(RRC_I_samples[0]) > abs(next_symbol)):
                            next_symbol = -1 * next_symbol
                            if(RRC_I_samples[0] > 0 and next_symbol < 0):        #LH
                                binary_data[0] = 0
                            elif(RRC_I_samples[0] < 0 and next_symbol > 0):        #HL
                                binary_data[0] = 1
                        else:
                            RRC_I_samples[0] = -1 * RRC_I_samples[0]
                            if(RRC_I_samples[0] > 0 and next_symbol < 0):        #LH
                                binary_data[0] = 0
                            elif(RRC_I_samples[0] < 0 and next_symbol > 0):        #HL
                                binary_data[0] = 1
                if(HH_LL_detected == 0):
                    for i in range(1,len(binary_data)):
                        if(RRC_I_samples[2*i] > 0 and RRC_I_samples[2*i-1] < 0):        #LH
                            binary_data[i] = 0

                        elif(RRC_I_samples[2*i] < 0 and RRC_I_samples[2*i-1] > 0):        #HL
                            binary_data[i] = 1
                        else:
                            error_counter = error_counter + 1
                            if(error_counter == e_count_limit):
                                #print("HHLL - 1: "+str(RRC_I_samples[2*i-1])+", 2: "+str(RRC_I_samples[2*i]))
                                HH_LL_detected = 1
                                bin_50_state = 0
                                next_symbol = RRC_I_samples[49]
                                break
                            else:
                                if(abs(RRC_I_samples[2*i]) > abs(RRC_I_samples[2*i-1])):
                                    RRC_I_samples[2*i-1] = -1 * RRC_I_samples[2*i-1]
                                    if(RRC_I_samples[2*i] > 0 and RRC_I_samples[2*i-1] < 0):        #LH
                                        binary_data[i] = 0
                                    elif(RRC_I_samples[2*i] < 0 and RRC_I_samples[2*i-1] > 0):        #HL
                                        binary_data[i] = 1
                                else:
                                    RRC_I_samples[2*i] = -1 * RRC_I_samples[2*i]
                                    if(RRC_I_samples[2*i] > 0 and RRC_I_samples[2*i-1] < 0):        #LH
                                        binary_data[i] = 0
                                    elif(RRC_I_samples[2*i] < 0 and RRC_I_samples[2*i-1] > 0):        #HL
                                        binary_data[i] = 1
                if(HH_LL_detected == 0):
                    bin_50_state = 1
                    next_symbol = RRC_I_samples[49]
        


          
        manchester_binary = np.zeros(len(binary_data))
        manchester_binary[0] = bool_to_int(bool(binary_data[0] == 1) ^ bool(manchester_binary_state == 1) )
        #print(manchester_binary[0])
        manchester_binary_processing_array = process_MBA(manchester_binary[0],manchester_binary_processing_array)
        for i in range(1,len(manchester_binary)):
            manchester_binary[i] = bool_to_int(bool(binary_data[i] == 1) ^ bool(binary_data[i-1] == 1))
            #print(manchester_binary[i])
            manchester_binary_processing_array = process_MBA(manchester_binary[i],manchester_binary_processing_array)

        manchester_binary_state = binary_data[len(binary_data) - 1]

        








        


        


        # # # # if(block_count == 1):
        # # # #     time = np.arange(0,round(len(RRC_I)),1,int)
        # # # #     ax96.clear()
        # # # #     ax96.plot(time, RRC_I)
        # # # #     ax96.set_ylabel('ax96.set_ylabel')
        # # # #     ax96.set_xlabel('ax96.set_xlabel')
        # # # #     ax96.set_title('ax96.set_title')

        # # # #     ax961.clear()
        # # # #     ax961.plot(time, RRC_Q)
        # # # #     ax961.set_ylabel('ax96.set_ylabel')
        # # # #     ax961.set_xlabel('ax96.set_xlabel')
        # # # #     ax961.set_title('ax96.set_title')

        # # # #     time = np.arange(0,round(len(ncoOut_I[0:10])),1,int)
        # # # #     ax962.clear()
        # # # #     ax962.plot(time, ncoOut_I[0:10])
        # # # #     ax962.plot(time, ncoOut_Q[0:10])
        # # # #     ax962.set_ylabel('ax96.set_ylabel')
        # # # #     ax962.set_xlabel('ax96.set_xlabel')
        # # # #     ax962.set_title('ax96.set_title')

        # # # #     ax963.clear()
        # # # #     ax963.plot(time, ncoOut_Q[0:10])
        # # # #     ax963.set_ylabel('ax96.set_ylabel')
        # # # #     ax963.set_xlabel('ax96.set_xlabel')
        # # # #     ax963.set_title('ax96.set_title')

        # # # #     fig1.savefig("../data/mrunobars5RRC_I" + ".png")


            

            
            
    
        # # # #     #np.savetxt('zero_array.out', zero_array, delimiter=',', fmt='%d  %1.5e')

        # # # #     indx = np.arange(0,RRC_I_samples.size,1,int)
        # # # #     indxd_data = np.column_stack((indx,RRC_I_samples))
    
        # # # #     np.savetxt('RRC_I_samples.out', indxd_data, delimiter=',', fmt='%d  %1.5e')

        # # # #     indx = np.arange(0,RRC_Q_samples.size,1,int)
        # # # #     indxd_data = np.column_stack((indx,RRC_Q_samples))
    
        # # # #     np.savetxt('RRC_Q_samples.out', indxd_data, delimiter=',', fmt='%d  %1.5e')




        # # # #     indx = np.arange(0,RDS_I_filt_ds.size,1,int)
        # # # #     indxd_data = np.column_stack((indx,RDS_I_filt_ds))
    
        # # # #     np.savetxt('RDS_I_filt_ds.out', indxd_data, delimiter=',', fmt='%d  %1.5e')

        # # # #     indx = np.arange(0,RDS_Q_filt_ds.size,1,int)
        # # # #     indxd_data = np.column_stack((indx,RDS_Q_filt_ds))
    
        # # # #     np.savetxt('RDS_Q_filt_ds.out', indxd_data, delimiter=',', fmt='%d  %1.5e')

        # # # #     indx = np.arange(0,RDS_I_filt.size,1,int)
        # # # #     indxd_data = np.column_stack((indx,RDS_I_filt))
    
        # # # #     np.savetxt('RDS_I_filt.out', indxd_data, delimiter=',', fmt='%d  %1.5e')

        # # # #     indx = np.arange(0,RDS_Q_filt.size,1,int)
        # # # #     indxd_data = np.column_stack((indx,RDS_Q_filt))
    
        # # # #     np.savetxt('RDS_Q_filt.out', indxd_data, delimiter=',', fmt='%d  %1.5e')

        # # # #     indx = np.arange(0,RRC_I.size,1,int)
        # # # #     indxd_data = np.column_stack((indx,RRC_I))
    
        # # # #     np.savetxt('RRC_I.out', indxd_data, delimiter=',', fmt='%d  %1.5e')

        # # # #     indx = np.arange(0,RRC_Q.size,1,int)
        # # # #     indxd_data = np.column_stack((indx,RRC_Q))
    
        # # # #     np.savetxt('RRC_Q.out', indxd_data, delimiter=',', fmt='%d  %1.5e')

        # # # #     indx = np.arange(0,ncoOut_I.size,1,int)
        # # # #     indxd_data = np.column_stack((indx,ncoOut_I))
    
        # # # #     np.savetxt('ncoOut_I.out', indxd_data, delimiter=',', fmt='%d  %1.5e')

        # # # #     indx = np.arange(0,ncoOut_Q.size,1,int)
        # # # #     indxd_data = np.column_stack((indx,ncoOut_Q))
    
        # # # #     np.savetxt('ncoOut_Q.out', indxd_data, delimiter=',', fmt='%d  %1.5e')

        # # # #     indx = np.arange(0,reco_filt.size,1,int)
        # # # #     indxd_data = np.column_stack((indx,reco_filt))
    
        # # # #     np.savetxt('reco_filt.out', indxd_data, delimiter=',', fmt='%d  %1.5e')

        # # # #     indx = np.arange(0,RDS_I.size,1,int)
        # # # #     indxd_data = np.column_stack((indx,RDS_I))
    
        # # # #     np.savetxt('RDS_I.out', indxd_data, delimiter=',', fmt='%d  %1.5e')

        # # # #     indx = np.arange(0,RDS_Q.size,1,int)
        # # # #     indxd_data = np.column_stack((indx,RDS_Q))
    
        # # # #     np.savetxt('RDS_Q.out', indxd_data, delimiter=',', fmt='%d  %1.5e')

        # # # #     indx = np.arange(0,tone_pass_coeff.size,1,int)
        # # # #     indxd_data = np.column_stack((indx,tone_pass_coeff))
    
        # # # #     np.savetxt('tone_pass_coeff.out', indxd_data, delimiter=',', fmt='%d  %1.5e')

        # # # #     indx = np.arange(0,stereo_pass_coeff.size,1,int)
        # # # #     indxd_data = np.column_stack((indx,stereo_pass_coeff))
    
        # # # #     np.savetxt('stereo_pass_coeff.out', indxd_data, delimiter=',', fmt='%d  %1.5e')

        # # # #     indx = np.arange(0,state_phase.size,1,int)
        # # # #     indxd_data = np.column_stack((indx,state_phase))
    
        # # # #     np.savetxt('state_phase.out', indxd_data, delimiter=',', fmt='%d  %1.5e')

        # # # #     indx = np.arange(0,state_lpf_16k.size,1,int)
        # # # #     indxd_data = np.column_stack((indx,state_lpf_16k))
    
        # # # #     np.savetxt('state_lpf_16k.out', indxd_data, delimiter=',', fmt='%d  %1.5e')

        # # # #     indx = np.arange(0,state_q_lpf_100k.size,1,int)
        # # # #     indxd_data = np.column_stack((indx,state_q_lpf_100k))
    
        # # # #     np.savetxt('state_q_lpf_100k.out', indxd_data, delimiter=',', fmt='%d  %1.5e')

        # # # #     indx = np.arange(0,state_i_lpf_100k.size,1,int)
        # # # #     indxd_data = np.column_stack((indx,state_i_lpf_100k))
    
        # # # #     np.savetxt('state_i_lpf_100k.out', indxd_data, delimiter=',', fmt='%d  %1.5e')


        # # # #     indx = np.arange(0,len(audio_data),1,int)
        # # # #     indxd_data = np.column_stack((indx,audio_data))
    
        # # # #     np.savetxt('audio_data.out', indxd_data, delimiter=',', fmt='%d  %1.5e')


        # # # #     indx = np.arange(0,len(audio_block),1,int)
        # # # #     indxd_data = np.column_stack((indx,audio_block))
    
        # # # #     np.savetxt('audio_block.out', indxd_data, delimiter=',', fmt='%d  %1.5e')

        # # # #     firstIQ = iq_data[0:block_size]

        # # # #     indx = np.arange(0,len(firstIQ),1,int)
        # # # #     indxd_data = np.column_stack((indx,firstIQ))
    
        # # # #     np.savetxt('iq.out', indxd_data, delimiter=',', fmt='%d  %1.5e')

        # # # #     isamples_c = iq_data[(block_count)*block_size:(block_count+1)*block_size:2]
        # # # #     qsamples_c = iq_data[(block_count)*block_size+1:(block_count+1)*block_size:2]

        # # # #     indx = np.arange(0,len(isamples_c),1,int)
        # # # #     indxd_data = np.column_stack((indx,isamples_c))
    
        # # # #     np.savetxt('isamples.out', indxd_data, delimiter=',', fmt='%d  %1.5e')

        # # # #     indx = np.arange(0,len(qsamples_c),1,int)
        # # # #     indxd_data = np.column_stack((indx,qsamples_c))
    
        # # # #     np.savetxt('qsamples.out', indxd_data, delimiter=',', fmt='%d  %1.5e')

        # # # #     indx = np.arange(0,len(i_ds),1,int)
        # # # #     indxd_data = np.column_stack((indx,i_ds))
    
        # # # #     np.savetxt('ids.out', indxd_data, delimiter=',', fmt='%d  %1.5e')

        # # # #     indx = np.arange(0,len(q_ds),1,int)
        # # # #     indxd_data = np.column_stack((indx,q_ds))
    
        # # # #     np.savetxt('qds.out', indxd_data, delimiter=',', fmt='%d  %1.5e')

        # # # #     indx = np.arange(0,len(fm_demod),1,int)
        # # # #     indxd_data = np.column_stack((indx,fm_demod))
    
        # # # #     np.savetxt('fm_demod.out', indxd_data, delimiter=',', fmt='%d  %1.5e')
        # # # #     #break

        # # # # # to save runtime select the range of blocks to log data
        # # # # # this includes both saving binary files as well plotting PSD
        # # # # # below we assume we want to plot for graphs for blocks 10 and 11
        # # # # if block_count >= 10 and block_count < 20:
        # # # #     # PSD after FM demodulation
        # # # #     ax0.clear()
        # # # #     ax0.psd(fm_demod, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
        # # # #     ax0.set_ylabel('PSD (dB/Hz)')
        # # # #     ax0.set_xlabel('Freq (kHz)')
        # # # #     ax0.set_title('Demodulated FM (block ' + str(block_count) + ')')
        # # # #     # output binary file name (where samples are written from Python)
        # # # #     fm_demod_fname = "../data/fm_demod_" + str(block_count) + ".bin"
        # # # #     # create binary file where each sample is a 32-bit float
        # # # #     fm_demod.astype('float32').tofile(fm_demod_fname)

        # # # #     # PSD after extracting mono audio
        # # # #     ax1.clear()
        # # # #     ax1.scatter(RRC_I_samples,RRC_Q_samples,s=10)
        # # # #     # ax1.psd(fm_demod, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
        # # # #     ax1.set_ylabel('PSD (dB/Hz)')
        # # # #     ax1.set_xlabel('Freq (kHz)')
        # # # #     ax1.set_title('Extracted tone')

        # # # #     # PSD after decimating mono audio
        # # # #     ax2.clear()
        # # # #     ax2.psd(extr_filt, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
        # # # #     ax2.set_ylabel('PSD (dB/Hz)')
        # # # #     ax2.set_xlabel('Freq (kHz)')
        # # # #     ax2.set_title('Extracted stereo')

        # # # #     # save figure to file
        # # # #     fig.savefig("../data/fmAudioRDS" + str(block_count) + ".png", dpi=300)


        # # # #     ax123.clear()
        # # # #     ax123.scatter(RRC_I_samples,RRC_Q_samples,s=10)
        # # # #     # ax1.psd(fm_demod, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
        # # # #     ax123.set_ylabel('Q')
        # # # #     ax123.set_xlabel('I')
        # # # #     ax123.set_title('Extracted Constellation')


        # # # #     fig2.savefig("../data/mrunoConstellation-"+"SampleSpacing-"+str(sample_spacing)+"-SampleShift-"+str(sample_shift)+"-Block-" +str(block_count) + ".png", dpi=500)

        # # # #     time = np.arange(0,round(len(RRC_I)),1,int)

        # # # #     axi.clear()
        # # # #     axi.plot(time, RRC_I)
        # # # #     # ax1.psd(fm_demod, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
        # # # #     axi.set_ylabel('I')
        # # # #     axi.set_xlabel('Time')
        # # # #     axi.set_title('I wave')

        # # # #     axq.clear()
        # # # #     axq.plot(time, RRC_Q)
        # # # #     # ax1.psd(fm_demod, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
        # # # #     axq.set_ylabel('Q')
        # # # #     axq.set_xlabel('Time')
        # # # #     axq.set_title('Q wave')

        # # # #     spaxi.clear()
        # # # #     spaxi.plot(time, RRC_I_sparse)
        # # # #     # ax1.psd(fm_demod, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
        # # # #     spaxi.set_ylabel('I')
        # # # #     spaxi.set_xlabel('Time')
        # # # #     spaxi.set_title('I wave')

        # # # #     spaxq.clear()
        # # # #     spaxq.plot(time, RRC_Q_sparse)
        # # # #     # ax1.psd(fm_demod, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
        # # # #     spaxq.set_ylabel('Q')
        # # # #     spaxq.set_xlabel('Time')
        # # # #     spaxq.set_title('Q wave')



        # # # #     iwave.savefig("../data/iwave"+ str(block_count) + ".png", dpi=300)

        # # # #     qwave.savefig("../data/qwave"+ str(block_count) + ".png", dpi=300)

        # # # #     spiwave.savefig("../data/spiwave"+ str(block_count) + ".png", dpi=300)

        # # # #     spqwave.savefig("../data/spqwave"+ str(block_count) + ".png", dpi=300)


        # # # #     indx = np.arange(0,RRC_I_samples.size,1,int)
        # # # #     indxd_data = np.column_stack((indx,RRC_I_samples))
    
        # # # #     np.savetxt('RRC_I_samples'+ str(block_count) + ".out", indxd_data, delimiter=',', fmt='%d  %1.5e')

        # # # #     indx = np.arange(0,RRC_Q_samples.size,1,int)
        # # # #     indxd_data = np.column_stack((indx,RRC_Q_samples))
    
        # # # #     np.savetxt('RRC_Q_samples'+ str(block_count) + ".out", indxd_data, delimiter=',', fmt='%d  %1.5e')

        block_count += 1

    print('Finished processing the raw I/Q samples')

    reshaped_data = np.reshape(audio_data_s, (-1, 2))
    # write audio data to a .wav file (assumes audio_data samples are -1 to +1)
    wavfile.write("../data/fmAudio2eroster.wav", int(audio_Fs), np.int16((reshaped_data/2)*32767))

    # uncomment assuming you wish to show some plots
    # plt.show()
