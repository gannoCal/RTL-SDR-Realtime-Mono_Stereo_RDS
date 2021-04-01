#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import numpy as np
import math, cmath

#
# you should use the demodulator based on arctan given below as a reference
#
# in order to implement your OWN FM demodulator without the arctan function,
# a very good and to-the-point description is given by Richard Lyons at:
#
# https://www.embedded.com/dsp-tricks-frequency-demodulation-algorithms/
#
# the demodulator boils down to implementing equation (13-117) from above, where
# the derivatives are nothing else but differences between consecutive samples
#
# needless to say, you should not jump directly to equation (13-117)
# rather try first to understand the entire thought process based on calculus
# identities, like derivative of the arctan function or derivatives of ratios
#

#
# use the four quadrant arctan function for phase detect between a pair of
# IQ samples; then unwrap the phase and take its derivative to demodulate
#
def fmDemodArctan(I, Q, prev_phase = 0.0):
#
# the default prev_phase phase is assumed to be zero, however
# take note in block processing it must be explicitly controlled

    # empty vector to store the demodulated samples
    fm_demod = np.empty(len(I))

    # iterate through each of the I and Q pairs
    for k in range(len(I)):

        # use the atan2 function (four quadrant version) to detect angle between
        # the imaginary part (quadrature Q) and the real part (in-phase I)
        current_phase = math.atan2(Q[k], I[k])

        # we need to unwrap the angle obtained in radians through arctan2
        # to deal with the case when the change between consecutive angles
        # is greater than Pi radians (unwrap brings it back between -Pi to Pi)
        [prev_phase, current_phase] = np.unwrap([prev_phase, current_phase])

        # take the derivative of the phase
        fm_demod[k] = current_phase - prev_phase

        # save the state of the current phase
        # to compute the next derivative
        prev_phase = current_phase

    # return both the demodulated samples as well as the last phase
    # (the last phase is needed to enable continuity for block processing)
    return fm_demod, prev_phase

# custom function for DFT that can be used by the PSD estimate
def DFT(x):

    # number of samples
    N = len(x)

    # frequency bins
    Xf = np.zeros(N, dtype='complex')

    # iterate through all frequency bins/samples
    for m in range(N):
        for k in range(N):
            Xf[m] += x[k] * cmath.exp(1j * 2 * math.pi * ((-k) * m) / N)

    # return the vector that holds the frequency bins
    return Xf

# custom function to estimate PSD based on the Bartlett method
# this is less accurate than the Welch method from matplotlib
# however, as the visual inspections confirm, the estimate gives
# the user a "reasonably good" view of the power spectrum
def estimatePSD(samples, NFFT, Fs):

    # rename the NFFT argument (notation consistent with matplotlib.psd)
    # to freq_bins (i.e., frequency bins for which we compute the spectrum)
    freq_bins = NFFT
    # frequency increment (or resolution)
    df = Fs/freq_bins

    # create the frequency vector to be used on the X axis
    # for plotting the PSD on the Y axis (only positive freq)
    freq = np.arange(0, Fs/2, df)

    # design the Hann window used to smoothen the discrete data in order
    # to reduce the spectral leakage after the Fourier transform
    hann = np.empty(freq_bins)
    for i in range(len(hann)):
        hann[i] = pow(math.sin(i*math.pi/freq_bins),2)

    # create an empty list where the PSD for each segment is computed
    psd_list = []

    # samples should be a multiple of frequency bins, so
    # the number of segments used for estimation is an integer
    # note: for this to work you must provide an argument for the
    # number of frequency bins not greater than the number of samples!
    no_segments = int(math.floor(len(samples)/float(freq_bins)))

    # iterate through all the segments
    for k in range(no_segments):

        # apply the hann window (using pointwise multiplication)
        # before computing the Fourier transform on a segment
        windowed_samples = samples[k*freq_bins:(k+1)*freq_bins] * hann

        # compute the Fourier transform using the built-in FFT from numpy
        Xf = np.fft.fft(windowed_samples, freq_bins)

        # note, you can check how MUCH slower is DFT vs FFT by replacing the
        # above function call with the one that is commented below
        #
        # Xf = DFT(windowed_samples)
        #
        # note, the slow impelementation of the Fourier transform is not as
        # critical when computing a static power spectra when troubleshooting
        # note: time permitting a custom FFT can be implemented

        # since input is real, we keep only the positive half of the spectrum
        # however, we will also add the signal energy of negative frequencies
        # to have a better a more accurate PSD estimate when plotting
        Xf = Xf[0:int(freq_bins/2)] # keep only positive freq bins
        psd_seg = 1/(Fs*freq_bins/2) * abs(Xf)**2 # compute signal power
        psd_seg = 2*psd_seg # add the energy from the negative freq bins

        # translate to the decibel (dB) scale
        for i in range(len(psd_seg)):
            psd_seg[i] = 10*math.log10(psd_seg[i])

        # append to the list where PSD for each segment is stored
        # in sequential order (first segment, followed by the second one, ...)
        psd_list.extend(psd_seg)

    # compute the estimate to be returned by the function through averaging
    psd_est = np.zeros(int(freq_bins/2))

    # iterate through all the frequency bins (positive freq only)
    # from all segments and average them (one bin at a time ...)
    for k in range(int(freq_bins/2)):
        # iterate through all the segments
        for l in range(no_segments):
            psd_est[k] += psd_list[k + l*int(freq_bins/2)]
        # compute the estimate for each bin
        psd_est[k] = psd_est[k] / no_segments

    # the frequency vector and PSD estimate
    return freq, psd_est





def fmPll(pllIn, freq, Fs, \
    ncoScale = 1.0, phaseAdjust = 0.0, normBandwidth = 0.01, \
    state_integrator = 0.0, \
    state_phaseEst = 0.0, \
    state_feedbackI = 1.0, \
    state_feedbackQ = 0.0, \
    state_ncoOut0 = 1.0, \
    state_trigOffset = 0):
        
            
    """
    pllIn              array of floats
                    input signal to the PLL (assume known frequency)

    freq             float
                    reference frequency to which the PLL locks

    Fs              float
                    sampling rate for the input/output signals

    ncoScale        float
                    frequency scale factor for the NCO output

    phaseAdjust        float
                    phase adjust to be added to the NCO only

    normBandwidth    float
                    normalized bandwidth for the loop filter
                    (relative to the sampling rate)

    state             to be added

    """

    # scale factors for proportional/integrator terms
    # these scale factors were derived assuming the following:
    # damping factor of 0.707 (1 over square root of 2)
    # there is no oscillator gain and no phase detector gain
    Cp = 2.666
    Ci = 3.555

    # gain for the proportional term
    Kp = (normBandwidth)*Cp
    # gain for the integrator term
    Ki = (normBandwidth*normBandwidth)*Ci

    # output array for the NCO
    ncoOut = np.empty(len(pllIn)+1)

    # initialize internal state
    integrator = state_integrator
    phaseEst = state_phaseEst
    feedbackI = state_feedbackI
    feedbackQ = state_feedbackQ
    ncoOut[0] = state_ncoOut0
    trigOffset = state_trigOffset
    # note: state saving will be needed for block processing

    for k in range(len(pllIn)):

        # phase detector
        errorI = pllIn[k] * (+feedbackI)  # complex conjugate of the
        errorQ = pllIn[k] * (-feedbackQ)  # feedback complex exponential

        # four-quadrant arctangent discriminator for phase error detection
        errorD = math.atan2(errorQ, errorI)

        # loop filter
        integrator = integrator + Ki*errorD

        # update phase estimate
        phaseEst = phaseEst + Kp*errorD + integrator
        # internal oscillator
        trigArg = 2*math.pi*(freq/Fs)*(trigOffset+k+1) + phaseEst
        feedbackI = math.cos(trigArg)
        feedbackQ = math.sin(trigArg)
        ncoOut[k+1] = math.cos(trigArg*ncoScale + phaseAdjust)



    state_integrator = integrator
    state_phaseEst = phaseEst
    state_feedbackI = feedbackI 
    state_feedbackQ = feedbackQ 
    state_ncoOut0 = ncoOut[len(ncoOut)-1]
    state_trigOffset = (trigOffset + len(pllIn))
    
    
    # for stereo only the in-phase NCO component should be returned
    # for block processing you should also return the state
    return ncoOut, \
    state_integrator, \
    state_phaseEst, \
    state_feedbackI, \
    state_feedbackQ, \
    state_ncoOut0, \
    state_trigOffset
    # for RDS add also the quadrature NCO component to the output



def impulseResponseRootRaisedCosine(Fs, N_taps):

    """
    Root raised cosine (RRC) filter

    Fs          sampling rate at the output of the resampler in the RDS path
                sampling rate must be an integer multipler of 2375
                this integer multiple is the number of samples per symbol

    N_taps      number of filter taps

    """

    # duation for each symbol - do NOT be changed for RDS!
    T_symbol = 1/2375.0

    # roll-off factor (greater than 0 and smaller than 1)
    beta = 0.90

    # the RRC inpulse response that will be computed in this function
    impulseResponseRRC = np.empty(N_taps)

    for k in range(N_taps):
        t = float((k-N_taps/2))/Fs
        #print(N_taps/2)
        # we ignore the 1/T_symbol scale factor
        if t == 0.0: impulseResponseRRC[k] = 1.0 + beta*((4/math.pi)-1)
        elif t == -T_symbol/(4*beta) or t == T_symbol/(4*beta):
            impulseResponseRRC[k] = (beta/np.sqrt(2))*(((1+2/math.pi)* \
                    (math.sin(math.pi/(4*beta)))) + ((1-2/math.pi)*(math.cos(math.pi/(4*beta)))))
        else: impulseResponseRRC[k] = (math.sin(math.pi*t*(1-beta)/T_symbol) +  \
                    4*beta*(t/T_symbol)*math.cos(math.pi*t*(1+beta)/T_symbol))/ \
                    (math.pi*t*(1-(4*beta*t/T_symbol)*(4*beta*t/T_symbol))/T_symbol)

    return impulseResponseRRC



if __name__ == "__main__":

    # do nothing when this module is launched on its own
    pass
