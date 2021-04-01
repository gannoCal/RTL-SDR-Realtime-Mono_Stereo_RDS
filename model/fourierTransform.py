import matplotlib.pyplot as plt
import numpy as np
import cmath, math
import os

print(os.listdir())
os.chdir(os.path.dirname(os.path.abspath(__file__)))

def SignalEnergy(X,inverse = False):
    Total = 0
    for x in X:
        Total = Total + math.pow(abs(x),2)
    if(inverse):
        Total = Total/len(X)
    return Total

def cFFT(x):
    N = len(x)
    m = 0
    k = 0
    X = np.zeros((N,), dtype=np.complex_)
    for m in range(0,N-1,1):
        X[m] = (0 + 0j)
        for k in range(0, N-1, 1):
            X[m] = X[m] + x[k]*np.exp((-2*1j*math.pi*k*m)/N)
    return X

def cIFFT(X):
    N = len(X)
    m = 0
    k = 0
    x = np.zeros((N,), dtype=np.complex_)
    for k in range(0,N-1,1):
        x[k] = (0 + 0j)
        for m in range(0, N-1, 1):
            x[k] = x[k] + X[m]*np.exp((2*1j*math.pi*k*m)/N)
    return x

def plotSpectrum(x, Fs, type = 'FFT'):

    n = len(x)             # length of the signal
    df = Fs/n              # frequency increment (width of freq bin)

    # compute Fourier transform, its magnitude and normalize it before plotting
    if type == 'FFT':
        Xfreq = cFFT(x)
        print(SignalEnergy(x), SignalEnergy(Xfreq,True))
        #print(Xfreq)
    XMag = abs(Xfreq)/n

    # Note: because x is real, we keep only the positive half of the spectrum
    # Note also: half of the energy is in the negative half (not plotted)
    XMag = XMag[0:int(n/2)]

    # freq vector up to Nyquist freq (half of the sample rate)
    freq = np.arange(0, Fs/2, df)

    fig, ax = plt.subplots()
    ax.plot(freq, XMag)
    ax.set(xlabel='Frequency (Hz)', ylabel='Magnitude',
        title='Frequency domain plot')
    # fig.savefig("freq.png")
    plt.show()

def plotTime(x, time):

    fig, ax = plt.subplots()
    ax.plot(time, x)
    ax.set(xlabel='Time (sec)', ylabel='Amplitude',
            title='Time domain plot')
    # fig.savefig("time.png")
    plt.show()

def generateSin(Fs, interval, frequency = 7.0, amplitude = 5.0, phase = 0.0):

    dt = 1.0/Fs                          # sampling period (increment in time)
    time = np.arange(0, interval, dt)    # time vector over interval

    # generate the sin signal
    x = amplitude*np.sin(2*math.pi*frequency*time+phase)

    return time, x

def generateSquare(Fs, interval , frequency = 4.0, amplitude = 5.0, phase = .05, duty = 0.6):

    dt = 1.0/Fs
    time = np.arange(0, interval, dt)

    # freq = np.sin(math.pi*time*frequency)/(math.pi*time*frequency)
    # freq[0] = 1
    # print(freq)
    # freq = np.zeros((len(time),), dtype=np.complex_)
    # freq[7] = 1
    # x = cIFFT(freq)
    print(1/frequency)
    localTime = time[0]
    x = np.zeros((len(time),), dtype=np.complex_)
    for i in range(0,len(time),1):
        localTime = (time[i]+phase)%(1/frequency)
        if(localTime < duty*(1/frequency)):
            x[i] = amplitude
        else:
            x[i] = 0
    
    return time, x

if __name__ == "__main__":

    Fs = 100.0          # sampling rate
    interval = 1.0      # set up to one full second

    # generate the user-defined sin function
    time, x = generateSquare(Fs, interval)
    # use np.random.randn() for randomization
    # we can owverwrie the default values
    # frequency =  8.0                     # frequency of the signal
    # amplitude =  3.0                     # amplitude of the signal
    # phase = 1.0                          # phase of the signal
    # time, x = generateSin(Fs, interval, frequency, amplitude, phase)

    # plot the signal in time domain
    plotTime(x, time)
    # plot the signal in frequency domain
    plotSpectrum(x, Fs, type = 'FFT')

    # compute the spectrum with your own DFT
    # you can use cmath.exp() for complex exponentials
    # plotSpectrum(x, Fs, type = 'your DFT name')

    # confirm DFT/IDFT correctness by checking if x == IDFT(DFT(x))
    # Note: you should also numerically check if the
    # signal energy in time and frequency domains is the same

    # generate randomized multi-tone signals
    # plot them in both time and frequency domain

    plt.show()
