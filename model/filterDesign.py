import matplotlib.pyplot as plt
from scipy import signal
import numpy as np
import os

print(os.listdir())
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# use generateSin/plotTime from the fourierTransform module
from fourierTransform import generateSin, plotTime

def cSinc(x):
	return np.sin(x*np.pi)/(np.pi*x)

def cFFT(x):
	N = len(x)
	m = 0
	k = 0
	X = np.zeros((N,), dtype=np.complex_)
	for m in range(0,N-1,1):
		X[m] = (0 + 0j)
		for k in range(0, N-1, 1):
			X[m] = X[m] + x[k]*np.exp((-2*1j*np.pi*k*m)/N)
	return X

def cIFFT(X):
    N = len(X)
    m = 0
    k = 0
    x = np.zeros((N,), dtype=np.complex_)
    for k in range(0,N-1,1):
        x[k] = (0 + 0j)
        for m in range(0, N-1, 1):
            x[k] = x[k] + X[m]*np.exp((2*1j*np.pi*k*m)/N)
    return x

def cLfilter2(b,a,x):
	fb = cFFT(b)
	fx = cFFT(x)
	fy = []
	for i in range(0,len(fb),1):
		fy.append(fx[i]*fb[i])
	return cIFFT(fy)




def cFirwin(N_Taps, cutoff = None, window='Hann', fc = None, fs = None):
	if(window == 'Hann' or window == 'hann'):
		if(fc != None and fs != None and cutoff == None):
			cutoff = fc/(fs/2)
			print('Using fs and fc')
		elif(fc != None and fs != None and cutoff != None):
			raise Exception('You Goofed')
		h = np.zeros((N_Taps,), dtype=np.complex_)
		i = 0
		for i in range(0,N_Taps, 1):
			if(i == (N_Taps-1)/2):
				h[i] = cutoff
			else:
				h[i] = cutoff * cSinc(cutoff*(i-(N_Taps - 1)/2))
			h[i] = h[i] * (np.sin(i*np.pi/N_Taps))**2
		# print(h)
		# print('\n')
		# print(signal.firwin(N_Taps,cutoff,window='hann'))
		l = signal.firwin(N_Taps,cutoff,window='hann')
		
		max = 0
		for k in range(0,N_Taps-1,1):
			if(h[k]-l[k] > max):
				max = h[k]-l[k]
		print(max)
		return h
	raise Exception('bad window')

def cFirwin2(N_Taps, cutoff = None, window='Hann', fc = None, fs = None, fo = 0,left = None, right = None):
	print()
	if(window == 'Hann' or window == 'hann'):
		if(fc != None and fs != None and cutoff == None):
			cutoff = fc/(fs/2)									
			print('Using fs and fc')
		elif(fc != None and fs != None and cutoff != None):
			raise Exception('You Goofed')
		h = np.zeros((N_Taps,), dtype=np.complex_)
		shift = 0
		if(left != None and right != None):
			diff = (right - left)
			cutoff = diff * (.1/.61)
			shift = (left + diff/2)*.16

		i = 0
		for i in range(0,N_Taps, 1):
			if(i == (N_Taps-1)/2):
				h[i] = cutoff * np.cos(shift*2*np.pi*(i-(N_Taps - 1)/2))
			else:
				h[i] = cutoff * cSinc(cutoff*(i-(N_Taps - 1)/2)) * np.cos(shift*2*np.pi*(i-(N_Taps - 1)/2))
			h[i] = h[i] * (np.sin(i*np.pi/N_Taps))**2
		# hn = np.zeros((N_Taps + 2*fo,), dtype=np.complex_)
		# for l in range(0,N_Taps,1):
		# 	hn[l+fo] = h[l]
		# h = hn
		print(h)
		# print('\n')
		# print(signal.firwin(N_Taps,cutoff,window='hann'))
		l = cFirwin(N_Taps,cutoff,window='hann')
		# max = 0
		# for k in range(0,N_Taps-1,1):
		# 	if(h[k]-l[k] > max):
		# 		max = h[k]-l[k]
		# print('ho',max,'hi')
		return h
	raise Exception('bad window')



def cLfilter(b,a,x):
	# if(len(b) > len(x)):
	# 	raise Exception('More Taps than Samples')
	#print(b)
	if(not isinstance(a,list)):
		a_temp = np.zeros((len(b),),dtype=float)
		for i in range(0,len(a_temp),1):
			a_temp[i] = float(a)
		a = a_temp
	y = np.zeros((len(x),),dtype=float)
	#print(a)
	for n in range(0,len(x),1):
		y[n] = 0
		for k in range(0,n,1):
				#print(n,k,n-k)
				if(k >= len(x) or k >= len(b)):
					break
				y[n] = y[n] + b[k]*x[n-k]
				#print(x[n-k],n-k)
			#print(y[n],n)
		y[n] = y[n]/a[0]
		
	return y

	
		

	

def freqzPlot(coeff, msg):

	# find the frequency response using freqz from SciPy:
	# https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.freqz.html
    w, h = signal.freqz(coeff)

	# plots the magnitude response where the x axis is normalized in rad/sample
	# Reminder: math.pi rad/sample is actually the Nyquist frequency
    fig, ax1 = plt.subplots()
    ax1.set_title('Digital filter frequency response (' + msg + ')')
    ax1.plot(w, 20 * np.log10(abs(h)), 'b')
    ax1.set_ylabel('Amplitude [dB]', color='b')
    ax1.set_xlabel('Frequency [rad/sample]')

	# uncomment the lines below if you wish to inspect the phase response
	# Note: as important as the phase response is, it is not critical
	# at this stage because we expect a linear phase in the passband

    ax2 = ax1.twinx()
    angles = np.unwrap(np.angle(h))
    ax2.plot(w, angles, 'g')
    ax2.set_ylabel('Angle (radians)', color='g')

def filterSin(Fs, Fc, coeff):

    # we can control the frequency relative to the filter cutoff
	time, x = generateSin(Fs, interval = 1.0, frequency = Fc * 1.2)
	plotTime(x, time)

    # use lfilter from SciPy for FIR filtering:
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.lfilter.html
	fx = cLfilter(coeff, 1.0, x)

    # you should cleary the effects (attenuation, delay) introduced by the filter
	plotTime(fx, time)

if __name__ == "__main__":

	Fs = 100.0           # sampling rate
	Fc = 15.0            # cutoff frequency
	N_taps = 90          # number of taps for the FIR

	# derive filter coefficients using firwin from Scipy:
	# https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.firwin.html
	# second argument is the normalized cutoff frequency, i.e., the
	# cutoff frequency divided by Nyquist frequency (half of sampling rate)
	firwin_coeff = cFirwin2(N_taps, Fc/(Fs/2), window=('hann') )
	

	

	#plt.plot(firwin_coeff)
	# plot the frequency response obtained through freqz
	#freqzPlot(firwin_coeff, 'firwin with ' + str(N_taps) + ' taps')

    # implement your own method for finding the coefficients for a low pass filter
    # my_own_coeff = ... provide the following arguments: Fc, Fs and N_taps
    # compare through visual inspection the frequency response against firwin
	# freqzPlot(my_own_coeff, 'my own FIR design with ' + str(N_taps) + ' taps')

	# you can confirm that a single tone has been filtered
	filterSin(Fs, Fc, firwin_coeff)

	plt.show()
