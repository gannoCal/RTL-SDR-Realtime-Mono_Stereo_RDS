import matplotlib.pyplot as plt
import numpy as np
from scipy.io import wavfile
from scipy import signal
import os

print(os.listdir())
os.chdir(os.path.dirname(os.path.abspath(__file__)))

def cSinc(x):
	return np.sin(x*np.pi)/(np.pi*x)


def cLfilter(b,a,x,bufferData = []):
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
	if(len(bufferData) != 0):
		special = 0
		for n in range(0,len(x),1):
			y[n] = 0
			special = 0
			for k in range(0,len(b),1):
					#print(n,k,n-k)
					if(k >= len(x) or k >= len(b)):
						break

					if(n-k >= 0):
						y[n] = y[n] + b[k]*x[n-k]
					elif (n-k < 0):
						y[n] = y[n] + b[k]*bufferData[len(bufferData) - 1 - special]
						special = special + 1
				#print(y[n],n)
			y[n] = y[n]/a[0]
	else:
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

def cfilter_block_processing(audio_data, \
							block_size, \
							audio_Fc, \
							audio_Fs, \
							N_taps):

	# derive filter coefficients
	firwin_coeff = cFirwin(N_taps, audio_Fc/audio_Fs)
	# we assume the data is stereo as in the audio test file
	filtered_data = np.empty(shape = audio_data.shape)
	# start at the first block (with relative position zero)
	position = 0
	buffer = np.zeros((block_size,),dtype=float)
	while True:

		# filter both left and right channels
		filtered_data[position:position+block_size, 0] = \
			cLfilter(firwin_coeff, 1.0, \
			audio_data[position:position+block_size, 0],bufferData = buffer)
		filtered_data[position:position+block_size, 1] = \
			cLfilter(firwin_coeff, 1.0, \
			audio_data[position:position+block_size, 1],bufferData = buffer)
		buffer = audio_data[position:position+block_size, 1]
		position += block_size
		if position > len(audio_data):
			break

	# to properly handle blocks you will need to use
	# the zi argument from lfilter from SciPy
	# explore SciPy, experiment, understand and learn!

	return filtered_data

def cfilter_single_pass(audio_data, audio_Fc, audio_Fs, N_taps):

	# derive filter coefficients
	firwin_coeff = cFirwin(N_taps, audio_Fc/audio_Fs)
	# we assume the data is stereo as in the audio test file
	filtered_data = np.empty(shape = audio_data.shape)
	# filter left channel
	filtered_data[:,0] = cLfilter(firwin_coeff, 1.0, audio_data[:,0])
	# filter stereo channel
	filtered_data[:,1] = cLfilter(firwin_coeff, 1.0, audio_data[:,1])

	return filtered_data




def filter_block_processing(audio_data, \
							block_size, \
							audio_Fc, \
							audio_Fs, \
							N_taps):

	# derive filter coefficients
	firwin_coeff = signal.firwin(N_taps, audio_Fc/audio_Fs)
	# we assume the data is stereo as in the audio test file
	filtered_data = np.empty(shape = audio_data.shape)
	# start at the first block (with relative position zero)
	position = 0

	while True:

		# filter both left and right channels
		filtered_data[position:position+block_size, 0] = \
			signal.lfilter(firwin_coeff, 1.0, \
			audio_data[position:position+block_size, 0])
		filtered_data[position:position+block_size, 1] = \
			signal.lfilter(firwin_coeff, 1.0, \
			audio_data[position:position+block_size, 1])

		position += block_size
		if position > len(audio_data):
			break

	# to properly handle blocks you will need to use
	# the zi argument from lfilter from SciPy
	# explore SciPy, experiment, understand and learn!

	return filtered_data

def filter_single_pass(audio_data, audio_Fc, audio_Fs, N_taps):

	# derive filter coefficients
	firwin_coeff = signal.firwin(N_taps, audio_Fc/audio_Fs)
	# we assume the data is stereo as in the audio test file
	filtered_data = np.empty(shape = audio_data.shape)
	# filter left channel
	filtered_data[:,0] = signal.lfilter(firwin_coeff, 1.0, audio_data[:,0])
	# filter stereo channel
	filtered_data[:,1] = signal.lfilter(firwin_coeff, 1.0, audio_data[:,1])

	return filtered_data

# audio test file from: https://www.videvo.net/royalty-free-music/
if __name__ == "__main__":

	# use use wavfile from scipy.io for handling .wav files
	audio_Fs, audio_data = wavfile.read("../data/audio_test.wav")
	print(' Audio sample rate = {0:f} \
		\n Number of channels = {1:d} \
		\n Numbef of samples = {2:d}' \
		.format(audio_Fs, audio_data.ndim, len(audio_data)))
	#audio_data = np.zeros(audio_data.shape,dtype=float)
	# you can control the cutoff frequency and number of taps
	single_pass_data = filter_single_pass(audio_data, \
						audio_Fc = 10e3, \
						audio_Fs = audio_Fs, \
						N_taps = 51)
	print(single_pass_data)
	# write filtered data back to a .wav file
	wavfile.write("../data/single_pass_filtered421.wav", \
	 			audio_Fs, \
				single_pass_data.astype(np.int16))
	#print(single_pass_data.astype(np.int16))
	# you can control also the block size
	block_processing_data = cfilter_block_processing(audio_data, \
						block_size = 1000, \
						audio_Fc = 10e3, \
						audio_Fs = audio_Fs, \
						N_taps = 51)

	wavfile.write("../data/block_processing_filtered421.wav", \
	 			audio_Fs, \
				block_processing_data.astype(np.int16))

	# it is suggested that you add plotting while troubleshooting
	# if you plot in the time domain, select a subset of samples,
	# from a particular channel (or both channels) e.g.,
	# audio_data[start:start+number_of_samples, 0]
	plt.show()
    
