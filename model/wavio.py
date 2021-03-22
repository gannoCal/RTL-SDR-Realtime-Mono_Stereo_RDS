import numpy as np
from scipy.io import wavfile

#
# since manipulating .wav files is not the objective of the SDR project and
# we are using them solely for "assessing" the outcome of the DSP tasks
# while troubleshooting, we will avoid processing any .wav files in C++,
# mainly because of the error prone nature of handling .wav file headers
#
# for the reason above, the Python script below can be used to parse/format
# .wav files to/from binary files where the sample representation is known
# (or better said agreed on) by both the Python script and the C++ program
#
# .wav files should be opened only in this Python script and samples written
# in binary (e.g., assuming 32-bit floating point for this example) should be
# read by the C++ program in binary format (raw data, no headers); subsequently,
# the C++ program should output the processed data also in binary formart,
# which can be read back by this Python script to be formatted properly with a
# a header into a .wav file that can then be used on a third part audio player
#

if __name__ == "__main__":

	# parse an audio file
	audio_Fs = int(48*10*10*10) #audio_data = wavfile.read("../data/blockfloat32filtered.wav")

	in_fname = "../data/new_binary69_test6.bin"

	float_data = np.fromfile(in_fname, dtype='float64')
	print(" Read binary data from \"" + in_fname + "\" in float64 format")

#	reshaped_data = np.reshape(float_data, (-1, 2))
# 	list1 = float_data.tolist()
# 	print(len(list1))
# 	for i in range(len(list1)):
# 		if(i > len(list1) - 1):
# 				break
# 		elif(list1[i] < 1e-5):
# 			list1.pop(i)
# 	print(len(list1))
# 	float_data = np.array(list1)


	wavfile.write("../data/new_binary69_test6.wav", \
				audio_Fs, \
				np.int16((float_data/2)*32767))

	# note: we can also dump audio data in other formats, if needed
	# audio_data.astype('int16').tofile('int16samples.bin')