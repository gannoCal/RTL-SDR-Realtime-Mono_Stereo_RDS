import math
import numpy as np
from scipy import signal

Fs_in = 2.4e6
Fs_out = 2.5e6

if __name__ == "__main__":

	in_fname = "../data/samples/my_samples_2.4_march16.raw"
	raw_data = np.fromfile(in_fname, dtype='uint8')
	# IQ data is normalized between -1 and +1
	iq_data = (raw_data - 128.0)/128.0
	print("Read raw RF data from \"" + in_fname + "\" in unsigned 8-bit format")

	expand = int(Fs_out) / math.gcd(int(Fs_in), int(Fs_out))
	decim = int(Fs_in) / math.gcd(int(Fs_in), int(Fs_out))

	# the resampling is done on the full recording
	resampled_i = signal.resample_poly(iq_data[0::2], expand, decim)
	resampled_q = signal.resample_poly(iq_data[1::2], expand, decim)

	min_i = np.min(resampled_i)
	max_i = np.max(resampled_i)
	min_q = np.min(resampled_q)
	max_q = np.max(resampled_q)

	print("min_i={0:.2f} max_i={1:.2f} min_q={2:.2f} max_q={3:.2f}"\
		.format(min_i, max_i, min_q, max_q))

	# create 8-bit integers from floats
	out_data = np.empty(2*len(resampled_i), dtype='uint8')
	for k in range(len(resampled_i)):
		rescale_i = int(int(math.floor((resampled_i[k]-min_i)*255.0/(max_i-min_i))))
		# if rescale_i < 0 or rescale_i > 255:
		# 	print("i/{0:d}: {1:d}" .format(k, rescale_i))
		out_data[2*k] = rescale_i

		rescale_q = int(int(math.floor((resampled_q[k]-min_q)*255.0/(max_q-min_q))))
		# if rescale_q < 0 or rescale_q > 255:
		# 	print("i/{0:d}: {1:d}" .format(k, rescale_q))
		out_data[2*k+1] = rescale_q

	# write resampled IQ data as 8-bit unsigned
	out_fname = "../data/samples_" + str(int(Fs_out/100e3)) + ".raw"
	out_data.astype('uint8').tofile(out_fname)
	print("Written resampled RF data from \"" + out_fname + "\" in unsigned 8-bit format")
