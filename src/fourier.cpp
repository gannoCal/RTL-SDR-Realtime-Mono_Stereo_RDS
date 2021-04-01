/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

// source code for Fourier-family of functions
#include "dy4.h"
#include "fourier.h"

// just DFT function (no FFT yet)
void DFT(const std::vector<float> &x, std::vector<std::complex<float> > &Xf) {
	Xf.resize(x.size(), static_cast<std::complex<float> >(0, 0));
	for (long m = 0; m < Xf.size(); m++) {
		for (long k = 0; k < x.size(); k++) {
				std::complex<float> expval(0, -2*PI*(k*m) / x.size());
				Xf[m] += x[k] * std::exp(expval);
		}
	}
}

// function to compute the magnitude values in a complex vector
void computeVectorMagnitude(const std::vector<std::complex<float>> &Xf, std::vector<float> &Xmag)
{
	// only the positive frequencies
	Xmag.resize(Xf.size(), static_cast<float>(0));
  for (auto i = 0; i < Xf.size(); i++) {
    Xmag[i] = std::abs(Xf[i])/Xf.size();
  }
}

void estimatePSD(std::vector<double> &Samples, int &NFFT_in, double &Fs, std::vector<double> &freq, std::vector<double> &psd_est){

	float freq_bins = (float) NFFT_in;

	float df = Fs/freq_bins;

	freq.resize(round(freq_bins), static_cast<float>(0));

	for(long ii = 0 ; ii < freq.size() ; ii++){
		freq[ii] = df*ii;
	}

	std::vector<float> hann;

	hann.resize(round(freq_bins), static_cast<float>(0));

	for(long ii = 0 ; ii < hann.size() ; ii++){
		hann[ii] = pow(sin(ii*PI/freq_bins),2);
	}

	std::vector<float> psd_list;

	int no_segments = (int)floor( ((float) Samples.size()) / freq_bins);

	for(long jj = 0 ; jj < no_segments ; jj++){

		std::vector<float> windowed_samples;

		windowed_samples.resize(round(freq_bins), static_cast<float>(0));

		for(long ii = 0 ; ii < windowed_samples.size() ; ii++){
			windowed_samples[ii] = Samples[jj*freq_bins + ii] * hann[ii];
		}

		std::vector<std::complex<float> > Xf_i;

		DFT(windowed_samples,Xf_i);

		std::vector<std::complex<float> > Xf;

		Xf.resize((int)(freq_bins/(float)2), static_cast<float>(0) );
		for( long ii = 0 ; ii < Xf.size() ; ii++){
			Xf[ii] = Xf_i[ii];
		}

		std::vector<float> Xmag;

		computeVectorMagnitude(Xf,Xmag);

		std::vector<float> PSD_seg;
		PSD_seg.resize(Xmag.size(), static_cast<float>(0));

		for(long ii = 0 ; ii < PSD_seg.size() ; ii++){
			PSD_seg[ii] = 10 * log10( (2 * (1 / (Fs*freq_bins/2)) * pow(Xmag[ii],2)) );
			
		}

		psd_list.insert(psd_list.end() , PSD_seg.begin() , PSD_seg.end());
	}

	psd_est.resize((int) (freq_bins / (float) 2 ),  static_cast<float>(0));

	for(long k = 0 ; k < ((int) (freq_bins / (float) 2 )) ; k++){
		for(long l = 0 ; l < no_segments ; l++ ){
			psd_est[k] = psd_est[k] + psd_list[k + l * ((int) (freq_bins / (float) 2 ))];

		}
		psd_est[k] = psd_est[k] / (float) no_segments;
		

	}
}
