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
void DFT(const std::vector<float> &x, std::vector<std::complex<float>> &Xf) {
	Xf.resize(x.size(), static_cast<std::complex<float>>(0, 0));
	for (auto m = 0; m < Xf.size(); m++) {
		for (auto k = 0; k < x.size(); k++) {
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

void estimatePSD(std::vector<float> &samples, const float Fs,std::vector<float> &freq,std::vector<float> &psd_est){
  int freq_bins = NFFT;
  int df = Fs/freq_bins;
	freq.resize(Fs/(2*df));
  for(auto i = 0; i <freq.size();i++ ){
    freq[i] = i*df/1000;
  }

  std::vector <float> hann (freq_bins);
  for(auto i = 0; i <hann.size();i++ ){
    hann[i] = pow(sin(i*PI/freq_bins),2);
  }

  std::vector <float> psd_list (1);

  int no_segments = floor(samples.size()/float(freq_bins));
// K LOOP
  for (auto k = 0; k < no_segments; k++){
    std::vector <float> windowed_samples (freq_bins); // fix this

		for(auto i = 0; i < windowed_samples.size(); i++){windowed_samples[i] = samples[k*freq_bins + i]*hann[i];}


     std::vector<std::complex<float>> Xf (freq_bins);
    DFT(windowed_samples,Xf);
		std::cout << Xf[4];
		std::vector <float> psd_seg (freq_bins/2);
		for(auto i = 0; i < psd_seg.size(); i++){psd_seg[i] = 1/(Fs*freq_bins/2) * pow(abs(Xf[i]),2);}
		for(auto i = 0; i < psd_seg.size(); i++){psd_seg[i] = 2*psd_seg[i];}
		for (auto k = 0; k < int (psd_seg.size()); k++){psd_seg[k] = 10*log10(psd_seg[k]);}

		psd_list.insert(psd_list.end(), psd_seg.begin(), psd_seg.end());

  }// END K LOOP


	psd_est.resize(int(freq_bins/2));

for(auto k = 0; k <int(freq_bins/2); k++ ){
  for(auto l = 0; l <no_segments; l++ ){
    psd_est[k] = psd_est[k] + psd_list[k + l*int(freq_bins/2)];
  }
  psd_est[k] = psd_est[k] / no_segments;
}

}
