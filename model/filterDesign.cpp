#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <complex>
#include <cmath>

#define PI 3.14159265358979323846

// function for DFT (reused from previous experiment)
void DFT(const std::vector<float> &x, std::vector<std::complex<float> > &Xf)
{
	printf("This is the zize %d\n\n",x.size());
	Xf.resize(x.size(), static_cast<std::complex<float> >(0, 0));
	for (auto m = 0; m < Xf.size(); m++) {
		for (auto k = 0; k < x.size(); k++) {
				std::complex<float> expval(0, -2*PI*(k*m) / x.size());
				Xf[m] += x[k] * std::exp(expval);
		}
	}
}

void cDFT(const std::vector<float> &x, std::vector<std::complex<float> > &Xf)
{
	Xf.resize(x.size()/4, static_cast<std::complex<float> >(0, 0));

	std::vector<std::complex<float> > Xf1;
	std::vector<std::complex<float> > Xf2;
	std::vector<std::complex<float> > Xf3;
	std::vector<std::complex<float> > Xf4;
	
	Xf1.resize(x.size()/4, static_cast<std::complex<float> >(0, 0));
	Xf2.resize(x.size()/4, static_cast<std::complex<float> >(0, 0));
	Xf3.resize(x.size()/4, static_cast<std::complex<float> >(0, 0));
	Xf4.resize(x.size()/4, static_cast<std::complex<float> >(0, 0));

	for (auto m = 0; m < Xf1.size(); m++) {
		for (auto k = 0; k < x.size()/4; k++) {
				std::complex<float> expval(0, -2*PI*(k*m) / x.size());
				Xf1[m] += x[k] * std::exp(expval);
		}
	}

	for (auto m = 0; m < Xf2.size(); m++) {
		for (auto k = x.size()/4; k < 2*x.size()/4; k++) {
				std::complex<float> expval(0, -2*PI*((k - x.size())*m) / x.size());
				Xf2[m] += x[k] * std::exp(expval);
		}
	}

	for (auto m = 0; m < Xf3.size(); m++) {
		for (auto k = 2*x.size()/4; k < 3*x.size()/4; k++) {
				std::complex<float> expval(0, -2*PI*((k - 2*x.size())*m) / x.size());
				Xf3[m] += x[k] * std::exp(expval);
		}
	}

	for (auto m = 0; m < Xf4.size(); m++) {
		for (auto k = 3*x.size()/4; k < x.size(); k++) {
				std::complex<float> expval(0, -2*PI*((k - 3*x.size())*m) / x.size());
				Xf4[m] += x[k] * std::exp(expval);
		}
	}

	for(auto m = 0; m < Xf.size(); m++) {
		Xf[m] = (Xf1[m] + Xf2[m] + Xf3[m] + Xf4[m])/(std::complex<float>)4;
	}
}

// function to generate a sine with N samples per second over interval
void generateSin(std::vector<float> &t, std::vector<float> &x, float Fs, float interval, float frequency = 7.0, float amplitude = 5.0, float phase = 0.0)
{
	// we do NOT allocate memory space explicitly
	// for the time (t) vector and sample (x) vector
	t.resize(0); x.resize(0);
	float dt = 1/Fs;
	for (auto i = 0.0; i < interval; i += dt) {
		// vector size increases when pushing new elements into it
		t.push_back(i);
		x.push_back(amplitude*std::sin(2*PI*frequency*i+phase));
	}
}

// function to mix an array of sines
void mixSin(const std::vector<std::vector<float> > &sv, std::vector<float> &mixed)
{
	// assumes at least one sine passed
	// assumes all input sines are of the same size
	for (auto i = 0.0; i < sv[0].size(); i ++) {
		float mixval = 0.0;
		// note: sv.size() returns the number of sines (or rows in 2D repr)
		// sv[0].size() returns the number of samples in a sine (or cols in 2D repr)
		for (auto k = 0; k < sv.size(); k++)
			mixval += sv[k][i];
		mixed.push_back(mixval);
	}
}

// function to print a real vector (reused from previous experiment)
void printRealVector(const std::vector<float> &x)
{
	std::cout << "Printing float vector of size " << x.size() << "\n";
	for (auto i = 0; i < x.size(); i++)
		std::cout << x[i] << " ";
	std::cout << "\n";
}

// function to print a complex vector (reused from previous experiment)
void printComplexlVector(const std::vector<std::complex<float> > &X)
{
	std::cout << "Printing complex vector of size " << X.size() << "\n";
	for (auto i = 0; i < X.size(); i++)
		std::cout << X[i] << " ";
	std::cout << "\n";
}

// function to record data in a format to be read by GNU plot
// the arguments are VERY specific to this usage in this experiment
// we have the time vector (t), a vector of sines (sv),
// input samples (x, i.e., mixed sines for this experiment),
// output samples (y, filtered input samples)
// frequency vectors (both Xf and Yf)
// note: the reference code does NOT do filtering, hence y and Yf are zeros by default
void plotMixedSinesSpectrum(const std::vector<float> &t, const std::vector<std::vector<float> > &sv, const std::vector<float> &x, const std::vector<float> &y, const std::vector<std::complex<float> > &Xf, const std::vector<std::complex<float> > &Yf)
{
	// write data in text format to be parsed by gnuplot
	const std::string filename = "../data/example.dat";
	std::fstream fd;  // file descriptor
	fd.open(filename, std::ios::out);
	fd << "#\tindex\tsine(0)\tsine(1)\tsine(2)\tdata in\tdata out\tspectrum in\tspectrum out\n";
	for (auto i = 0; i < t.size(); i++) {
		fd << "\t " << i << "\t";
		for (auto k = 0; k < sv.size(); k++)
			fd << std::fixed << std::setprecision(3) << sv[k][i] << "\t ";
		fd << x[i] << "\t "<< y[i] << "\t\t ";
		fd << std::abs(Xf[i])/Xf.size() << "\t\t " << std::abs(Yf[i])/Yf.size() <<"\n";
	}
	std::cout << "Generated " << filename << " to be used by gnuplot\n";
	fd.close();
}

// function to compute the impulse response "h" based on the sinc function
// see pseudocode from previous lab that was implemented in Python
void impulseResponseLPF(float Fs, float Fc, unsigned short int num_taps, std::vector<float> &h)
{
	// allocate memory for the impulse response
	h.resize(num_taps, 0.0);

	// the rest of the code in this function is to be completed by you
	// based on your understanding and the Python code from the first lab
}

// function to compute the filtered output "y" by doing the convolution
// of the input data "x" with the impulse response "h"; this is based on
// your Python code from the take-home exercise from the previous lab
void convolveFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h)
{
	// allocate memory for the output (filtered) data
	y.resize(x.size()+h.size()-1, 0.0);

	// the rest of the code in this function is to be completed by you
	// based on your understanding and the Python code from the first lab
}

int main()
{

	float Fs = 1024.0;                    // samples per second
	float interval = 1.0;                 // number of seconds
	unsigned short int num_taps = 101;    // number of filter taps
	float Fc = 20.0;                      // cutoff frequency (in Hz)

	// declare a vector of vectors for multiple sines
	std::vector<std::vector<float> > sv;
	// declare time and sine vectors
	std::vector<float> t, sine;
	// note: there is no explicit memory allocation through vector resizing
	// vector memory space will increase via the push_back method

	// generate and store the first tone
	// check the function to understand the order of arguments
	generateSin(t, sine, Fs, interval, 10.0, 5.0, 0.0);
	sv.push_back(sine);
	// generate and store the second tone
	generateSin(t, sine, Fs, interval, 40.0, 2.0, 0.0);
	sv.push_back(sine);
	// generate and store the third tone
	generateSin(t, sine, Fs, interval, 50.0, 3.0, 0.0);
	sv.push_back(sine);

	// declare the mixed sine vector and mix the three tones
	std::vector<float> x;
	mixSin(sv, x);
	// printRealVector(x);

	// declare a vector of complex values for DFT; no memory is allocated for it
	std::vector<std::complex<float> > Xf;
	cDFT(x, Xf);
	// printComplexlVector(Xf);

	// generate the impulse response h
	// convolve it with the input data x
	// in order to produce the output data y
	std::vector<float> h;              // impulse response
	impulseResponseLPF(Fs, Fc, num_taps, h);
	std::vector<float> y;              // filter out
	convolveFIR(y, x, h);

	// compute DFT of the filtered data
	std::vector<std::complex<float> > Yf;
	cDFT(y, Yf);

	printf("Yo mans are guessing the zize %d\n\n",Yf.size());

	// prepare the data for gnuplot
	plotMixedSinesSpectrum(t, sv, x, y, Xf, Yf);

	// naturally, you can comment the line below once you are comfortable to run gnuplot
	std::cout << "Run: gnuplot -e 'set terminal png size 1024,768' example.gnuplot > ../data/example.png\n";

	return 0;
}
