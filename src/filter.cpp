/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"
#include <math.h>
#define PI 3.14159265358979323846
// function to compute the impulse response "h" based on the sinc function
void impulseResponseLPF(double Fs, double Fc, unsigned short int num_taps, std::vector<double> &h, double decim)
{
	
	double cutoff = Fc/((Fs/decim)/2);
	// allocate memory for the impulse response
	h.resize(num_taps, 0.0);
	for(auto i = 0 ; i < num_taps ; i++){
		if(i == (num_taps-1)/2){
			h[i] = cutoff;
		}else{
			h[i] = cutoff * sin( PI * cutoff *( i-(num_taps-1)/2 ) ) / ( PI * cutoff *( i-(num_taps-1)/2 ) );
		}
		h[i] = h[i] * (sin(i * PI / num_taps)*sin(i * PI / num_taps));
		//printf("h[%d] = %f\n",i,h[i]);
	}
	
}


void impulseResponseBPF(double Fs, double Fb,double Fe, unsigned short int num_taps, std::vector<double> &h, double decim)
{
	
	double center = ((Fe+Fb)/2.0)/((Fs/decim)/2.0);
    double pass = (Fe-Fb)/((Fs/decim)/2.0);
	// allocate memory for the impulse response
	h.resize(num_taps, 0.0);
	for(auto i = 0 ; i < num_taps ; i++){
		if(i == (num_taps-1)/2){
			h[i] = pass;
		}else{
			h[i] = pass * sin( PI * (pass/2) *( i-(num_taps-1)/2 ) ) / ( PI * (pass/2) *( i-(num_taps-1)/2 ) );
		}
        h[i] = h[i] * cos(i*PI*center);
		h[i] = h[i] * (sin(i * PI / num_taps)*sin(i * PI / num_taps));
		//printf("h[%d] = %f\n",i,h[i]);
	}
	
}


// function to compute the filtered output "y" by doing the convolution
// of the input data "x" with the impulse response "h"
/*void convolveFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h)
{
	// allocate memory for the output (filtered) data
	y.resize(x.size()+h.size()-1, 0.0);

	// the rest of the code in this function is to be completed by you
	// based on your understanding and the Python code from the first lab
	auto max_size = x.size();
	if (h.size() > max_size)
	{
		max_size = h.size();
	}
	for (auto n = 0; n < y.size(); n++)
	{
		for (auto m = 0; m < h.size(); m++)
		{
			if ((n-m) >= 0 || (n-m) < max_size)
			{
				y[n] += x[n-m] * h[m];
			}
		}
	}
}*/


void convolveFIR_N_dec(const int step_size, std::vector<double> &y, const std::vector<double> &x, const std::vector<double> &h, std::vector<double> &state )
{
	auto max_size = x.size();
	if (h.size() > max_size)
	{
		max_size = h.size();
	}
	long special = 0;
	for (auto n = 0; n < y.size(); n++)
	{
		special = 0;
		y[n] = 0;
		for (auto m = 0; m < h.size(); m++)
		{
			if ((step_size*n-m) >= 0 && (step_size*n-m) < max_size)
			{
				y[n] += x[step_size*n-m] * h[m];
			}else if((step_size*n-m) < 0 && state.size() > 0){
				y[n] += state[state.size() - 1 - special] * h[m];
				special++;
			}
			

		}
	}
	for(auto ii = 0 ; ii < state.size(); ii++){
		state[ii] = x[(x.size()) - state.size() + ii];
	}
}

void convolveFIR_N_dec_RDS(const int step_size, std::vector<double> &y, const std::vector<double> &x, const std::vector<double> &h, std::vector<double> &state )
{
    y.resize(trunc(x.size()/step_size));
    auto max_size = x.size();
	if (h.size() > max_size)
	{
		max_size = h.size();
	}
	long special = 0;
	for (auto n = 0; n < y.size(); n++)
	{
		special = 0;
		y[n] = 0;
		for (auto m = 0; m < h.size(); m++)
		{
			if ((step_size*n-m) >= 0 && (step_size*n-m) < max_size)
			{
				y[n] += x[step_size*n-m] * h[m];
			}else if((step_size*n-m) < 0 && state.size() > 0){
				y[n] += state[state.size() - 1 - special] * h[m];
				special++;
			}
			

		}
	}
	for(auto ii = 0 ; ii < state.size(); ii++){
		state[ii] = x[(x.size()) - state.size() + ii];
	}
}

void convolve_UPSAMPLE_N_dec_New(int step_size,int upsample_size, std::vector<double> &y, const std::vector<double> &x, const std::vector<double> &h, std::vector<double> &state)
{
    y.resize(trunc(upsample_size*x.size()/step_size));
    auto max_size = x.size();
	if (h.size() > max_size)
	{
		max_size = h.size();
	}
    int r0;
    int limit;
    int u = upsample_size;
    int d = step_size;
    int test = 0;
    int testcnt = 0;
	for (auto n = 0; n < y.size(); n++){
        r0 = h.size()%upsample_size;
        limit = (n%upsample_size < r0) ? trunc(((double)h.size())/ (double)upsample_size)+1 : trunc(((double)h.size())/ (double)upsample_size);
        y[n] = 0;
        test = 0;
        for(auto k = 0 ; k < limit ; k++){
            if (trunc(n*d/u) - k >= 0 && trunc(n*d/u) - k < max_size){
                y[n] += h[ (n*d)%u + u*k] * x[trunc(n*d/u) - k];
                test = 1;
            }else if(trunc(n*d/u) - k < 0 && state.size() > 0){
                y[n] += h[ (n*d)%u + u*k] * state[state.size() + trunc(n*d/u) - k];
                test = 1;
            }
        }
        if(test == 1)
            testcnt++;
        
        
    }
	for(auto ii = 0 ; ii < state.size(); ii++){
		state[ii] = x[(x.size()) - state.size() + ii];
	}

}


void fmDemodArctanBlock(std::vector<double> &fm_demod,std::vector<double> &I, std::vector<double> &Q,std::vector<double> &prev_phase){
	fm_demod.resize(I.size(), 0.0);
	double thetadelta = 0, a, b, c, current_phase;
	for(auto n = 0; n < I.size(); n++){
		a = b =c = current_phase = 0;
		if(n == 0){
			a = I[n]*(Q[n]-prev_phase[0]);		//prev phase is never being stored
			b = Q[n]*(I[n]-prev_phase[1]);
			c = (I[n]*I[n] + Q[n]*Q[n]);
			thetadelta = (a-b)/c;
		}
		else{
			a = I[n]*(Q[n]-Q[n-1]);
			b = Q[n]*(I[n]-I[n-1]);
			c = (I[n]*I[n] + Q[n]*Q[n]);

			thetadelta = (a-b)/c;
		}
        
		if(!std::isnan(thetadelta)){
		fm_demod[n] = thetadelta;
		}else{
		fm_demod[n] = (a-b)*2;	
		}
	}
	prev_phase.resize(2);
	prev_phase[0] = Q[Q.size() - 1];
	prev_phase[1] = I[I.size() - 1];
}


void convolve_UPSAMPLE_N_dec(int step_size,int upsample_size, std::vector<double> &y, const std::vector<double> &x, const std::vector<double> &h, std::vector<double> &state)
{
	auto max_size = x.size();
	if (h.size() > max_size)
	{
		max_size = h.size();
	}
	long special = 0;
	for(auto phase = 0; phase < upsample_size; phase++){

			for (auto n = phase; n < y.size(); n = n + upsample_size)
			{
				int x_count = n;
				special = 0;
				y[n] = 0;
				for (auto m = 0; m < h.size(); m = m + 1)
				{
					if (((step_size%upsample_size)*n-m) >= 0 && ((step_size%upsample_size)*n-m) < max_size)
						{
							y[n] += x[(step_size%upsample_size)*n-m] * h[m];
						}else if(((step_size%upsample_size)*n-m) < 0 && state.size() > 0){
							y[n] += state[state.size() - 1 - special] * h[m];
							special++;
						}

					}
			}
	}

	for(auto ii = 0 ; ii < state.size(); ii++){
		state[ii] = x[(x.size()) - state.size() + ii];
	}

}

void impulseResponseRootRaisedCosine(double Fs,int N_taps,std::vector<double> &impulseResponseRRC){
// duation for each symbol - do NOT be changed for RDS!
	double T_symbol = 1/2375.0;

	// roll-off factor (greater than 0 and smaller than 1)
	double beta = 0.90;

	// the RRC inpulse response that will be computed in this function
	impulseResponseRRC.resize(N_taps,0);

	for (auto k = 0 ; k <(N_taps) ; k++){
		double t = (double)(((double)k-(double)N_taps/2.0))/Fs;
		// we ignore the 1/T_symbol scale factor
		if (t == 0.0){ 
            impulseResponseRRC[k] = 1.0 + beta*((4.0/PI)-1.0);
            }
		else if (t == -T_symbol/(4.0*beta) || t == T_symbol/(4.0*beta)){
			impulseResponseRRC[k] = (beta/sqrt(2.0))*(((1.0+2.0/PI)* 
					(sin(PI/(4.0*beta)))) + ((1.0-2.0/PI)*(cos(PI/(4.0*beta)))));
                    }
		else{ impulseResponseRRC[k] = (sin(PI*t*(1.0-beta)/T_symbol) +  
					4.0*beta*(t/T_symbol)*cos(PI*t*(1.0+beta)/T_symbol))/ 
					(PI*t*(1.0-(4.0*beta*t/T_symbol)*(4.0*beta*t/T_symbol))/T_symbol);
                    }
    }
	
}

void impulseResponseRootRaisedCosine2(double Fs, int N_taps, std::vector<double> &impulseResponseRRC)
{
	double T_s = 1.0/2375.0;
	float beta = 0.9;
	impulseResponseRRC.resize(N_taps);
	auto t = 0;
	for(int k = 0; k < N_taps;k++)
	{
		t = float((k-N_taps/2))/Fs;

		if(t == 0.0)
		{
			impulseResponseRRC[k] = 1.0 + beta*((4/PI)-1);
		}
		else if(t == -T_s/(4*beta) || t == T_s/(4*beta))
		{
			impulseResponseRRC[k] = (beta/sqrt(2))*(((1+2/PI)*(sin(PI/(4*beta)))) + ((1-2/PI)*(cos(PI/(4*beta)))));
		}
		else
		{
			impulseResponseRRC[k] = (sin(PI*t*(1-beta)/T_s) +  4*beta*(t/T_s)*cos(PI*t*(1+beta)/T_s))/ (PI*t*(1-(4*beta*t/T_s)*(4*beta*t/T_s))/T_s);
		}
	}
}

void process_MBA(int &new_bit, std::vector<int> &MBA, int &previous_match, int &is_nSync, 
int &nSync_Hit_counter, int &nSync_Flop_counter, int &allowed_nSync_Flops, 
std::vector<std::vector<int>> &found_array,
std::vector< std::vector<int> > &parityArray, std::vector<int> &syndrome){
    MBA.insert(MBA.begin(), new_bit);
    MBA.erase(MBA.end()-1);

    if(is_nSync == 0){
        Gfield_mult_reverse(MBA,parityArray,syndrome);
        previous_match = 0;
        int i=0;
        if (syndrome[i+0] == 1 &&
        syndrome[i+1] == 1 &&
        syndrome[i+2] == 1 &&
        syndrome[i+3] == 1 &&
        syndrome[i+4] == 0 &&
        syndrome[i+5] == 1 &&
        syndrome[i+6] == 1 &&
        syndrome[i+7] == 0 &&
        syndrome[i+8] == 0 &&
        syndrome[i+9] == 0) {
            //Do append (Push????)
            std::vector<int> ele(2);
            ele[0] = (int)'A';
            ele[1] = 0;
            found_array.insert(found_array.end(),ele);
            previous_match=(int)'D';
            //std::cerr << ("Found code A") << "\n";
        }
        else if (syndrome[i+0] == 1 &&
        syndrome[i+1] == 1 &&
        syndrome[i+2] == 1 &&
        syndrome[i+3] == 1 &&
        syndrome[i+4] == 0 &&
        syndrome[i+5] == 1 &&
        syndrome[i+6] == 0 &&
        syndrome[i+7] == 1 &&
        syndrome[i+8] == 0 &&
        syndrome[i+9] == 0) {
            //Do append
            std::vector<int> ele(2);
            ele[0] = (int)'B';
            ele[1] = 0;
            found_array.insert(found_array.end(),ele);
            previous_match=(int)'A';
            //std::cerr << ("Found code B")<< "\n";
        }
        else if (syndrome[i+0] == 1 &&
        syndrome[i+1] == 0 &&
        syndrome[i+2] == 0 &&
        syndrome[i+3] == 1 &&
        syndrome[i+4] == 0 &&
        syndrome[i+5] == 1 &&
        syndrome[i+6] == 1 &&
        syndrome[i+7] == 1 &&
        syndrome[i+8] == 0 &&
        syndrome[i+9] == 0) {
            //Do append
            std::vector<int> ele(2);
            ele[0] = (int)'C';
            ele[1] = 0;
            found_array.insert(found_array.end(),ele);
            previous_match=(int)'B';
            //std::cerr << ("Found code C")<< "\n";
        }
        else if (syndrome[i+0] == 1 &&
        syndrome[i+1] == 1 &&
        syndrome[i+2] == 1 &&
        syndrome[i+3] == 1 &&
        syndrome[i+4] == 0 &&
        syndrome[i+5] == 0 &&
        syndrome[i+6] == 1 &&
        syndrome[i+7] == 1 &&
        syndrome[i+8] == 0 &&
        syndrome[i+9] == 0) {
            //Do append
            std::vector<int> ele(2);
            ele[0] = (int)'C';
            ele[1] = 0;
            found_array.insert(found_array.end(),ele);
            previous_match=(int)'B';
            //std::cerr << ("Found code C")<< "\n";
        }
        else if (syndrome[i+0] == 1 &&
        syndrome[i+1] == 0 &&
        syndrome[i+2] == 0 &&
        syndrome[i+3] == 1 &&
        syndrome[i+4] == 0 &&
        syndrome[i+5] == 1 &&
        syndrome[i+6] == 1 &&
        syndrome[i+7] == 0 &&
        syndrome[i+8] == 0 &&
        syndrome[i+9] == 0) {
            //Do append
            std::vector<int> ele(2);
            ele[0] = (int)'D';
            ele[1] = 0;
            found_array.insert(found_array.end(),ele);
            previous_match=(int)'C';
            //std::cerr << ("Found code D")<< "\n";
        }

        std::vector<int> del_array;
        for (auto i=0; i<found_array.size()/*Size of first dimension of found array*/; i++){
            if(found_array[i][0] == previous_match && found_array[i][1] == 26){
                is_nSync = 1;
                if(previous_match == (int)'A')
                    previous_match = (int)'B';
                else if(previous_match == (int)'B')
                    previous_match = (int)'C';
                else if(previous_match == (int)'C')
                    previous_match = (int)'D';
                else if(previous_match == (int)'D')
                    previous_match = (int)'A';
                found_array.clear();
                std::cerr << ("Synchronized...")<< "\n";
                break;
            }else if (found_array[i][0] != previous_match && found_array[i][1] >= 26)
                del_array.insert(del_array.begin(),i);
            else
                found_array[i][1] = found_array[i][1] + 1;
        }

        if(is_nSync == 0){  
            for( auto i = 0; i < del_array.size() ; i++){
                found_array.erase(found_array.begin() + del_array[i]);
            }
        }
    }else{
        if(nSync_Hit_counter == 26){
            int codeFound = 0;
            Gfield_mult_reverse(MBA,parityArray,syndrome);

            int i=0;
            if (syndrome[i+0] == 1 &&
            syndrome[i+1] == 1 &&
            syndrome[i+2] == 1 &&
            syndrome[i+3] == 1 &&
            syndrome[i+4] == 0 &&
            syndrome[i+5] == 1 &&
            syndrome[i+6] == 1 &&
            syndrome[i+7] == 0 &&
            syndrome[i+8] == 0 &&
            syndrome[i+9] == 0) {
                codeFound = 1;
                if(previous_match != (int)'D'){
                    if(nSync_Flop_counter >= allowed_nSync_Flops){
                        std::cerr << ("De-synchronized...") << "\n";
                        is_nSync = 0;
                        nSync_Hit_counter = 0;
                        nSync_Flop_counter = 0;
                    }else{
                        std::cerr << ("Expected Code D") << "\n";
                        nSync_Flop_counter = nSync_Flop_counter + 1;
                    }

                }else{
                    std::cerr << ("Found code A") << "\n";
                    nSync_Flop_counter = 0;
                }
                if(previous_match == (int)'A')
                    previous_match = (int)'B';
                else if(previous_match == (int)'B')
                    previous_match = (int)'C';
                else if(previous_match == (int)'C')
                    previous_match = (int)'D';
                else if(previous_match == (int)'D')
                    previous_match = (int)'A';
                
            }
            else if (syndrome[i+0] == 1 &&
            syndrome[i+1] == 1 &&
            syndrome[i+2] == 1 &&
            syndrome[i+3] == 1 &&
            syndrome[i+4] == 0 &&
            syndrome[i+5] == 1 &&
            syndrome[i+6] == 0 &&
            syndrome[i+7] == 1 &&
            syndrome[i+8] == 0 &&
            syndrome[i+9] == 0) {
                codeFound = 1;
                if(previous_match != (int)'A'){
                  if(nSync_Flop_counter >= allowed_nSync_Flops){
                        std::cerr << ("De-synchronized...") << "\n";
                        is_nSync = 0;
                        nSync_Hit_counter = 0;
                        nSync_Flop_counter = 0;
                    }else{
                        std::cerr << ("Expected Code A") << "\n";
                        nSync_Flop_counter = nSync_Flop_counter + 1;
                    }

                }else{
                    std::cerr << ("Found code B") << "\n";
                    nSync_Flop_counter = 0;
                }
                if(previous_match == (int)'A')
                    previous_match = (int)'B';
                else if(previous_match == (int)'B')
                    previous_match = (int)'C';
                else if(previous_match == (int)'C')
                    previous_match = (int)'D';
                else if(previous_match == (int)'D')
                    previous_match = (int)'A';

            }
            else if (syndrome[i+0] == 1 &&
            syndrome[i+1] == 0 &&
            syndrome[i+2] == 0 &&
            syndrome[i+3] == 1 &&
            syndrome[i+4] == 0 &&
            syndrome[i+5] == 1 &&
            syndrome[i+6] == 1 &&
            syndrome[i+7] == 1 &&
            syndrome[i+8] == 0 &&
            syndrome[i+9] == 0) {
                codeFound = 1;
                if(previous_match != (int)'B'){
                    if(nSync_Flop_counter >= allowed_nSync_Flops){
                        std::cerr << ("De-synchronized...") << "\n";
                        is_nSync = 0;
                        nSync_Hit_counter = 0;
                        nSync_Flop_counter = 0;
                    }else{
                        std::cerr << ("Expected Code B") << "\n";
                        nSync_Flop_counter = nSync_Flop_counter + 1;
                    }

                }else{
                    std::cerr << ("Found code C") << "\n";
                    nSync_Flop_counter = 0;
                }
                if(previous_match == (int)'A')
                    previous_match = (int)'B';
                else if(previous_match == (int)'B')
                    previous_match = (int)'C';
                else if(previous_match == (int)'C')
                    previous_match = (int)'D';
                else if(previous_match == (int)'D')
                    previous_match = (int)'A';

            }
            else if (syndrome[i+0] == 1 &&
            syndrome[i+1] == 1 &&
            syndrome[i+2] == 1 &&
            syndrome[i+3] == 1 &&
            syndrome[i+4] == 0 &&
            syndrome[i+5] == 0 &&
            syndrome[i+6] == 1 &&
            syndrome[i+7] == 1 &&
            syndrome[i+8] == 0 &&
            syndrome[i+9] == 0) {
                codeFound = 1;
                if(previous_match != (int)'B'){
                    if(nSync_Flop_counter >= allowed_nSync_Flops){
                        std::cerr << ("De-synchronized...") << "\n";
                        is_nSync = 0;
                        nSync_Hit_counter = 0;
                        nSync_Flop_counter = 0;
                    }else{
                        std::cerr << ("Expected Code B") << "\n";
                        nSync_Flop_counter = nSync_Flop_counter + 1;
                    }

                }else{
                    std::cerr << ("Found code C'") << "\n";
                    nSync_Flop_counter = 0;
                }
                if(previous_match == (int)'A')
                    previous_match = (int)'B';
                else if(previous_match == (int)'B')
                    previous_match = (int)'C';
                else if(previous_match == (int)'C')
                    previous_match = (int)'D';
                else if(previous_match == (int)'D')
                    previous_match = (int)'A';

            }
            else if (syndrome[i+0] == 1 &&
            syndrome[i+1] == 0 &&
            syndrome[i+2] == 0 &&
            syndrome[i+3] == 1 &&
            syndrome[i+4] == 0 &&
            syndrome[i+5] == 1 &&
            syndrome[i+6] == 1 &&
            syndrome[i+7] == 0 &&
            syndrome[i+8] == 0 &&
            syndrome[i+9] == 0) {
                codeFound = 1;
                if(previous_match != (int)'C'){
                    if(nSync_Flop_counter >= allowed_nSync_Flops){
                        std::cerr << ("De-synchronized...") << "\n";
                        is_nSync = 0;
                        nSync_Hit_counter = 0;
                        nSync_Flop_counter = 0;
                    }else{
                        std::cerr << ("Expected Code C") << "\n";
                        nSync_Flop_counter = nSync_Flop_counter + 1;
                    }

                }else{
                    std::cerr << ("Found code D") << "\n";
                    nSync_Flop_counter = 0;
                }
                if(previous_match == (int)'A')
                    previous_match = (int)'B';
                else if(previous_match == (int)'B')
                    previous_match = (int)'C';
                else if(previous_match == (int)'C')
                    previous_match = (int)'D';
                else if(previous_match == (int)'D')
                    previous_match = (int)'A';

            }

            if(codeFound == 0){
                if(nSync_Flop_counter >= allowed_nSync_Flops){
                        std::cerr << ("De-synchronized...") << "\n";
                        is_nSync = 0;
                        nSync_Hit_counter = 0;
                        nSync_Flop_counter = 0;
                }else{
                    nSync_Hit_counter = 0;
                    nSync_Flop_counter = nSync_Flop_counter + 1;

                    if(previous_match == (int)'A')
                        previous_match = (int)'B';
                    else if(previous_match == (int)'B')
                        previous_match = (int)'C';
                    else if(previous_match == (int)'C')
                        previous_match = (int)'D';
                    else if(previous_match == (int)'D')
                        previous_match = (int)'A';
                    std::cerr << ("Expected Code ") << (char)previous_match << "\n";
                }
            }
            nSync_Hit_counter = 0;

        }

        nSync_Hit_counter = nSync_Hit_counter + 1;
    }
}

void Gfield_mult_reverse(std::vector<int> &MBA, std::vector< std::vector<int> > &parityArray, std::vector<int> &syndrome){
    syndrome.clear();
    syndrome.resize(10,0.0);
    int working_column;
    for(auto i = 0 ; i < 10 ; i++){
        working_column = (MBA[25] ==1 && parityArray[0][i] == 1) ? 1 : 0;
        for(auto j = 1; j <26 ; j++){
            working_column = (((MBA[25-j] == 1 && parityArray[j][i] == 1) ? 1 : 0) != working_column ) ?  1 : 0;
        }
        syndrome[i] = working_column;
    }
}

void Gfield_mult(std::vector<int> &MBA, std::vector< std::vector<int> > &parityArray, std::vector<int> &syndrome){
    syndrome.clear();
    syndrome.resize(10,0.0);
    int working_column;
    for(auto i = 0 ; i < 10 ; i++){
        working_column = (MBA[0] ==1 && parityArray[0][i] == 1) ? 1 : 0;
        for(auto j = 1; j <26 ; j++){
            working_column = (((MBA[j] == 1 && parityArray[j][i] == 1) ? 1 : 0) != working_column ) ?  1 : 0;
        }
        syndrome[i] = working_column;
    }
}
