/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_FILTER_H
#define DY4_FILTER_H

// add headers as needed
#include <iostream>
#include <vector>

// declaration of a function prototypes
void impulseResponseLPF(double, double, unsigned short int, std::vector<double> &,double);
void impulseResponseBPF(double, double,double, unsigned short int, std::vector<double> &,double);
//void convolveFIR(std::vector<float> &, const std::vector<float> &, const std::vector<float> &);
void convolveFIR_N_dec(const int, std::vector<double> &, const std::vector<double> &, const std::vector<double> &, std::vector<double> &);
void convolveFIR_N_dec_RDS(const int, std::vector<double> &, const std::vector<double> &, const std::vector<double> &, std::vector<double> &);
void fmDemodArctanBlock(std::vector<double> &,std::vector<double> &, std::vector<double> &,std::vector<double> &);
void convolve_UPSAMPLE_N_dec(int,int, std::vector<double> &, const std::vector<double> &, const std::vector<double> &, std::vector<double> &);
void convolve_UPSAMPLE_N_dec_New(int step_size,int upsample_size, std::vector<double> &y, const std::vector<double> &x, const std::vector<double> &h, std::vector<double> &state);
void impulseResponseRootRaisedCosine(double ,int ,std::vector<double>&);
void impulseResponseRootRaisedCosine2(double ,int ,std::vector<double>&);
void process_MBA(int &new_bit, std::vector<int> &MBA, int &previous_match, int &is_nSync, 
int &nSync_Hit_counter, int &nSync_Flop_counter, int &allowed_nSync_Flops, 
std::vector<std::vector<int>> &found_array,
std::vector< std::vector<int> > &parityArray, std::vector<int> &syndrome);
void Gfield_mult_reverse(std::vector<int> &MBA, std::vector< std::vector<int> > &parityArray, std::vector<int> &syndrome);
void Gfield_mult(std::vector<int> &MBA, std::vector< std::vector<int> > &parityArray, std::vector<int> &syndrome);

#endif // DY4_FILTER_H
