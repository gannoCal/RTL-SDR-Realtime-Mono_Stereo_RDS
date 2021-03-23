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
void fmDemodArctanBlock(std::vector<double> &fm_demod,std::vector<double> &I, std::vector<double> &Q,std::vector<double> &prev_phase);
#endif // DY4_FILTER_H
