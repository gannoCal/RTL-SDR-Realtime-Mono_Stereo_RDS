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
void resampler(const int , const int , std::vector<double> &, const std::vector<double> &, const std::vector<double> &, std::vector<double> &);
void CDR(std::vector<double>, int , int );
void Manchester_and_differntial(std::vector<double> , int , std::vector<double> , int);
impulseResponseRootRaisedCosine(double , int , std::vector<double> );

#endif // DY4_FILTER_H
