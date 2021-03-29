/*
Comp Eng 3DY4 (Computer Systems Integration Project)
Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_RDS_H
#define DY4_RDS_H

// add headers as needed
#include <iostream>
#include <vector>
#include "dy4.h"
#include "filter.h"
#include <unistd.h>
#include <math.h>
#include <valarray>
// declaration of a function prototypes
void resampler(const int , const int , std::vector<double> &, const std::vector<double> &, const std::vector<double> &, std::vector<double> &);
void CDR(std::vector<double>, int , int );
void Manchester_and_differntial(std::vector<double> , int , std::vector<double>);
void impulseResponseRootRaisedCosine(double , int , std::vector<double> );

#endif // DY4_FILTER_H
