/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "logfunc.h"
#include <tgmath.h>
#define PI 3.14159265358979323846

/*

IN THE VAR prevstate, each postion strores the following

prevstate[0] = integrator;
prevstate[1] = phaseEst;
prevstate[2] = feedbackI;
prevstate[3] = feedbackQ;
prevstate[4] = ncoOut[ncoOut - 1];

TO BE USED IN THE NEXT STATE

*/

// function to generate a vector whose value is equal to its index
// this is useful when plotting a vector because we use the index on the X axis
void fmPll(std::vector<double> &pllIn, const double &freq,const double &Fs,const double &ncoScale, const double &phaseAdjust,const double &normBandwidth,std::vector<double> &ncoOut, std::vector<double> &prevstate) {
    double Cp = 2.666;
    double Ci = 3.555;
    double Kp = (normBandwidth)*Cp;
    double Ki = (normBandwidth*normBandwidth)*Ci;
    ncoOut.resize(pllIn.size(),0);
// implement state saving later begin
    double integrator = prevstate[0];
    double phaseEst = prevstate[1];
    double feedbackI = prevstate[2];
    double feedbackQ = prevstate[3];
    ncoOut[0] = prevstate[4];
    double trigOffset = prevstate[5];
    // implement state saving later end
    double trigArg;

    double errorI;
    double errorQ;
    double errorD;

    for(auto k = 0; k < pllIn.size(); k++){
        errorI = pllIn[k] * (feedbackI);
        errorQ = pllIn[k] * (-1*feedbackQ);

        errorD =  atan2(errorQ,errorI);

        integrator = integrator + Ki*errorD;

        phaseEst = phaseEst + Kp*errorD + integrator;

        trigArg = 2*PI*(freq/Fs)*(trigOffset+k) + phaseEst;
        feedbackI = cos(trigArg);
        feedbackQ = sin(trigArg);
        ncoOut[k] = cos(trigArg*ncoScale + phaseAdjust);
    }

    prevstate[0] = integrator;
    prevstate[1] = phaseEst;
    prevstate[2] = feedbackI;
    prevstate[3] = feedbackQ;
    prevstate[4] = ncoOut[ncoOut.size() - 1];
    prevstate[5] = (trigOffset + pllIn.size());
}


void fmPll_RDS(std::vector<double> &pllIn, const double &freq,const double &Fs,const double &ncoScale, const double &phaseAdjust,const double &normBandwidth,std::vector<double> &ncoOut, std::vector<double> &prevstate) {
    double Cp = 2.666;
    double Ci = 3.555;
    double Kp = (normBandwidth)*Cp;
    double Ki = (normBandwidth*normBandwidth)*Ci;
    ncoOut.resize(pllIn.size()+1,0);
// implement state saving later begin
    double integrator = prevstate[0];
    double phaseEst = prevstate[1];
    double feedbackI = prevstate[2];
    double feedbackQ = prevstate[3];
    ncoOut[0] = prevstate[4];
    double trigOffset = prevstate[5];
    // implement state saving later end
    double trigArg;

    double errorI;
    double errorQ;
    double errorD;

    for(auto k = 0; k < pllIn.size(); k++){
        errorI = pllIn[k] * (feedbackI);
        errorQ = pllIn[k] * (-1*feedbackQ);

        errorD =  atan2(errorQ,errorI);

        integrator = integrator + Ki*errorD;

        phaseEst = phaseEst + Kp*errorD + integrator;

        trigArg = 2*PI*(freq/Fs)*(trigOffset+k+1) + phaseEst;
        feedbackI = cos(trigArg);
        feedbackQ = sin(trigArg);
        ncoOut[k+1] = cos(trigArg*ncoScale + phaseAdjust);
    }

    prevstate[0] = integrator;
    prevstate[1] = phaseEst;
    prevstate[2] = feedbackI;
    prevstate[3] = feedbackQ;
    prevstate[4] = ncoOut[ncoOut.size() - 1];
    prevstate[5] = (trigOffset + pllIn.size());
}
