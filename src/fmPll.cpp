/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "logfunc.h"
#include <cmath>
#include "fmPll.h"

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
void fmPll(std::vector<float> &pllIn, const float &freq,const float &Fs,const float &ncoScale, const float &phaseAdjust,const float &normBandwidth,std::vector<float> &ncoOut, std::vector<float> &prevstate) {
	float Cp = 2.666;
	float Ci = 3.555;
	float Kp = (normBandwidth)*Cp;
	float Ki = (normBandwidth*normBandwidth)*Ci;
	ncoOut.resize(pllIn.size()+1,0);
// implement state saving later begin
	float integrator = prevstate[0];
	float phaseEst = prevstate[1];
	float feedbackI = prevstate[2];
	float feedbackQ = prevstate[3];
	ncoOut[0] = prevstate[4];
	// implement state saving later end
	for(auto k = 0; k < pllIn.size(); k++){
		float errorI = pllIn[k] * (feedbackI);
		float errorQ = pllIn[k] * (-1*feedbackQ);

		float errorD =  atan2(errorQ,errorI);

		integrator = integrator + Ki*errorD;

		phaseEst = phaseEst + Kp*errorD + integrator;

		double trigArg = 2*PI*(freq/Fs)*(k+1) + phaseEst;
		feedbackI = cos(trigArg);
		feedbackQ = sin(trigArg);
		ncoOut[k+1] = cos(trigArg*ncoScale + phaseAdjust);
	}

	prevstate[0] = integrator;
	prevstate[1] = phaseEst;
	prevstate[2] = feedbackI;
	prevstate[3] = feedbackQ;
	prevstate[4] = ncoOut[ncoOut.size() - 1];

}
