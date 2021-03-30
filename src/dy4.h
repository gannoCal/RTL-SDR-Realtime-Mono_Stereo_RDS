/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_DY4_H
#define DY4_DY4_H

// some general and reusable stuff
// our beloved PI constant
#define PI 3.14159265358979323846
#define QUEUE_bLoCks 200
#define Taps 151
#include <queue>
#include <deque>
#include <iostream>
#include <thread>
#include <mutex>
#include <condition_variable>

// although we use DFT (no FFT ... yet), the number of points for a
// Fourier transform is defined as NFFT (same as matplotlib)
#define NFFT 512

#endif // DY4_DY4_H
