#ifndef FFT_H
#define FFT_H 1

#include "TH1D.h"
#include "TVirtualFFT.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TraceAnalysis.h"
#include "TFile.h"
#include "TKey.h"
#include <math.h>

void FFT(TFile* root_results);

void TraceFFT(std::string s_input, TFile* root_results);
 #endif
