#ifndef TRACE_ANALYSIS_H
#define TRACE_ANALYSIS_H 1

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

#include "TH1.h"
//#include "TGo4EventProcessor.h"

class TraceAnalysis {

	public:
		TraceAnalysis() {;}
		TraceAnalysis(std::string s_config);


		friend std::ostream& operator<<(std::ostream& os, TraceAnalysis& sta);

		double fep_low, fep_high, t_align, fep_amplitude;
		unsigned int real_trace_length, ct_low, ct_high, t_bl_high, avg_counter, lin_step, w_avg;

		//For summing all pulses need to be aligned for proper averaging and this is done the most accurate through some kind of interpolation (see Validation of pulse shapes for inspiration)
		void AlignRise(TH1* h_to_align, TH1* h_aligned, int t_start, int t_align);

		//Baseline restore
		void BaselineRestore(TH1* h_to_restore, TH1* h_bl_restore);

		//Calculate baseline standard deviation ?

		//Get height of pulses (important to set times) (store in histogram). Could vary between coordinates, i.e. how to write trees? Maybe only store coordinate specific properties to the trees!? This avoids unnecessary processing.

		// MWD, is necessary to compute the amplitude of the pulse accurately. How to deal with pile-ups=flag?

		//Obtain T10, T30, T90
		std::vector<UShort_t> CalculateTXX(TH1* h_bl_restored);

		//Sum traces and make new average
		void UpdateAvgTrace(TH1* old_avg_trace, TH1* new_trace);

		//Linear interpolation for more accurate time alignment and calculations
		void LinearInterpolation(TH1* h_to_interpolate, TH1* h_interpolated, unsigned int step);

		//Forward moving average with window length w
		void ForwardMovingAvg(TH1* h_to_move_avg, TH1* h_mov_avged, unsigned int w);

		//Trace confidence bands for each coordinate (there has got to be an uncertainty!)

		//superimposition of all aligned traces?
		//maybe this should be done in the subsequent analysis steps

		void DeconvolveHist(TH1* blr_trace, TH1* to_deconv, double tau);
		void CopyContent(TH1* orig, TH1* copy);

		ClassDefNV(TraceAnalysis, 1);

};

#endif //TRACE_ANALYSIS_H
