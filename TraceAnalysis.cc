#include "TraceAnalysis.h"

//#include "TFeb3Proc.h"

#include <cstdlib>

//TraceAnalysis::TraceAnalysis(std::string s_config) :  fep_amplitude(0), avg_counter(0), lin_step(0), w_avg(0) {
//
//	std::ifstream file(s_config.c_str());
//	std::string line, word;
//
//	while(getline(file, line)) {
//		std::istringstream iss(line);
//
//		if(iss.str().empty() || iss.str().find_first_not_of(' ') == std::string::npos) {
//			iss.clear();
//			iss.str("");
//			continue;
//		}
//		if(iss.str()[0]=='#') {
//			iss.clear();
//			iss.str("");
//			continue;
//		}
//		iss >> word;
//		if(word == "fep_gate") iss >> fep_low >> fep_high;
//		else if(word == "real_trace_length") iss >> real_trace_length;
//		else if(word == "comp_trace_gate") iss >> ct_low >> ct_high;
//		else if(word == "t_align") iss >> t_align;
//		else if(word == "t_baseline") iss >> t_bl_high;
//		else if(word == "fep_amplitude") iss >> fep_amplitude;
//		else if(word == "lin_step") iss >> lin_step;
//		else if(word == "w_avg") iss >> w_avg;
//	}
//
//	try {
//		std::cout << *this;
//	}
//	catch(...) {
//		std::cout << "Parameters for Trace Analysis not read correctly!" << std::endl;
//		abort();
//	}
//
//}

void TraceAnalysis::BaselineRestore(TH1* h_to_restore, TH1* h_bl_restored) {
	//double baseline = static_cast<double>( h_to_restore->Integral(0, t_bl_high) / t_bl_high ); //Why not working????

	double baseline = 0;
	for(unsigned int i = 0; i < t_bl_high; ++i) {
		baseline += h_to_restore->GetBinContent(i)/static_cast<double>( t_bl_high );
	}
	//std::cout << "Baseline = " << baseline << std::endl;

	for(int i = 0; i < h_to_restore->GetNbinsX(); ++i) {
		h_bl_restored->SetBinContent(i, h_to_restore->GetBinContent(i) - baseline);
	}
}

void TraceAnalysis::UpdateAvgTrace(TH1* old_avg_trace, TH1* new_trace) {
	old_avg_trace->Add(new_trace);
	++avg_counter;
}

//std::map<std::string, UShort_t> TraceAnalysis::CalculateTXX(TH1* h_bl_restored) {
std::vector<UShort_t> TraceAnalysis::CalculateTXX(TH1* h_bl_restored) {

	//std::map<std::string, UShort_t> txx;
	std::vector<UShort_t> txx;

	//it is set manually from file (for debug only)
	if(fep_amplitude > 0) {
		int hxx = 10;
		double limit = hxx*0.01*fep_amplitude;
		UShort_t i = 0;
		for( ;i < h_bl_restored->GetNbinsX(); ++i) {
			if(h_bl_restored->GetBinContent(i) > limit) {
				//txx[std::to_string(hxx)] = i;
				txx.push_back(i);
				hxx += 10;
				if(hxx == 100) break;
				limit = hxx*0.01*fep_amplitude;
			}
		}
		if(i == h_bl_restored->GetNbinsX()) {
			std::cout << "	#ERROR: Could not determine t10, t30 and t90!" << std::endl;
			//exit(1);
		}

	}

	return txx;
}

// Takes the histogram trace which should be linearly interpolated and writes it to another one which needs to be created before. Step is the number of bins per original bin.
void TraceAnalysis::LinearInterpolation(TH1* h_to_interpolate, TH1* h_interpolated, unsigned int step) {
	if(step > 1) {
		for(int i = 0; i < h_to_interpolate->GetNbinsX()-1; ++i) {
			//setting first bin within 2 original bins
			h_interpolated->SetBinContent( i*step, h_to_interpolate->GetBinContent(i) );
			for(unsigned int j = i*step+1; j < step*(i+1); ++j) {
				h_interpolated->SetBinContent(j, (h_to_interpolate->GetBinContent(i+1) - h_to_interpolate->GetBinContent(i)) * ( (j%step) * 1./step ) + h_to_interpolate->GetBinContent(i));
				//std::cout << "j,  lin_step " << j << ", " << (j%lin_step) * 1./lin_step << std::endl;
			}
		}
		//setting last bin within 2 original bins
		h_interpolated->SetBinContent(h_interpolated->GetNbinsX(), h_to_interpolate->GetBinContent(h_to_interpolate->GetNbinsX()) );
	}
}

//Performs a forward moving average with window length w
void TraceAnalysis::ForwardMovingAvg(TH1* h_to_move_avg, TH1* h_mov_avged, unsigned int w) {

	if(w > 1) {
		double avg = 0;
		double norm = 1./w;
		for(unsigned int i = 0; i < (unsigned int) h_to_move_avg->GetNbinsX(); ++i) {
			//adding
			avg += h_to_move_avg->GetBinContent(i)*norm;
			//subtract
			if(i >= w) avg -= h_to_move_avg->GetBinContent(i-w)*norm;
			//inserting
			if(i >= w-1) h_mov_avged->SetBinContent(i-w+1, avg);
		}
	}

}

void TraceAnalysis::AlignRise(TH1* h_to_align, TH1* h_aligned, int t_start, int t_align) {

	int shift = t_start-t_align;
	Int_t i = 0;
	for( ;i < h_to_align->GetNbinsX(); ++i) {
		if(i > shift) {
			h_aligned->SetBinContent(i-shift, h_to_align->GetBinContent(i));
		}
		else {
			h_aligned->SetBinContent(i, 0);
		}
	}

}

void TraceAnalysis::DeconvolveHist(TH1* blr_trace, TH1* to_deconv, double tau)
{
	for(int j = 0; j < blr_trace->GetNbinsX(); j++)
      {
        if(j == 0)
          {
            to_deconv->SetBinContent(j, blr_trace->GetBinContent(j));
          }
        else
          {
            to_deconv->SetBinContent(j, to_deconv->GetBinContent(j-1)+blr_trace->GetBinContent(j)-blr_trace->GetBinContent(j-1)*(1-(1./tau)));
          }
			}
}

void TraceAnalysis::CopyContent(TH1* orig, TH1* copy)
{
	Int_t nbins1 = copy->GetNbinsX();
	Int_t nbins2 = orig->GetNbinsX();
	Int_t nbins = nbins1;
	if(nbins2 < nbins) nbins = nbins2;
	for(int j = 0; j < nbins; j++)
	{
		copy->SetBinContent(j, orig->GetBinContent(j));
	}

}

//std::ostream& operator<<(std::ostream& os, TraceAnalysis& sta) {
//	os << "	-----------------------------------------" << std::endl;
//	os << "	#TraceAnalysis:" << std::endl;
//	os << "		Full energy peak gate: [" << sta.fep_low << "," << sta.fep_high << "]" << std::endl;
//	os << "		Real trace length: " << sta.real_trace_length << std::endl;
//	os << "		Compressed trace gate: [" << sta.ct_low << "," << sta.ct_high << "]" << std::endl;
//	os << "		Alignment time: " << sta.t_align << std::endl;
//	os << "		Baseline calculated on the basis of: [0," << sta.t_bl_high << "]" << std::endl;
//	os << "		Full energy peak amplitude: " << sta.fep_amplitude << std::endl;
//	os << "		Linear interpolation: " << sta.lin_step << std::endl;
//	os << "	-----------------------------------------" << std::endl;
//	return os;
//}
