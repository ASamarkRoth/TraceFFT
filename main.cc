// std includes
#include <iostream>
#include <algorithm>
#include <chrono>

//ROOT
#include "TFile.h"
#include "TApplication.h"
#include "TSystem.h"
#include "TROOT.h"

//extras
#include "yaml-cpp/yaml.h"

#include "FFT.h"

using namespace std;

int main(int argc, char *argv[])
{
	std::vector<std::string> s_input;

	TFile* root_results;

	root_results = new TFile("results.root", "RECREATE");

	//should be able to read in many files simultaneously
	if(argc > 1) {
		for(int i = 1; i < argc; ++i) {
			s_input.push_back(argv[1]);
		}
	}
	//else {
	//	std::cout << "Please provide a correct input file name for the run" << std::endl;
	//	abort();
	//}


	TApplication theApp("tapp", &argc, argv);
	if(s_input.size() > 0)
	{
		// Here loop over the number of input files
		for(auto& s_i : s_input)
		{
			TraceFFT(s_i, root_results);

		}
	}
	else
	{
		cout << "Running default simulation which has a sample periodicity of 32, this usually takes a while ..." << endl;
		FFT(root_results);
	}

	cout << "DONE! \n 	Results can be found in the file results.root" << endl;

	gROOT->ProcessLine(".L complete_analysis.root");
	theApp.Run();


	return 0;
}
