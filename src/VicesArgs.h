#ifndef __VICESARGS_H__
#define __VICESARGS_H__

#include <getopt.h>
#include <iostream>

using namespace std;

class VicesArgs {
	public:
	
	char * report_list_fname = NULL;
	char * output_fname = NULL;
	double AF_threshold = 0.1;
	double contam_threshold = 0.005;
	char * sample_list_fname = NULL;
	bool AF_only = false;
	char const * colsnp_name = "SNP Name";
	char const * colab1_name = "Allele1 - AB";
	char const * colab2_name = "Allele2 - AB";
	char const * colBint_name = "B Allele Freq";
	int n_markers_limit = -1;
	#ifdef _OPENMP
		int threads = 1;
	#endif

	VicesArgs(int argc, char **argv);

	void usage(void);

};

#endif
