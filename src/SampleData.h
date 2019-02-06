#ifndef __SAMPLEDATA_H__
#define __SAMPLEDATA_H__

#include <string>
#include <vector>
#include <numeric>
#include "Report.h"

using namespace std;

class SampleData {
	public:
	
		string sampleID;
		int sample_index;

		vector<int> Geno;
		vector<double> Bint;
		
		vector<int> n_byGeno = {0,0,0};
		vector<double> sum_AF_byGeno = {0,0,0};
		vector<double> sum_Bint_byGeno = {0,0,0};
		double sum_AF_sq = 0;
		double sum_AFxBint = 0;
		double alpha_AF_estimate = -1;

		vector<int> donor_indeces;
		// int n_incl_donors = 2;
		int n_donors = 0;
		int n_recipients = 0;
		int original_to_reordered_index = -1;
		
		vector<int> n_byRGeno_byDGeno;
		vector<double> sum_AF_byRGeno_byDGeno;
		vector<double> sum_Bint_byRGeno_byDGeno;
		vector<double> sum_AF_sq_Recipient_Donor;
		vector<double> sum_AFxBint_Recipient_Donor;
		
		vector< vector<double> > alpha_AF_donor_estimate;
		
		vector<double> regression_result;
		
		SampleData(Report * report);
		
		void add_donor(SampleData &d_sample);
		bool remove_donor(int d_ind, vector<SampleData> &sample);
		void print_result(ostream& ofs, vector<SampleData> const &sample) const;
		string get_sampleID(void) const;

};

#endif
