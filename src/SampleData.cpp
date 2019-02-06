#include "SampleData.h"

using namespace std;

SampleData::SampleData(Report * report){
	sampleID = report->sampleID;
	sample_index = report->sample_index;
}

void SampleData::add_donor(SampleData &d_sample){
	donor_indeces.push_back(d_sample.sample_index);
	// donor_indeces.erase(donor_indeces.begin() + d_ind);
	n_donors++;
	d_sample.n_recipients++;
}

bool SampleData::remove_donor(int d_ind, vector<SampleData> &sample){
	SampleData * d_sample = &sample[donor_indeces[d_ind]];
	donor_indeces.erase(donor_indeces.begin() + d_ind);
	n_donors--;
	d_sample->n_recipients--;
	
	return n_donors >= 0 && d_sample->n_recipients >= 0;
}

void SampleData::print_result(ostream& ofs, vector<SampleData> const &sample) const {
	vector<double>::const_iterator result_it = regression_result.begin();
	ofs << get_sampleID() << '\t' << accumulate(regression_result.begin(), regression_result.end(), 0.0) << '\t';
	
	if(regression_result.size() > n_donors){
		ofs << "AF:" << *result_it;
		result_it++;
	}
	
	for(int d_ind = 0; d_ind < n_donors; d_ind++){
		result_it + d_ind != regression_result.begin() && ofs << "; ";
		ofs << "Donor_" << sample[donor_indeces[d_ind]].get_sampleID() << ':' << *(result_it + d_ind); 
	}
	
	ofs << endl;
}

string SampleData::get_sampleID(void) const {
	if(sampleID.empty())
		return to_string(sample_index);
	else
		return sampleID;
}
