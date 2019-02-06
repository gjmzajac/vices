#ifndef __REGRESSION_H__
#define __REGRESSION_H__

#include <vector>
#include <armadillo>
#include "SampleData.h"
#include "VicesArgs.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace arma;

double calc_alpha_AF_estimate(SampleData const &sample);
vector<double> calc_alpha_AF_donor_estimates(SampleData const &recipient, SampleData const &donor);
vector<double> calc_alpha_donor_estimate(SampleData const &recipient, SampleData const &donor);
bool no_missing_genos(vector<SampleData> const &sample, int r_ind, int marker_index);
void prep_donorSearch(vector<SampleData> &sample, 
	vector<double> const &AF, 
	int n_reports, 
	int n_markers, 
	int n_contam, 
	vector<int> &reordered_to_original_index);
vector<double> regression_AF_multiple_donors_model_selection(
	const int n_markers, 
	vector<SampleData> &sample, 
	int r_ind, 
	vector<double> const &AF, 
	VicesArgs const &a);
vector<double> regression_selector(const int n_markers, 
	vector<SampleData> &sample, 
	int r_ind, 
	vector<double> const &AF, 
	VicesArgs const &args);
void resize_arrays(int i, int n_contam, int n_reports, SampleData &sample);
void store_donor_info(SampleData &r_samp, SampleData &d_samp, double AF, int m_ind, int n_contam);
mat refit_no_AF(mat const &XtX, vec const &XtY, int n_donors);
void fill_XtX_XtY(vector<SampleData> const &sample, mat &XtX, vec &XtY, int r_ind, int m_ind, double AF);

#endif
