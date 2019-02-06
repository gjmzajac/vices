#include "regression.h"

using namespace std;
using namespace arma;

double calc_alpha_AF_estimate(SampleData const &sample){
	//load the data into variables with shorter and more readable names
	double nAA = sample.n_byGeno[0];
	double nAB = sample.n_byGeno[1];
	double nBB = sample.n_byGeno[2];
	double sum_AF_AA = sample.sum_AF_byGeno[0];
	double sum_AF_AB = sample.sum_AF_byGeno[1];
	double sum_AF_BB = sample.sum_AF_byGeno[2];
	double sum_Bint_AA = sample.sum_Bint_byGeno[0];
	double sum_Bint_AB = sample.sum_Bint_byGeno[1];
	double sum_Bint_BB = sample.sum_Bint_byGeno[2];
	double sum_AF_sq = sample.sum_AF_sq;
	double sum_AFxBint = sample.sum_AFxBint;

	//this is a special case of multiple regression and so the formula for 
	//the coefficient becomes more simple than the general one
	try{
		double determinant = 
			-nAB * (nBB * (sum_AF_AA * sum_AF_AA - nAA * sum_AF_sq) + nAA * sum_AF_BB * sum_AF_BB) - 
			nAA * nBB * sum_AF_AB * sum_AF_AB;
		if (determinant == 0)
			throw runtime_error("Division by zero");
		return (-nAB * nBB * sum_AF_AA * sum_Bint_AA - 
				nAA * nBB * sum_AF_AB * sum_Bint_AB - 
				nAA * nAB * sum_AF_BB * sum_Bint_BB + 
				nAA * nAB * nBB * sum_AFxBint) / determinant;
	}
	catch(exception& e){
		cerr << "error: matrix singular when regressing sample " << sample.get_sampleID() << " intensities on AF" << endl;
		return NAN;
	}
}

bool no_missing_genos(vector<SampleData> const &sample, int r_ind, int marker_index){
	for(int j=0; j < sample[r_ind].n_donors; j++){
		int d_ind = sample[r_ind].donor_indeces[j];
		if(sample[d_ind].Geno[marker_index] < 0)
			return false;
	}
	return true;
}

void prep_donorSearch(vector<SampleData> &sample, 
		vector<double> const &AF, 
		int n_reports, 
		int n_markers, 
		int n_contam, 
		vector<int> &reordered_to_original_index){
	//load the arrays
	for(int i = 0; i < n_reports; i++){
		int orig_i = reordered_to_original_index[i];
		resize_arrays(i, n_contam, n_reports, sample[orig_i]);
	}

	//store the info to regress on AF and every possible donor
	#pragma omp parallel for
	for(int i=0; i < n_reports; i++){
		int orig_i = reordered_to_original_index[i];
		for(int k=0; k < n_markers; k++){
			if(sample[orig_i].Geno[k] >= 0){
				for(int j=0; j < i && j < n_contam; j++){
					int orig_j = reordered_to_original_index[j];
					if(sample[orig_j].Bint[k] >= 0)
						store_donor_info(sample[orig_j], sample[orig_i], AF[k], k, n_contam);
				}
			}
		}
	}
}

vector<double> regression_AF_multiple_donors_model_selection(
		const int n_markers, 
		vector<SampleData> &sample, 
		int r_ind, 
		vector<double> const &AF, 
		VicesArgs const &a){
	SampleData * r_samp = &sample[r_ind];
	//Multiple regression of the form Y~X*B can be fit by
	//B = (X^t * X)^-1 * X^t * Y
	//here we denote as coef = (XtX)^-1 * XtY
	//allocate XtX and XtY matrices for regression
	mat XtX(4 + r_samp->n_donors, 4 + r_samp->n_donors, fill::zeros);
	vec XtY(4 + r_samp->n_donors, fill::zeros);

	//fill in XtY vector and the lower triangle of the matrix XtX
	for(int i=0; i < n_markers; i++)
		fill_XtX_XtY(sample, XtX, XtY, r_ind, i, AF[i]);

	//fill in the upper triangle
	XtX = symmatl(XtX);

	//fit the regression coefficients
	mat coef;
	try {
		coef = solve(XtX, XtY, solve_opts::no_approx);
	}
	catch(exception& e) {
		cerr << "error: check for duplicate samples" << endl;
		return vector<double> (r_samp->n_donors + 1, NAN);
	}

	int min_d_ind = conv_to<int>::from(index_min(coef.rows(4, coef.n_rows - 1)));
	double min_d_alpha = conv_to<double>::from(coef.row(4 + min_d_ind));
	if(min_d_alpha < a.contam_threshold){
		r_samp->remove_donor(min_d_ind, sample);

		return regression_selector(n_markers, sample, r_ind, AF, a);
	}
	else if(coef[3] <= 0)
		coef = refit_no_AF(XtX, XtY, r_samp->n_donors);
	
	//return the regression results for all donors and AF
	return conv_to< vector<double> >::from(coef.rows(3, coef.n_rows - 1));
}

vector<double> calc_alpha_AF_donor_estimates(SampleData const &recipient, SampleData const &donor){
	int r_ind = recipient.original_to_reordered_index;
	int d_ind = donor.original_to_reordered_index;
	SampleData const * u_samp;
	int l_ind;
	int l_AB;
	int l_BB;
	int u_AB;
	int u_BB;

	if(r_ind < d_ind){
		u_samp = &donor;
		l_ind = r_ind;
		l_AB = 1;
		l_BB = 2;
		u_AB = 3;
		u_BB = 6;
	}
	else if(r_ind > d_ind){
		u_samp = &recipient;
		l_ind = d_ind;
		l_AB = 3;
		l_BB = 6;
		u_AB = 1;
		u_BB = 2;
	}
	else{
		return {NAN, NAN};
	}
	double nAA = 
		u_samp->n_byRGeno_byDGeno[9*l_ind] + 
		u_samp->n_byRGeno_byDGeno[l_AB + 9*l_ind] + 
		u_samp->n_byRGeno_byDGeno[l_BB + 9*l_ind];
	double nAB = 
		u_samp->n_byRGeno_byDGeno[u_AB + 9*l_ind] + 
		u_samp->n_byRGeno_byDGeno[l_AB + u_AB + 9*l_ind] + 
		u_samp->n_byRGeno_byDGeno[l_BB + u_AB + 9*l_ind];
	double nBB = 
		u_samp->n_byRGeno_byDGeno[u_BB + 9*l_ind] + 
		u_samp->n_byRGeno_byDGeno[l_AB + u_BB + 9*l_ind] + 
		u_samp->n_byRGeno_byDGeno[l_BB + u_BB + 9*l_ind];
	double sum_AF_AA = 
		u_samp->sum_AF_byRGeno_byDGeno[9*l_ind] + 
		u_samp->sum_AF_byRGeno_byDGeno[l_AB + 9*l_ind] + 
		u_samp->sum_AF_byRGeno_byDGeno[l_BB + 9*l_ind];
	double sum_AF_AB = 
		u_samp->sum_AF_byRGeno_byDGeno[ u_AB + 9*l_ind] + 
		u_samp->sum_AF_byRGeno_byDGeno[l_AB + u_AB + 9*l_ind] + 
		u_samp->sum_AF_byRGeno_byDGeno[l_BB + u_AB + 9*l_ind];
	double sum_AF_BB = 
		u_samp->sum_AF_byRGeno_byDGeno[ u_BB + 9*l_ind] + 
		u_samp->sum_AF_byRGeno_byDGeno[l_AB + u_BB + 9*l_ind] + 
		u_samp->sum_AF_byRGeno_byDGeno[l_BB + u_BB + 9*l_ind];
	double sum_AF_sq = u_samp->sum_AF_sq_Recipient_Donor[l_ind];
	double sum_Gd_AA = 
		((double) u_samp->n_byRGeno_byDGeno[l_AB + 9*l_ind]) * 0.5 + 
		u_samp->n_byRGeno_byDGeno[l_BB + 9*l_ind];
	double sum_Gd_AB = 
		((double) u_samp->n_byRGeno_byDGeno[l_AB + u_AB + 9*l_ind]) * 0.5 + 
		u_samp->n_byRGeno_byDGeno[l_BB + u_AB + 9*l_ind];
	double sum_Gd_BB = 
		((double) u_samp->n_byRGeno_byDGeno[l_AB + u_BB + 9*l_ind]) * 0.5 + 
		u_samp->n_byRGeno_byDGeno[l_BB + u_BB + 9*l_ind];
	double sum_AFxGd = 
		(u_samp->sum_AF_byRGeno_byDGeno[l_AB + 9*l_ind] + 
			u_samp->sum_AF_byRGeno_byDGeno[l_AB + u_AB + 9*l_ind] + 
			u_samp->sum_AF_byRGeno_byDGeno[l_AB + u_BB + 9*l_ind]) * 0.5 + 
		u_samp->sum_AF_byRGeno_byDGeno[l_BB + 9*l_ind] + 
		u_samp->sum_AF_byRGeno_byDGeno[l_BB + u_AB + 9*l_ind] + 
		u_samp->sum_AF_byRGeno_byDGeno[l_BB + u_BB + 9*l_ind];
	double sum_Gd_sq = 
		((double) (u_samp->n_byRGeno_byDGeno[l_AB + 9*l_ind] + 
			u_samp->n_byRGeno_byDGeno[l_AB + u_AB + 9*l_ind] + 
			u_samp->n_byRGeno_byDGeno[l_AB + u_BB + 9*l_ind])) * 0.25 + 
		((double) (u_samp->n_byRGeno_byDGeno[l_BB + 9*l_ind] + 
			u_samp->n_byRGeno_byDGeno[l_BB + u_AB + 9*l_ind] + 
			u_samp->n_byRGeno_byDGeno[l_BB + u_BB + 9*l_ind]));
	double sum_Bint_AA = 
		recipient.sum_Bint_byRGeno_byDGeno[9*d_ind] + 
		recipient.sum_Bint_byRGeno_byDGeno[ 3 + 9*d_ind] + 
		recipient.sum_Bint_byRGeno_byDGeno[ 6 + 9*d_ind];
	double sum_Bint_AB = 
		recipient.sum_Bint_byRGeno_byDGeno[1 + 9*d_ind] + 
		recipient.sum_Bint_byRGeno_byDGeno[1 + 3 + 9*d_ind] + 
		recipient.sum_Bint_byRGeno_byDGeno[1 + 6 + 9*d_ind];
	double sum_Bint_BB = 
		recipient.sum_Bint_byRGeno_byDGeno[2 + 9*d_ind] + 
		recipient.sum_Bint_byRGeno_byDGeno[2 + 3 + 9*d_ind] + 
		recipient.sum_Bint_byRGeno_byDGeno[2 + 6 + 9*d_ind];
	double sum_AFxBint = recipient.sum_AFxBint_Recipient_Donor[d_ind];
	double sum_GdxBint = 
		(recipient.sum_Bint_byRGeno_byDGeno[ 3 + 9*d_ind] + 
			recipient.sum_Bint_byRGeno_byDGeno[1 + 3 + 9*d_ind] + 
			recipient.sum_Bint_byRGeno_byDGeno[2 + 3 + 9*d_ind]) * 0.5 + 
		recipient.sum_Bint_byRGeno_byDGeno[ 6 + 9*d_ind] + 
		recipient.sum_Bint_byRGeno_byDGeno[1 + 6 + 9*d_ind] + 
		recipient.sum_Bint_byRGeno_byDGeno[2 + 6 + 9*d_ind];

	//this is a special case of multiple regression and so the formula for the coefficient becomes more simple than the general one
	double determinant = 
		(- nBB * ( nAB * (- nAA * sum_AF_sq * sum_Gd_sq + nAA * sum_AFxGd * sum_AFxGd + 
				sum_AF_AA * sum_AF_AA * sum_Gd_sq - 
				2 * sum_AF_AA * sum_Gd_AA * sum_AFxGd + 
				sum_AF_sq * sum_Gd_AA * sum_Gd_AA ) + 
			sum_Gd_AB * sum_Gd_AB * ( nAA * sum_AF_sq - sum_AF_AA * sum_AF_AA ) + 
			2 * sum_AF_AB * sum_Gd_AB * ( sum_AF_AA * sum_Gd_AA - nAA * sum_AFxGd ) + 
			sum_AF_AB * sum_AF_AB * ( nAA * sum_Gd_sq - sum_Gd_AA * sum_Gd_AA )) + 
		nAB * ( sum_Gd_BB * sum_Gd_BB * ( sum_AF_AA * sum_AF_AA - nAA * sum_AF_sq ) + 
			sum_AF_BB * ( 2 * nAA * sum_Gd_BB * sum_AFxGd - 2 * sum_AF_AA * sum_Gd_AA * sum_Gd_BB ) + 
			sum_AF_BB * sum_AF_BB * ( sum_Gd_AA * sum_Gd_AA - nAA * sum_Gd_sq )) + 
		nAA * ( sum_AF_BB * sum_Gd_AB - sum_AF_AB * sum_Gd_BB ) * ( sum_AF_BB * sum_Gd_AB - sum_AF_AB * sum_Gd_BB ) );
	//cout << determinant << endl;
	if (determinant == 0)
		return { NAN, NAN};
	else{ 
		//regression parameter estimate for the AF
		return {(( nBB * sum_AF_AA * sum_Gd_AB * sum_Gd_AB - 
				nBB * sum_AF_AB * sum_Gd_AA * sum_Gd_AB + 
				nAB * sum_AF_AA * sum_Gd_BB * sum_Gd_BB - 
				nAB * sum_AF_BB * sum_Gd_AA * sum_Gd_BB + 
				nAB * nBB * sum_Gd_AA * sum_AFxGd - 
				nAB * nBB * sum_AF_AA * sum_Gd_sq ) * sum_Bint_AA +
			( nBB * sum_AF_AB * sum_Gd_AA * sum_Gd_AA - 
				nBB * sum_AF_AA * sum_Gd_AB * sum_Gd_AA + 
				nAA * sum_AF_AB * sum_Gd_BB * sum_Gd_BB - 
				nAA * sum_AF_BB * sum_Gd_AB * sum_Gd_BB + 
				nAA * nBB * sum_Gd_AB * sum_AFxGd - 
				nAA * nBB * sum_AF_AB * sum_Gd_sq ) * sum_Bint_AB +
			( nAB * sum_AF_BB * sum_Gd_AA * sum_Gd_AA - 
				nAB * sum_AF_AA * sum_Gd_BB * sum_Gd_AA + 
				nAA * sum_AF_BB * sum_Gd_AB * sum_Gd_AB - 
				nAA * sum_AF_AB * sum_Gd_AB * sum_Gd_BB + 
				nAA * nAB * sum_Gd_BB * sum_AFxGd - 
				nAA * nAB * sum_AF_BB * sum_Gd_sq ) * sum_Bint_BB +
			(- nAB * nBB * sum_Gd_AA * sum_Gd_AA - 
				nAA * nBB * sum_Gd_AB * sum_Gd_AB - 
				nAA * nAB * sum_Gd_BB * sum_Gd_BB + 
				nAA * nAB * nBB * sum_Gd_sq ) * sum_AFxBint +
			( nAB * nBB * sum_AF_AA * sum_Gd_AA + 
				nAA * nBB * sum_AF_AB * sum_Gd_AB + 
				nAA * nAB * sum_AF_BB * sum_Gd_BB - 
				nAA * nAB * nBB * sum_AFxGd ) * sum_GdxBint ) / determinant,
		//regression parameter estimate for the donor genotype
			(( nBB * sum_Gd_AA * sum_AF_AB * sum_AF_AB - 
				nBB * sum_AF_AA * sum_Gd_AB * sum_AF_AB + 
				nAB * sum_AF_BB * sum_AF_BB * sum_Gd_AA - 
				nAB * nBB * sum_AF_sq * sum_Gd_AA - 
				nAB * sum_AF_AA * sum_AF_BB * sum_Gd_BB + 
				nAB * nBB * sum_AF_AA * sum_AFxGd ) * sum_Bint_AA +
			( nBB * sum_Gd_AB * sum_AF_AA * sum_AF_AA - 
				nBB * sum_AF_AB * sum_Gd_AA * sum_AF_AA + 
				nAA * sum_AF_BB * sum_AF_BB * sum_Gd_AB - 
				nAA * nBB * sum_AF_sq * sum_Gd_AB - 
				nAA * sum_AF_AB * sum_AF_BB * sum_Gd_BB + 
				nAA * nBB * sum_AF_AB * sum_AFxGd ) * sum_Bint_AB + 
			( nAB * sum_Gd_BB * sum_AF_AA * sum_AF_AA - 
				nAB * sum_AF_BB * sum_Gd_AA * sum_AF_AA - 
				nAA * sum_AF_AB * sum_AF_BB * sum_Gd_AB + 
				nAA * sum_AF_AB * sum_AF_AB * sum_Gd_BB - 
				nAA * nAB * sum_AF_sq * sum_Gd_BB + 
				nAA * nAB * sum_AF_BB * sum_AFxGd ) * sum_Bint_BB + 
			( nAB * nBB * sum_AF_AA * sum_Gd_AA + 
				nAA * nBB * sum_AF_AB * sum_Gd_AB + 
				nAA * nAB * sum_AF_BB * sum_Gd_BB - 
				nAA * nAB * nBB * sum_AFxGd ) * sum_AFxBint +
			(- nAB * nBB * sum_AF_AA * sum_AF_AA - 
				nAA * nBB * sum_AF_AB * sum_AF_AB - 
				nAA * nAB * sum_AF_BB * sum_AF_BB + 
				nAA * nAB * nBB * sum_AF_sq ) * sum_GdxBint) / determinant};
	}
}

vector<double> calc_alpha_donor_estimate(SampleData const &recipient, SampleData const &donor){
	int r_ind = recipient.original_to_reordered_index;
	int d_ind = donor.original_to_reordered_index;
	SampleData const * u_samp;
	int l_ind;
	int l_AB;
	int l_BB;
	int u_AB;
	int u_BB;

	if(r_ind < d_ind){
		u_samp = &donor;
		l_ind = r_ind;
		l_AB = 1;
		l_BB = 2;
		u_AB = 3;
		u_BB = 6;
	}
	else if(r_ind > d_ind){
		u_samp = &recipient;
		l_ind = d_ind;
		l_AB = 3;
		l_BB = 6;
		u_AB = 1;
		u_BB = 2;
	}
	else{
		return {NAN};
	}
	double nAA = 
		u_samp->n_byRGeno_byDGeno[9*l_ind] + 
		u_samp->n_byRGeno_byDGeno[l_AB + 9*l_ind] + 
		u_samp->n_byRGeno_byDGeno[l_BB + 9*l_ind];
	double nAB = 
		u_samp->n_byRGeno_byDGeno[u_AB + 9*l_ind] + 
		u_samp->n_byRGeno_byDGeno[l_AB + u_AB + 9*l_ind] + 
		u_samp->n_byRGeno_byDGeno[l_BB + u_AB + 9*l_ind];
	double nBB = 
		u_samp->n_byRGeno_byDGeno[u_BB + 9*l_ind] + 
		u_samp->n_byRGeno_byDGeno[l_AB + u_BB + 9*l_ind] + 
		u_samp->n_byRGeno_byDGeno[l_BB + u_BB + 9*l_ind];
	double sum_Gd_AA = 
		((double) u_samp->n_byRGeno_byDGeno[l_AB + 9*l_ind]) * 0.5 + 
		u_samp->n_byRGeno_byDGeno[l_BB + 9*l_ind];
	double sum_Gd_AB = 
		((double) u_samp->n_byRGeno_byDGeno[l_AB + u_AB + 9*l_ind]) * 0.5 + 
		u_samp->n_byRGeno_byDGeno[l_BB + u_AB + 9*l_ind];
	double sum_Gd_BB = 
		((double) u_samp->n_byRGeno_byDGeno[l_AB + u_BB + 9*l_ind]) * 0.5 + 
		u_samp->n_byRGeno_byDGeno[l_BB + u_BB + 9*l_ind];
	double sum_Bint_AA = 
		recipient.sum_Bint_byRGeno_byDGeno[9*d_ind] + 
		recipient.sum_Bint_byRGeno_byDGeno[3 + 9*d_ind] + 
		recipient.sum_Bint_byRGeno_byDGeno[6 + 9*d_ind];
	double sum_Bint_AB = 
		recipient.sum_Bint_byRGeno_byDGeno[1 + 9*d_ind] + 
		recipient.sum_Bint_byRGeno_byDGeno[1 + 3 + 9*d_ind] + 
		recipient.sum_Bint_byRGeno_byDGeno[1 + 6 + 9*d_ind];
	double sum_Bint_BB = 
		recipient.sum_Bint_byRGeno_byDGeno[2 + 9*d_ind] + 
		recipient.sum_Bint_byRGeno_byDGeno[2 + 3 + 9*d_ind] + 
		recipient.sum_Bint_byRGeno_byDGeno[2 + 6 + 9*d_ind];
	double sum_Gd_sq = 
		((double) (u_samp->n_byRGeno_byDGeno[l_AB + 9*l_ind] + 
			u_samp->n_byRGeno_byDGeno[l_AB + u_AB + 9*l_ind] + 
			u_samp->n_byRGeno_byDGeno[l_AB + u_BB + 9*l_ind])) * 0.25 + 
		((double) (u_samp->n_byRGeno_byDGeno[l_BB + 9*l_ind] + 
			u_samp->n_byRGeno_byDGeno[l_BB + u_AB + 9*l_ind] + 
			u_samp->n_byRGeno_byDGeno[l_BB + u_BB + 9*l_ind]));
	double sum_GdxBint = 
		(recipient.sum_Bint_byRGeno_byDGeno[3 + 9*d_ind] + 
			recipient.sum_Bint_byRGeno_byDGeno[1 + 3 + 9*d_ind] + 
			recipient.sum_Bint_byRGeno_byDGeno[2 + 3 + 9*d_ind]) * 0.5 + 
		recipient.sum_Bint_byRGeno_byDGeno[6 + 9*d_ind] + 
		recipient.sum_Bint_byRGeno_byDGeno[1 + 6 + 9*d_ind] + 
		recipient.sum_Bint_byRGeno_byDGeno[2 + 6 + 9*d_ind];

	//this is a special case of multiple regression
	//so the formula for the coefficient becomes more simple than the general one
	double determinant = 
		-nAB *(nBB *(sum_Gd_AA * sum_Gd_AA - nAA * sum_Gd_sq) + nAA * sum_Gd_BB * sum_Gd_BB) - 
		nAA * nBB * sum_Gd_AB * sum_Gd_AB;
	//cout << determinant << endl;
	if (determinant == 0)
		return {NAN};
	else
		return {(-nAB * nBB * sum_Gd_AA * sum_Bint_AA - 
				nAA * nBB * sum_Gd_AB * sum_Bint_AB - 
				nAA * nAB * sum_Gd_BB * sum_Bint_BB + 
				nAA * nAB * nBB *       sum_GdxBint) / determinant};

}

vector<double> regression_selector(const int n_markers, 
		vector<SampleData> &sample, 
		int r_ind, 
		vector<double> const &AF, 
		VicesArgs const &args){
	SampleData * r_samp = &sample[r_ind];

	if(r_samp->n_donors == 0)
		return {r_samp->alpha_AF_estimate};
	else if(r_samp->n_donors == 1){
		int d_ind = r_samp->donor_indeces[0];
		if(r_samp->alpha_AF_donor_estimate[d_ind][0] <= 0) //calc 1 donor no AF
			return calc_alpha_donor_estimate(*r_samp, sample[d_ind]);
		else //stored AF + 1 donor
			return r_samp->alpha_AF_donor_estimate[d_ind];
	}
	else //multiple regression
		return regression_AF_multiple_donors_model_selection(n_markers, sample, r_ind, AF, args);
}

void resize_arrays(int i, int n_contam, int n_reports, SampleData &sample){
	if(i < n_contam){
		sample.n_byRGeno_byDGeno.resize(9*i, 0);
		sample.sum_AF_byRGeno_byDGeno.resize(9*i, 0);
		sample.sum_AF_sq_Recipient_Donor.resize(i, 0);
		sample.sum_Bint_byRGeno_byDGeno.resize(9*n_reports, 0);
		sample.sum_AFxBint_Recipient_Donor.resize(n_reports, 0);
	}
	else {
		sample.n_byRGeno_byDGeno.resize(9*n_contam, 0);
		sample.sum_AF_byRGeno_byDGeno.resize(9*n_contam, 0);
		sample.sum_AF_sq_Recipient_Donor.resize(n_contam, 0);
	}
}

void store_donor_info(SampleData &r_samp, SampleData &d_samp, double AF, int m_ind, int n_contam){
	int r_ind = r_samp.original_to_reordered_index;
	int d_ind = d_samp.original_to_reordered_index;
	int RGeno = r_samp.Geno[m_ind];
	int DGeno = d_samp.Geno[m_ind];
	double RBint = r_samp.Bint[m_ind];
	double DBint = d_samp.Bint[m_ind];

	d_samp.n_byRGeno_byDGeno[DGeno + 3*RGeno + 9*r_ind]++;
	d_samp.sum_AF_byRGeno_byDGeno[DGeno + 3*RGeno + 9*r_ind] += AF;
	d_samp.sum_AF_sq_Recipient_Donor[r_ind] += AF * AF;
	r_samp.sum_Bint_byRGeno_byDGeno[RGeno + 3*DGeno + 9*d_ind] += RBint;
	r_samp.sum_AFxBint_Recipient_Donor[d_ind] += AF * RBint;

	if(d_ind < n_contam && DBint >= 0){
		d_samp.sum_Bint_byRGeno_byDGeno[DGeno + 3*RGeno + 9*r_ind] += DBint;
		d_samp.sum_AFxBint_Recipient_Donor[r_ind] += AF * DBint;
	}
}

mat refit_no_AF(mat const &XtX, vec const &XtY, int n_donors){
	uvec exclude_AF(3 + n_donors);
	exclude_AF(0) = 0;
	exclude_AF(1) = 1;
	exclude_AF(2) = 2;
	for(int i=0; i < n_donors; i++)
		exclude_AF(3+i) = 4+i;

	return solve(XtX.submat(exclude_AF, exclude_AF), XtY.elem(exclude_AF));
}

void fill_XtX_XtY(vector<SampleData> const &sample, mat &XtX, vec &XtY, int r_ind, int m_ind, double AF){
	SampleData const * r_samp = &sample[r_ind];
	int RGeno = r_samp->Geno[m_ind];
	double RBint = r_samp->Bint[m_ind];

	if(RBint >= 0 && no_missing_genos(sample, r_ind, m_ind)){
		XtX.diag()(RGeno)++; //n_byGeno
		XtX(3, RGeno) += AF; //sum_AF_byGeno
		XtX(3, 3) += AF * AF; //sum_AF_sq

		XtY(RGeno) += RBint; //sum_Bint_byGeno
		XtY(3) += AF * RBint; //sum_AFxBint

		for(int d1_ind = 0; d1_ind < r_samp->n_donors; d1_ind++){
			double D1Geno = sample[r_samp->donor_indeces[d1_ind]].Geno[m_ind];
			XtX(4+d1_ind, RGeno) += 0.5 * D1Geno; //sum_DGeno_byRGeno
			XtX(4+d1_ind, 3) += AF * 0.5 * D1Geno; //sum_AFxDGeno
			XtY(4+d1_ind) += RBint * 0.5 * D1Geno; //sum_BintxDGeno
			for(int d2_ind = 0; d2_ind <= d1_ind; d2_ind++){
				double D2Geno = sample[r_samp->donor_indeces[d2_ind]].Geno[m_ind];
				XtX(4+d1_ind, 4+d2_ind) += 0.25 * D1Geno * D2Geno; //sum_DGenoxDGeno
			}
		}
	}
}
