/*VICES v1.0 by Gregory Zajac (c) 2018*/

#include <iostream>
#include <fstream>
#include "InputFile.h"
#include <string>
#include <vector>
#include <algorithm>
#include "regression.h"
#include "readReports.h"
#include "Report.h"
#include "VicesArgs.h"
#include "SampleData.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

int main(int argc, char **argv){
	cout << "VICES v1.0: Verify Intensity Contamination from Estimated Sources" << endl;
	cout << "(c) 2019 - Gregory Zajac, Goncalo Abecasis" << endl;
	
	VicesArgs args (argc, argv);
	
	#ifdef _OPENMP
        omp_set_num_threads(args.threads);
	#endif

	vector<double> AF;
	int n_contam = 0;
	int n_uncontam = 0;
	int tot_markers = 0;
	int n_markers_used = 0;
	ofstream ofs (args.output_fname);
	if (!ofs){
		cerr << "error: Could not open output file: " << args.output_fname << endl;
		exit(1);
	}

	vector<Report*> report;

	//open reports, read in data, and calculate terms for AF regression
	const int n_reports = openReports(args, report);
	cout << "Reading " << n_reports << " report files..." << endl;
	vector<SampleData> sample (report.begin(), report.end());
	try{
		tot_markers = loadData_calcAFregTerms(report, sample, AF, n_markers_used, n_reports, args);
	}
	catch(exception& e){
		cerr << "error: " << e.what() << endl;
		exit(1);
	}
	closeReports(report, n_reports);
	cout << "Done" << endl;

	cout << "Read " << tot_markers << " Markers in the report files" << endl;
	
	cout << "Using " << n_markers_used << " markers ";
	cout << "with MAF >= " << args.AF_threshold << endl;
	
	vector<int> reordered_to_original_index(n_reports, -1);
	
	//calculate the AF estimates
	for(int i=0; i<n_reports; i++){
		sample[i].alpha_AF_estimate = calc_alpha_AF_estimate(sample[i]);

		if(sample[i].alpha_AF_estimate >= args.contam_threshold){
			sample[i].original_to_reordered_index = n_contam;
			reordered_to_original_index[n_contam++] = i;
		}
		else{
			reordered_to_original_index[n_reports - (++n_uncontam)] = i;
			sample[i].original_to_reordered_index = n_reports - n_uncontam;
		}
	}
	
	cout << "Initial estimates based on AFs complete" << endl;
	cout << "Contamination above " << args.contam_threshold;
	cout << " detected in " << n_contam << " samples" << endl;
	if(!args.AF_only){
		cout << "Starting donor search..." << endl;

		//search for candidate donors
		//regress intensities of every contaminated sample on
		//the genotype calls of all other samples
		prep_donorSearch(sample, AF, n_reports, n_markers_used, n_contam, reordered_to_original_index);
		
		for(int i=0; i<n_contam; i++){
			int orig_i = reordered_to_original_index[i];
			sample[orig_i].alpha_AF_donor_estimate.resize(n_reports, {NAN,NAN});
			for(int j=0; j<n_reports; j++){
				sample[orig_i].alpha_AF_donor_estimate[j] = calc_alpha_AF_donor_estimates(sample[orig_i], sample[j]);

				if(sample[orig_i].alpha_AF_donor_estimate[j][1] >= args.contam_threshold){
					sample[orig_i].add_donor(sample[j]);
				}
			}
		}
		cout << "Done" << endl;

		//calculate final estimates with all estimated donors and print
		//intensities of contaminated samples are regressed on genotypes of
		//all putative donors to filter out false donors
		cout << "Pruning estimated donors and calculating final estimates..." << endl;
	}
	
	if(args.sample_list_fname == NULL)
		ofs << "Recipient_Index\tEstimated_contamination\tSources" << endl;
	else
		ofs << "Recipient_ID\tEstimated_contamination\tSources" << endl;

	#pragma omp parallel for
	for(int i=0; i < n_reports; i++)
		sample[i].regression_result = regression_selector(n_markers_used, sample, i, AF, args);

	for(int i=0; i < n_reports; i++)
		sample[i].print_result(ofs, sample);

	cout << "Done" << endl;
	cout << "Results written to " << realpath(args.output_fname, NULL) << endl;
	
	return 0;
}
