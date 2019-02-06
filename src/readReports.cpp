#include "readReports.h"

using namespace std;

int openReports(VicesArgs const & args, vector<Report*> &report){
	string report_file_path;

	//open the list containing file names of reports
	IFILE report_files_list = NULL;
	IFILE sample_list = NULL;
	report_files_list = ifopen(args.report_list_fname, "r", InputFile::DEFAULT);
	if(report_files_list == NULL){
		cerr << "error: Could not open input file: " << args.report_list_fname << endl;
		exit(1);
	}
	
	//get the sample names if available
	if(args.sample_list_fname != NULL)
		sample_list = ifopen(args.sample_list_fname, "r", InputFile::DEFAULT);

	while(!report_files_list->ifeof() && report_files_list >> report_file_path){
		//open each file and get to the start of the data

		Report * ptrReport = NULL;
		try{
			ptrReport = new Report(report_file_path, args);
		}
		catch(exception& e){
			cerr << "error: " << e.what() << report_file_path << endl;
			exit(1);
		}
		sample_list != NULL && !sample_list->ifeof() && sample_list >> ptrReport->sampleID;

		report.push_back( ptrReport );
	}

	//close the file connection
	report_files_list->ifclose();
	delete report_files_list;
	report_files_list = NULL;

	return Report::num_samples;
}

void closeReports(vector<Report*> &report, int num_reports){
	for(int i=0; i<num_reports; i++){
		delete report[i];
		report[i] = NULL;
	}
}

int loadData_calcAFregTerms(vector<Report*> &report, 
		vector<SampleData> &sample, 
		vector<double> &AF, 
		int &n_markers_used,
		int num_reports, VicesArgs const &a){
	int i=0;
	int nGeno=0;
	int tot_markers=0;
	int sumGeno=0;
	
	while(report[i]->readLine() && (a.n_markers_limit < 0 || n_markers_used < a.n_markers_limit)){
		
		report[i]->readSnpname();
		if (report[0]->currSnpname != report[i]->currSnpname)
			throw runtime_error("Report file " + string(report[i]->fname) + " has different snp name(s)");
		
		//read Geno, Bint data from the file
		report[i]->currBint = -1;
		
		if(report[i]->readGeno() && report[i]->currGeno >= 0){
			nGeno += 2;
			sumGeno += report[i]->currGeno;
			report[i]->readBint();
		}
		i++;

		//store Geno, Bint, etc into the sample's object
		if(i >= num_reports){
			if(nGeno > 0){
				double currAF = ((double) sumGeno)/((double) nGeno);

				if(currAF >= a.AF_threshold && currAF <= (1 - a.AF_threshold)){
					AF.push_back(currAF);
					n_markers_used++;
					for(int j = 0; j < num_reports; j++)
						store_variant_info(report[j], sample[j], currAF);
				}
			}
			i=0;
			nGeno=0;
			sumGeno=0;
			tot_markers++;
		}
	}
	
	if(i != 0)
		throw runtime_error("Report file " + string(report[i]->fname) + " has different length");
	else if (report[i]->ifile->ifeof()) {
		for(int i=0; i<num_reports; i++){
			if(!report[i]->ifile->ifeof())
				throw runtime_error("Report file " + string(report[i]->fname) + " has different length");
		}
	}
	return tot_markers;
}

void store_variant_info(Report * report, SampleData &sample, double currAF){
	int Geno = report->currGeno;
	double Bint = report->currBint;
	
	sample.Geno.push_back(Geno);
	sample.Bint.push_back(Bint);
	
	if(Bint >= 0){
		sample.n_byGeno[Geno]++;
		sample.sum_AF_byGeno[Geno] += currAF;
		sample.sum_Bint_byGeno[Geno] += Bint;
		sample.sum_AF_sq += currAF * currAF;
		sample.sum_AFxBint += currAF * Bint;
	}
}
