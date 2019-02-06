#ifndef __READREPORTS_H__
#define __READREPORTS_H__

#include <string>
#include "InputFile.h"
#include <vector>
#include <algorithm>
#include "Report.h"
#include "VicesArgs.h"
#include "SampleData.h"

using namespace std;

int openReports(VicesArgs const &args, vector<Report*> &report);
void closeReports(vector<Report*> &report, int num_reports);
int loadData_calcAFregTerms(vector<Report*> &report, 
	vector<SampleData> &sample, 
	vector<double> &AF, 
	int &n_markers_used,
	int num_reports, VicesArgs const &a);
void store_variant_info(Report *report, SampleData &sample, double currAF);

#endif
