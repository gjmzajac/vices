#ifndef __REPORT_H__
#define __REPORT_H__

#include <string>
#include "InputFile.h"
#include <algorithm>
#include "VicesArgs.h"

using namespace std;

class Report {
	public:
	
		string fname;
		string sampleID;
		IFILE ifile;
		string currLine;
		int sample_index = 0;
		static int num_samples;
		int currLineNo = 0;
		int colsnp;
		int colab1;
		int colab2;
		int colBint;
		int currGeno = -4;
		double currBint = -1;
		string currSnpname;

		Report(string report_file_path, VicesArgs const &a);
		
		bool readLine(void);
		bool readSnpname(void);
		bool readGeno(void);
		bool readBint(void);
		int pos_nth_occurance_of_char(string const &str, int n, char c);
		double nth_field_to_double(string const &str, int n, char delim);
		string nth_field_to_string(string const &str, int n, char delim);
		int find_col(string const &header, string colname);
		int reportLine_Geno(string const &line, int colab1, int colab2);

		~Report();
};

#endif
