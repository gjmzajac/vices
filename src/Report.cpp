#include "Report.h"

using namespace std;

int Report::num_samples = 0;

Report::Report(string report_file_path, VicesArgs const &a){

	fname = report_file_path;
	if ((ifile = ifopen(report_file_path.c_str(), "r", InputFile::DEFAULT)) == NULL)
		throw runtime_error("could not read file");
	
	while(readLine() && currLine.find("[Data]"));
	
	if(readLine()){ 
		colsnp = find_col(currLine, a.colsnp_name);
		colab1 = find_col(currLine, a.colab1_name);
		colab2 = find_col(currLine, a.colab2_name);
		colBint = find_col(currLine, a.colBint_name);
	}
	else
		throw runtime_error("could not read header line in file ");
	sample_index = num_samples++;
}

bool Report::readLine(void){
	if(!ifile->ifeof() && ifile >> currLine)
		return true;
	else
		return false;
}

bool Report::readSnpname(void){
	currSnpname = nth_field_to_string(currLine, colsnp, '\t');
	return true;
}

bool Report::readGeno(void){
	currGeno = reportLine_Geno(currLine, colab1, colab2);
	return true;
}

bool Report::readBint(void){
	currBint = nth_field_to_double(currLine, colBint, '\t');
	return true;
}

int Report::pos_nth_occurance_of_char(string const &str, int n, char c){
		if(n == 0)
			return 0;

		int count_tab = 0;
		for(int pos = 0; pos < str.length(); pos++)
			if(str[pos] == c && ++count_tab == n)
				return pos+1;

		return -1;
}

double Report::nth_field_to_double(string const &str, int n, char delim){
	int pos1=pos_nth_occurance_of_char(str, n, delim);
	if(pos1 < 0)
		return -1;
	
	for (int pos2 = pos1 + 1; pos2 < str.length(); pos2++){
		if(isspace(str[pos2]))
			return(stod(str.substr(pos1, pos2 - pos1)));
		else if(!(isdigit(str[pos2]) || str[pos2] == '.'))
			return -1;
	}
	return -1;
}

string Report::nth_field_to_string(string const &str, int n, char delim){
	int pos1=pos_nth_occurance_of_char(str, n, delim);
	if(pos1 < 0)
		return "";
	
	for (int pos2 = pos1 + 1; pos2 < str.length(); pos2++){
		if(str[pos2] == delim || pos2 == str.length()-1)
			return(str.substr(pos1, pos2 - pos1));
	}
	return "";
}

int Report::reportLine_Geno(string const &line, int colab1, int colab2){
	int Geno = -4;
	const int a1 = line[pos_nth_occurance_of_char(line, colab1, '\t')];
	if (a1 == 'A')
		Geno += 2;
	else if (a1 == 'B')
		Geno += 3;

	const int a2 = line[pos_nth_occurance_of_char(line, colab2, '\t')];
	if (a2 == 'A')
		Geno += 2;
	else if (a2 == 'B')
		Geno += 3;

	return Geno;
}

int Report::find_col(string const &header, string col_name){
	string field;

	for(int col=0; !(field = nth_field_to_string(header, col, '\t')).empty(); col++){
		if(field == col_name)
			return col;
	}
	
	throw runtime_error("could not find column named \"" + col_name + "\" in file ");
	return -1;
}

Report::~Report(){
	ifile->ifclose();
	delete ifile;
	ifile = NULL;
}
