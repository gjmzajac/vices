#include "VicesArgs.h"

using namespace std;

VicesArgs::VicesArgs(int argc, char **argv){
	
	#ifdef _OPENMP
		static const char *optString = "r:o:f:c:s:an:1:2:b:m:t:h?";
	#else
		static const char *optString = "r:o:f:c:s:an:1:2:b:m:h?";
	#endif

	static const struct option longOpts[] = {
		{ "report-list", required_argument, NULL, 'r' },
		{ "out", required_argument, NULL, 'o' },
		{ "maf-threshold", required_argument, NULL, 'f' },
		{ "contam-threshold", required_argument, NULL, 'c' },
		{ "sample-list", required_argument, NULL, 's' },
		{ "af-only", no_argument, NULL, 'a' },
		{ "snp-name-col", required_argument, NULL, 'n' },
		{ "allele1-col", required_argument, NULL, '1' },
		{ "allele2-col", required_argument, NULL, '2' },
		{ "b-allele-intensity-col", required_argument, NULL, 'b' },
		{ "num-markers", required_argument, NULL, 'm' },
		#ifdef _OPENMP
			{ "threads", required_argument, NULL, 't' },
		#endif
		{ "help", no_argument, NULL, 'h' },
		{ NULL, no_argument, NULL, 0 }
	};

	int c;

	while ((c = getopt_long(argc, argv, optString, longOpts, NULL)) >= 0){
		switch (c) {
			case 'r': report_list_fname = optarg; break;
			case 'o': output_fname = optarg; break;
			case 'f': AF_threshold = stod(optarg); break;
			case 'c': contam_threshold = stod(optarg); break;
			case 's': sample_list_fname = optarg; break;
			case 'a': AF_only = true; break;
			case 'n': colsnp_name = optarg; break;
			case '1': colab1_name = optarg; break;
			case '2': colab2_name = optarg; break;
			case 'b': colBint_name = optarg; break;
			case 'm': n_markers_limit = stoi(optarg); break;
			#ifdef _OPENMP
				case 't': threads = stoi(optarg); break;
			#endif
			case 'h':
			case '?': usage(); break;
			default: break;
		}
	}
	if(output_fname == NULL || report_list_fname == NULL)
		usage();
}

void VicesArgs::usage(void){
	cout << endl;
	cout << "This program performs joint estimation of DNA contamination and its sources in genotyping arrays" << endl;
	cout << "Usage: vices -r <reports_list.txt> -o <contam_estimates.txt> [-f <0.1>] [-c <0.005>] [-s <sample_ids.txt>] [-a] [-n <\"SNP Name\">] [-1 <\"Allele1 - AB\">] [-2 <\"Allele2 - AB\">] [-b <\"B Allele Freq\">] [-m <-1>]";
	#ifdef _OPENMP
		cout << " [-t <1>]";
	#endif
	cout << endl;
	cout << "Options:" << endl;
	cout << "   -r, --report-list <file>              File with paths to Illumina report files" << endl;
	cout << "   -o, --output <file>                   Write output to a file [standard output]" << endl;
	cout << "   -f, --maf-threshold <float>           Min minor allele frequency for markers" << endl;
	cout << "   -c, --contam-threshold <float>        Threshold for estimating donor samples" << endl;
	cout << "   -s, --sample-list <file>              File with sample ids for report files" << endl;
	cout << "   -a, --af-only                         Specify analysis with AF only. No donor estimation" << endl;
	cout << "   -n, --snp-name-col <string>           Name of report file column containing SNP names" << endl;
	cout << "   -1, --allele1-col <string>            Name of report file column containing allele 1" << endl;
	cout << "   -2, --allele2-col <string>            Name of report file column containing allele 2" << endl;
	cout << "   -b, --b-allele-intensity-col <string> Name of report file column containing B allele intensity" << endl;
	cout << "   -m, --num-markers <int>               Maximum number of markers for contamination estimation" << endl;
	#ifdef _OPENMP
		cout << "   -t, --threads <int>                   Number of threads for parallel computation" << endl;
	#endif
	cout << "   -h, --help                            This help page" << endl;
	cout << endl;

	exit(1);
}
