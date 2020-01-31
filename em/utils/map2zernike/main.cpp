#include <vector>
#include <sstream>
#include <iostream>
#include <istream>
#include <fstream>
#include <stdexcept>
#include <cassert>

#include "gzstream/gzstream.h"
#include "Grid.h"
#include "ZernikeDescriptor.h"
#include "util.h"

using namespace std;

bool VERBOSE = false;

static string formatFilenamePrefix(string s, string c){
	if(s.find("%s") == std::string::npos){
		return s + c;
	}else{
		char buf[1024];
		sprintf(buf, s.c_str(), c.c_str());
		return string(buf);
	}
}

static bool endsWith(const string &a, const string &b){
	return a.size() >= b.size() && a.compare(a.size()-b.size(), b.size(), b) == 0;
}

template<class T>
static void saveInvariants(ZernikeDescriptor<T> &zd, const string &path)
{
	std::ofstream outfile;
	outfile.exceptions(ios::badbit|ios::failbit);
	try{
		outfile.open(path.c_str());
	}catch(std::exception &e){
		throw std::runtime_error("Output file " + path + " could not be opened for writing.");
	}

	try{
		outfile << zd.size() << std::endl;
		for(auto &inv : zd){
			outfile << (inv/10) << std::endl;
		}
	}catch(std::exception &e){
		throw std::runtime_error("Write to output file " + path + " failed.");
	}
}

int main(int argc, char** argv)
{
	assert(sizeof(char) == 1);
	assert(sizeof(int) == 4);
	assert(sizeof(float) == 4);
	assert(sizeof(double) == 8);
	assert(sizeof(long) == 8);

	Binomial<double>::computePascalsTriangle(60);

	string programName(argv[0]);
	int maxOrder = 20; // default to max stable order
	double contour = 1.0; //default to assuming preprocessed map
	string contourString = "1.0";
	string inFileName;
	string outFilePrefix;
	bool isGzipped = false;
	vector<double> contours;
	vector<string> contourStrings;

	for(int i=1; i<argc; i++){
		string cur(argv[i]);

		if(cur == "-n" || cur == "--order"){
			if((i+1) < argc){
				i++;
				if(!(istringstream(argv[i]) >> maxOrder)){
					cerr << "[error] Expected integer in [1..20] after '" << cur << "', not '" << argv[i] << "'." << endl;
					exit(1);
				}
				if(maxOrder<1 || maxOrder>20){
					cerr << "[error] Expected integer in [1..20] after '" << cur << "', not '" << argv[i] << "'." << endl;
					exit(1);
				}
			}else{
					cerr << "[error] Expected integer in [1..20] after '" << cur << "'." << endl;
					exit(1);
			}
		}else if(cur == "-c" || cur == "--contour"){
			if((i+1) < argc){
				i++;
				if(!(istringstream(argv[i]) >> contour)){
					cerr << "[error] Expected real contour value after '" << cur << "'." << endl;
					exit(1);
				}
				contourString = argv[i];
			}else{
				cerr << "[error] Expected real contour value after '" << cur << "'." << endl;
				exit(1);
			}
		}else if(cur == "--contours"){
			if((i+1) < argc){
				i++;
				stringstream ss,tt;
				ss << argv[i];
				tt << argv[i];
				double tmp;
				string stmp;
				while(ss >> tmp){
					contours.push_back(tmp);
					tt >> stmp;
					contourStrings.push_back(stmp);
				}
				if(ss.get() != EOF){
					cerr << "[error] Expected space-delimited list of contours after '" << cur << "'." << endl;
					exit(1);
				}
			}else{
				cerr << "[error] Expected space-delimited list of contours after '" << cur << "'." << endl;
				exit(1);
			}
		}else if(cur == "-p" || cur == "--prefix"){
			if((i+1) < argc){
				i++;
				outFilePrefix = argv[i];
				if(outFilePrefix == ""){
					cerr << "[error] Expected a valid filename prefix after '" << cur << "'." << endl;
					exit(1);
				}
			}else{
				cerr << "[error] Expected a valid filename prefix after '" << cur << "'." << endl;
				exit(1);
			}
		}else if(cur == "-h" || cur == "--help" || cur == "-?" || cur == "/?"){
			argc = -1; // print usage
			break;
		}else if(cur == "-g" || cur == "--gzip"){
			isGzipped = true;
		}else if(cur == "-v" || cur == "--verbose"){
			VERBOSE = true;
		}else{
			if(inFileName != ""){
				cerr << programName << ": ambiguous input file, is it '" << inFileName << "' or '" << cur << "'?" << endl;
				exit(1);
			}
			inFileName = cur;
		}
	}
	

	// if no args were given (or user asked for help)
	if(argc < 2){ 
		cerr << "Computes 3D Zernike descriptors a la EM-SURFER." << endl << endl;

		cerr << "Usage: " << programName << " [OPTION]... FILE" << endl << endl;

		cerr << "FILE:" << endl;
		/*cerr << "  If FILE is -, the input map will be read from stdin." << endl;*/
		//cerr << "  Otherwise " << programName << " expects a valid filename. Supported file formats are" << endl;
		cerr << "  " << programName << " expects FILE to be a valid filename. Supported file formats are" << endl;
		cerr << "  Situs, CCP4, and MRC2000. Input may be gzipped. Format is deduced from FILE's" << endl;
		cerr << "  extension. Use the Situs tool map2map to inspect maps that fail to read." << endl << endl;

		cerr << "OPTION:" << endl;
		cerr << "  -n, --order N        Maximum order of invariant to generate. " << endl;
		cerr << "                       Restricted to N in [1..20]. Defaults to 20." << endl;
		cerr << "  -c, --contour C      The level of isosurface to generate invariants for." << endl;
		cerr << "                       C is a real number. Defaults to 1." << endl;
		cerr << "  --contours L         Specify multiple contours." << endl;
		cerr << "                       L is a list of reals delimited by spaces." << endl;
		cerr << "  -p, --prefix S       The output file will be named 'S.inv'. Defaults to FILE." << endl;
		cerr << "                       If S contains '%s', it will be replaced by the contour " << endl;
		cerr << "                       level. If multiple contours are specified but S does not" << endl;
		cerr << "                       contain '%s', the contour will be appended to S." << endl;
		cerr << "  -g, --gzip           Set this flag to force reading input as gzipped." << endl;
		cerr << "  -h, --help, -?, /?   Displays this." << endl << endl;

		cerr << "EXAMPLES:" << endl;
		cerr << "  " << programName << " protein.situs -c 2.75" << endl;
		cerr << "  " << programName << " protein.map.gz -n 11 --contours \"23.2 31.1 45.9\" -p \"protein(%s)\"" << endl << endl;

		return argc==-1 ? 0 : 1;
	}

	// default prefix
	if(outFilePrefix == ""){
		outFilePrefix = inFileName;
	}

	int type = TYPE_UNKNOWN;
	if(endsWith(inFileName, ".gz")){
		isGzipped = true;
		if(endsWith(inFileName, ".map.gz") || endsWith(inFileName, ".ccp4.gz") || endsWith(inFileName, ".mrc.gz")){
			type = TYPE_MRC;
		}else if(endsWith(inFileName, ".situs.gz")){
			type = TYPE_SITUS;
		}else if(endsWith(inFileName, ".grid.gz")){
			type = TYPE_DUMB_GRID;
		}
	}else{
		if(endsWith(inFileName, ".map") || endsWith(inFileName, ".ccp4") || endsWith(inFileName, ".mrc")){
			type = TYPE_MRC;
		}else if(endsWith(inFileName, ".situs")){
			type = TYPE_SITUS;
		}else if(endsWith(inFileName, ".grid")){
			type = TYPE_DUMB_GRID;
		}
	}

	if(type == TYPE_UNKNOWN){
		cerr << "[warning] Failed to deduce map type. Assuming default (Situs)." << endl;
		type = TYPE_SITUS;
	}


	istream *f = NULL;

	try{
		if(isGzipped){
			igzstream *ff = new igzstream();

			ff->exceptions(ios::badbit|ios::failbit);
			ff->open(inFileName.c_str(), ios::in|ios::binary);
			
			f = ff;
		}else{
			ifstream *ff = new ifstream();

			ff->exceptions(ios::badbit|ios::failbit);
			ff->open(inFileName.c_str(), ios::in|ios::binary);

			f = ff;
		}
	}catch(std::exception &e){
		cerr << "[error] Input file " << inFileName << " could not be opened (" << e.what() << ")." <<  endl;
		exit(1);
	}

	try{
		Grid<double> g(*f, type);

		if(VERBOSE) cerr << "[info] Successfully read the input file." << endl;
		
		bool hasPS = false;
		if(outFilePrefix.find("%s") != std::string::npos){
			hasPS = true;
		}
		if(contours.size() == 0){
			Grid<double> pass(g,contour);
			if(VERBOSE) cerr << "[info] Contoured map to " << contourString << "." << endl;
			ZernikeDescriptor<double> zd(pass, maxOrder);
			if(hasPS){
				saveInvariants(zd, formatFilenamePrefix(outFilePrefix, contourString)+".inv");
			}else{
				saveInvariants(zd, outFilePrefix+".inv");
			}
		}else if(contours.size() == 1){
			contour = contours[0];
			Grid<double> pass(g,contour);
			if(VERBOSE) cerr << "[info] Contoured map to " << contourStrings[0] << "." << endl;
			ZernikeDescriptor<double> zd(pass, maxOrder);
			if(hasPS){
				saveInvariants(zd, formatFilenamePrefix(outFilePrefix, contourStrings[0])+".inv");
			}else{
				saveInvariants(zd, outFilePrefix+".inv");
			}
		}else{
			// we have a list
			if(VERBOSE) cerr << "[info] Multiple contours were given." << endl;
			auto cs = contourStrings.begin();
			for(auto &cont : contours){		
				Grid<double> pass(g,cont);
				if(VERBOSE) cerr << "[info] Contoured map to " << *cs << "." << endl;
				ZernikeDescriptor<double> zd(pass, maxOrder);

				if(hasPS){
					saveInvariants(zd, formatFilenamePrefix(outFilePrefix, *cs)+".inv");
				}else{
					saveInvariants(zd, outFilePrefix + *cs + ".inv");
				}
				cs++;
			}
		}
	}catch(bad_alloc &e){
		cerr << "[error] Ran out of memory (" << e.what() << ")." << endl;
		exit(1);
	}catch(runtime_error &e){
		cerr << "[error] An error occurred (" << e.what() << ")." << endl;
		exit(1);
	}catch(exception &e){
		cerr << "[error] An unknown fatal error occurred (" << e.what() << ")." << endl;
		exit(1);
	}

	delete f;

	return 0;
}
