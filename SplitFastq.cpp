#include <string>
#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <unistd.h>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <regex>
#include <cstdlib>
#include "Read.hpp"
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <zlib.h>

namespace io = boost::iostreams;


using namespace std;

int main(int argc, char** argv)
{
	bool exit = false;
    int opt;
    char group_level = 's';
    string fastqFile, outputFileName, taxonomy_file,line;
    while ((opt = getopt(argc,argv,"i:o:m:x:")) != EOF)
    {
        switch(opt)
        {
            case 'i': fastqFile = optarg; break;
            case 'o': outputFileName = optarg; break;
            case 'm': taxonomy_file = optarg; break;
            case 'x': group_level = optarg[0]; break;
            case '?': exit = true; break;
            default: exit=true;
        }
    }
    if (optind < argc) {
        exit = true;
        cerr << "Non-option argument: ";
        while (optind < argc)
            cerr << argv[optind++];
        cerr << endl;
    }
    if(argc == 1){
        exit = true;
    }
    if(exit){
        cerr << "\nUsage: SplitFastq\n -i InputFastq\n -o Output Path and Prefix\n -m mpa taxonomy file\n (optional) -x group_level [s: species (default), g: genus, f: family, o: order, c: class, p: phylum, d: domain]\n\n";
        return 0;
    }

    cerr << "Reading in MPA File" << endl;
    ifstream inputFile_mpa(taxonomy_file);
    map<string, Read> read_classification;
    set<string> classification_set;
    while (getline(inputFile_mpa, line))
    {   
        istringstream lin(line);
        string id, classification, level;
        lin >> id;
        getline(lin, classification);
        stringstream  data(classification);
        vector<string> result;
        while(getline(data,level,'|'))
        {
            result.push_back(level);
        }
        string tmp = "";
        if(read_classification.count(id) == 0){
        	Read tmp_read = Read(id, 0);
            read_classification.insert(make_pair(id, tmp_read));
            read_classification.at(id).setTaxonomy(result);
            tmp = read_classification.at(id).getClassificationForFileName(group_level);
        }
        if(tmp != ""){
        	classification_set.insert(tmp);
        }
    }

    // Read in the input fastq
    io::filtering_istream in_filter;
    in_filter.push(io::gzip_decompressor());
    in_filter.push(io::file_source(fastqFile));

	//ifstream inputFile_filter(file_name);
    set<string> to_drop;
    cerr << "Reading in fastq" << endl;
    int count = 0;
    int seen = 0;
    // Need a list of output handles to write too
    map<string, string> output_handles;
    // Create a handle for each classificaion in classification_set
    for (set<string>::iterator it = classification_set.begin(); it != classification_set.end(); ++it)
    {
      	output_handles.insert(make_pair(*it, outputFileName + "_" + *it +".fasta.gz"));
    }

    string id, seq, plus, qual;
	while (getline(in_filter, line, '\n'))
	{
		if(count == 0){
			istringstream lin(line);
			// Only want the read ID, not anything else in the header
			lin >> id;
			count++;
			continue;
		} else if(count == 1){
			seq	= line;
			count++;
			continue;
		} else if(count == 2){
			plus = line;
			count++;
			continue;
		} else if(count == 3){
			qual = line;
			count = 0;
			seen++;
		}
		if(seen % 10000 == 0){
			cerr << "Processed " << seen << " reads" << endl;
		}
		// Need to get the species for this read now that we have read it in properly
		// ignore any read that doesn't have a classification
		if(read_classification.count(id.substr(1)) == 0){
			continue;
		}
		string level = read_classification.at(id.substr(1)).getClassificationForFileName(group_level);
		if(level != ""){
			// Classifeid up to the group level
			// Can write to the exact handle
			gzFile tmp_file = gzopen(output_handles[level].c_str(),"ab");
			string val = ">" + id.substr(1)+"\n";
			gzprintf(tmp_file,val.c_str());
			val = seq+"\n";
			gzprintf(tmp_file,val.c_str());
			gzclose_w(tmp_file);

		} else{
			// Need to check which classifications our read should be counted against
			for (map<string, string>::iterator it = output_handles.begin(); it != output_handles.end(); ++it)
    		{
				if(read_classification.at(id.substr(1)).parentLevel(it->first)){
					// Read is a possibly at this classification so write it to file
					gzFile tmp_file = gzopen(it->second.c_str(),"ab");
					string val = ">" + id.substr(1)+"\n";
					gzprintf(tmp_file,val.c_str());
					val = seq+"\n";
					gzprintf(tmp_file,val.c_str());
					gzclose_w(tmp_file);
				}
			}
		}
	}

	cerr << "Completing Output" << endl;
	ofstream tmp;
    tmp.open(outputFileName + "_classification_list.txt");
	for (map<string, string>::iterator it = output_handles.begin(); it != output_handles.end(); ++it)
    {
    	tmp << outputFileName << "_" << it->first << "\n";
    	//string str = "gzip "; 
    	//str = str + it->second; 
    	//const char *command = str.c_str(); 
    	//system(command);
    }
    tmp.close();

	return 0;
}
