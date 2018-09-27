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

    string id, seq;
    map<string, string> read_seq;
    map<string, string> read_levels;
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
			count++;
			continue;
		} else if(count == 3){
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
        read_seq.insert(make_pair(id.substr(1), seq));
		string level = read_classification.at(id.substr(1)).getClassificationForFileName(group_level);
        read_levels.insert(make_pair(id.substr(1), level));
    }
    cerr << "Outputting Reads" << endl;
    map<string, set<string>> read_groups;
    for(map<string, string>::iterator it = output_handles.begin(); it != output_handles.end(); ++it)
    {   
        set<string> tmp;
        read_groups.insert(make_pair(it->first, tmp));
        cerr << it->first << endl;
        for(map<string, string>::iterator it2 = read_levels.begin(); it2 != read_levels.end(); ++it2){
            if(it2->second == it->first || (it2->second == "" && read_classification.at(it2->first).parentLevel(it->first))){
                read_groups[it->first].insert(it2->first);    
            }
        }
    }
    for(map<string, string>::iterator it = output_handles.begin(); it != output_handles.end(); ++it)
    {   
        gzFile tmp_file = gzopen(it->second.c_str(),"ab");
        for(set<string>::iterator it2 = read_groups[it->first].begin(); it2 != read_groups[it->first].end(); ++it2){
            string val = ">" + *it2 +"\n";
            gzprintf(tmp_file,val.c_str());
            val = read_seq[*it2]+"\n";
            gzprintf(tmp_file,val.c_str());
        }
        gzclose_w(tmp_file);
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
