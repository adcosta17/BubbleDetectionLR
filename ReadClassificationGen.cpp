#include <iostream>
#include <unistd.h>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <regex>
#include <cstdlib>
#include <string>
#include <vector>
#include <map>
#include <set>

bool endsWith(const std::string &mainStr, const std::string &toMatch)
{
	if(mainStr.size() >= toMatch.size() &&
			mainStr.compare(mainStr.size() - toMatch.size(), toMatch.size(), toMatch) == 0)
			return true;
		else
			return false;
}

using namespace std;

string longestCommonSubstring(const string& str1, const string& str2, char delim)
{
	if(str1.empty() || str2.empty())
  	{
    	return 0;
  	}
  	int *curr = new int [str2.size()];
  	int *prev = new int [str2.size()];
  	int *swap = NULL;
  	int maxSubstr = 0;
   	string longest;
  	for(unsigned int i = 0; i < str1.size(); ++i)
  	{
    	for(unsigned int j = 0; j < str2.size(); ++j)
    	{
      		if(str1[i] != str2[j]){
        		curr[j] = 0;
      		} else {
        		if(i == 0 || j == 0){
          			curr[j] = 1;
        		} else {
          			curr[j] = 1 + prev[j-1];
        		}
          		if(maxSubstr < curr[j]) {
          			maxSubstr = curr[j];
             		longest.clear();
        		}
          		if (maxSubstr == curr[j])
          		{
            		longest += str1.substr(i - maxSubstr + 1, i + 1);
          		}
      		}
    	}
    	swap=curr;
    	curr=prev;
    	prev=swap;
  	}
  	delete [] curr;
  	delete [] prev;
  	string val = longest.substr(0, maxSubstr);
  	string d(1, delim);
  	if(endsWith(val, d)){
  		return val.substr(0, val.size()-1);
  	}
  	if(val.find("root") == 0 || val.find("d__") == 0){
  		size_t found = val.find_last_of(d);
  		return val.substr(0, val.size()-(val.size() - found));
  	}
  	if(delim == ';'){
  		return "\troot";
  	} else {
  		return "\t";
  	}
}

int main(int argc, char** argv)
{
	bool exit = false;
	int opt;
	string outputFileName, inputFile, kraken_classificaiton, kraken_mpa_classificaiton;
	while ((opt = getopt(argc,argv,"o:i:k:m:")) != EOF)
    {
        switch(opt)
        {
            case 'i': inputFile = optarg; break;
            case 'o': outputFileName = optarg; break;
            case 'k': kraken_classificaiton = optarg; break;
            case 'm': kraken_mpa_classificaiton = optarg; break;
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
        cerr << "\nUsage: ReadClassificationGen\n -i classification_list.txt\n -o Output Path and Prefix\n -k kraken_classificaiton_map\n -m kraken_mpa_classificaiton_map\n\n";
        return 0;
    }

    map<string, vector<string>> read_classifications_kraken;
    map<string, vector<string>> read_classifications_mpa;
    map<string, string> kcm;
    map<string, string> kmm;

    ifstream krakenClassification(kraken_classificaiton);
    string line;
    while (getline(krakenClassification, line))
    {
        istringstream lin(line);
        string id, classification;
        lin >> id;
        getline(lin, classification);
        if(kcm.count(id) == 0){
        	kcm.insert(make_pair(id, classification));
        }
    }

    ifstream krakenMpaClassification(kraken_mpa_classificaiton);
    while (getline(krakenMpaClassification, line))
    {
        istringstream lin(line);
        string id, classification;
        lin >> id;
        getline(lin, classification);
        if(kmm.count(id) == 0){
        	kmm.insert(make_pair(id, classification));
        }
    }

	ifstream inputFileClassification(inputFile);
    while (getline(inputFileClassification, line))
    {
        istringstream lin(line);
        string id, classification;
        lin >> id >> classification;
        if(read_classifications_kraken.count(id) == 0){
        	vector<string> tmp;
        	read_classifications_kraken.insert(make_pair(id, tmp));
        }
        read_classifications_kraken[id].push_back(kcm[classification]);
        if(read_classifications_mpa.count(id) == 0){
        	vector<string> tmp;
        	read_classifications_mpa.insert(make_pair(id, tmp));
        }
        read_classifications_mpa[id].push_back(kmm[classification]);
    }

    string kfile_name = outputFileName + "_kraken.txt";
    string mpafile_name = outputFileName + "_mpa.txt";
    std::ofstream krakenOutput;
  	krakenOutput.open(kfile_name);
  	std::ofstream mpaOutput;
  	mpaOutput.open(mpafile_name);
    for(map<string, vector<string>>::iterator it = read_classifications_kraken.begin(); it != read_classifications_kraken.end(); ++it){
    	string read = it->first;
    	// For each read see if it has a single classification, or multiple
    	// If single output it to the respective file
    	if(it->second.size() == 1){
    		krakenOutput << read << it->second[0] << "\n";
    		mpaOutput << read << read_classifications_mpa[it->first][0] << "\n";
    		continue;
    	}
    	// Otherwise compute the longest common substring between the classifications
    	string k_string = it->second[0];
    	string m_string = read_classifications_mpa[it->first][0];
    	for(std::size_t i = 1; i < it->second.size(); i++){
    		k_string = longestCommonSubstring(k_string, it->second[i], ';');
    		m_string = longestCommonSubstring(m_string, read_classifications_mpa[it->first][i], '|');
    	}
    	krakenOutput << read << k_string << "\n";
    	mpaOutput << read << m_string << "\n";
    }
    krakenOutput.close();
    mpaOutput.close();


	return 0;
}
