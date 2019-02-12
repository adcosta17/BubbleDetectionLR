#include <iostream>
#include <unistd.h>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <regex>
#include <cstdlib>
#include <sys/stat.h>
#include <dirent.h>
#include "MatchUtils.hpp"

bool sortSecond(const std::pair<std::string, int>& a, const std::pair<std::string, int>& b){
	return (a.second > b.second);
}

bool is_file(const char* path) {
    struct stat buf;
    stat(path, &buf);
    return S_ISREG(buf.st_mode);
}

bool is_dir(const char* path) {
    struct stat buf;
    stat(path, &buf);
    return S_ISDIR(buf.st_mode);
}

using namespace std;

int main(int argc, char** argv)
{
    bool exit = false;
    bool tax = false;
    bool col = false;
    bool chimeric = false;
    bool coverage = false;
    bool mpa = false;
    bool binned = false;
    bool use_ng50 = false;
    bool collapse_contigs = false;
    bool is_binned = false;
    int opt, iterations = 10, fuzz = 1000, threshold = 5, genome_size = 0;
    string pafFile, outputFileName, taxonomy_file, colour_file, chimeric_read_file, coverage_file, mpa_file;
    while ((opt = getopt(argc,argv,"p:o:i:f:t:r:s:c:h:m:g:l:b:")) != EOF)
    {
        switch(opt)
        {
            case 'p': pafFile = optarg; break;
            case 'o': outputFileName = optarg; break;
            case 'r': taxonomy_file = optarg; tax = true; break;
            case 's': colour_file = optarg; col = true; break;
            case 'h': chimeric_read_file = optarg; chimeric = true; break;
            case 'm': mpa_file = optarg; mpa = true; break;
            case 'i': iterations = atoi(optarg); break;
            case 'f': fuzz = atoi(optarg); break;
            case 't': threshold = atoi(optarg); break;
            case 'c': coverage_file = optarg; coverage = true; break;
            case 'g': genome_size = atoi(optarg); use_ng50 = true; break;
            case 'l': collapse_contigs = true; break;
            case 'b': is_binned = true; break;
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
        cerr << "Usage: BubbleDetect\n -p Paf_Input_File.paf or directory containing paf files\n -o Output Path and Prefix\n (optional) -i Iterations# [10]\n (optional) -f fuzz[1000]\n (optional) -t threshold[5]\n (optional) -r read classification file\n (optional) -s species to colour map\n (optional) -h chimeric read map\n (optional) -c read coverage map\n (optional) -m kraken mpa classification file\n (optional) -g total estimated genome size (used for NG50 rather than N50 Calculation)\n";
        return 0;
    }

    /*
    string tmp = "5e6696ec-5b69-4cb7-b0b0-33cccfc6714b";
    cout << tmp << endl;
    tmp = MatchUtils::get_hex_string(tmp);
    cout << tmp << endl;
    cout << "done1" << endl;
    tmp = MatchUtils::get_read_string(tmp);
    cout << tmp << endl;
    cout << "done" << endl;
    */

    // Read in and Parse input file
	map<string, vector<Match> > all_matches;
    map<string, vector<Match> > raw_matches;
    set<string> read_ids;
    map<string, int> read_lengths;
    set<string> chimeric_reads;
    map<string, Read> read_classification;
    map<string, float> read_coverage;
    map<string, string> read_levels;
    map<string, string> read_full_taxonomy;
    map<string, string> read_lowest_taxonomy;
    map<string, float> classification_avg_coverage;
    map<string, float> per_species_coverage;
    map<string, set<string>> read_groups;

    if(chimeric){
        cerr << "Reading in Chimeric Read List" << endl;
        ifstream inputFileChimeric(chimeric_read_file);
        string line;
        while (getline(inputFileChimeric, line))
        {
            istringstream lin(line);
            string id, classification;
            lin >> classification >> id;
            if(classification == "Chimeric"){
                chimeric_reads.insert(id);
            }
        }
    }
    
    if(coverage){
        if(!is_file(coverage_file.c_str()) && !is_dir(coverage_file.c_str())){
            cerr << "Coverage File provided is neither file nor directory" << endl;
            return 0;
        } else if(is_file(coverage_file.c_str())){
            if(mpa){
                // need the mpa to bin reads based on their classification
                cerr << "Reading in MPA File" << endl;
                ifstream inputFile_mpa(mpa_file);
                map<string, Read> mpa_read_classification;
                set<string> classification_set;
                map<string, int> classification_count;
                classification_set.insert("");
                set<string> tmp_set;
                read_groups.insert(make_pair("", tmp_set));
                string line;
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
                    if(mpa_read_classification.count(id) == 0){
                        Read tmp_read = Read(id, 0);
                        mpa_read_classification.insert(make_pair(id, tmp_read));
                        mpa_read_classification.at(id).setTaxonomy(result);
                    }
                    vector<char> groups;
                    groups.push_back('d');
                    groups.push_back('p');
                    groups.push_back('c');
                    groups.push_back('o');
                    groups.push_back('f');
                    groups.push_back('g');
                    groups.push_back('s');
                    read_levels.insert(make_pair(id, tmp));
                    read_groups[tmp].insert(id);
                    for(int k = 0; k < groups.size(); k++){
                        tmp = mpa_read_classification.at(id).getClassificationForFileName(groups[k]);
                        if(tmp != ""){
                            read_levels[id] = tmp;
                            if(read_groups.count(tmp) == 0){
                                read_groups.insert(make_pair(tmp, tmp_set));
                            }
                            read_groups[tmp].insert(id);
                        }
                    }
                }
                cerr << "Reading in Read Coverage List" << endl;
                ifstream inputFileCoverage(coverage_file);
                map<string, int> num_bases;
                map<string, float> tmp_read_coverage;
                while (getline(inputFileCoverage, line, '\n'))
                {
                    istringstream lin(line);
                    string id;
                    int len;
                    float cov;
                    lin >> id >> len >> cov;
                    if(cov < 1){
                    	cov = 1.0;
                    }
                    read_coverage.insert(make_pair(id, cov));
                    tmp_read_coverage.insert(make_pair(id, cov*len));
                    for(map<string, set<string>>::iterator it = read_groups.begin(); it != read_groups.end(); ++it)
                    {
                        if(it->second.count(id) > 0){
                            if(num_bases.count(it->first) == 0){
                                num_bases.insert(make_pair(it->first, 0));
                            }
                            num_bases[it->first] += len;
                        }
                    }
                }
                for(map<string, set<string>>::iterator it = read_groups.begin(); it != read_groups.end(); ++it)
                {
                    per_species_coverage.insert(make_pair(it->first, 0.0));
                    for (map<string, float>::iterator it2 = tmp_read_coverage.begin(); it2 != tmp_read_coverage.end(); ++it2)
                    {
                        per_species_coverage[it->first] += it2->second/num_bases[it->first];
                    }
                }
            }
        } else if(is_dir(coverage_file.c_str())){
            // Need to compute the per read coverage for each species
            binned = true;
            set<string> cov_files;
            DIR *dir;
            struct dirent *ent;
            string suffix = ".read-cov.txt";
            string combined_file = "";
            string combined_suffix = ".combined.read-cov.txt";
            if ((dir = opendir (coverage_file.c_str())) != NULL) {
              /* print all the files and directories within directory */
              while ((ent = readdir (dir)) != NULL) {
                string name(ent->d_name);
                //cerr << name << endl;
                if (name.find(suffix) != string::npos){
                    if(name.find(combined_suffix) != string::npos){
                        combined_file = coverage_file + "/" + name;
                    } else {
                        cov_files.insert(coverage_file + "/" + name);  
                        //cerr << file_name + "/" + name << endl;
                    }
                }
              }
              closedir (dir);
            }
            cerr << "Found " << cov_files.size() << " Coverage Files" << endl;

            // Take the list of paf_files and then for each of them read in the file
            for (set<string>::iterator it = cov_files.begin(); it != cov_files.end(); ++it) {
                map<string, float> tmp_read_coverage;
                int num_bases = 0;
                string cov_fh = *it;
                ifstream inputFileCoverage(cov_fh);
                size_t found = cov_fh.find_last_of("/");
                string name(cov_fh.substr(found+1).substr(0,cov_fh.substr(found+1).size()-11));
                string line;
                while (getline(inputFileCoverage, line, '\n'))
                {
                    istringstream lin(line);
                    string id;
                    int len;
                    float cov;
                    if(cov < 1){
                    	cov = 1.0;
                    }
                    lin >> id >> len >> cov;
                    tmp_read_coverage.insert(make_pair(id, cov*len));
                    num_bases += len;
                    //cout << read_coverage[id] << endl;
                }
                per_species_coverage.insert(make_pair(name, 0.0));
                for (map<string, float>::iterator it2 = tmp_read_coverage.begin(); it2 != tmp_read_coverage.end(); ++it2)
                {
                    per_species_coverage[name] += it2->second/num_bases;
                }
            }
            // Read in combined PAF classification file from same directory and populate read_coverage with it
            cerr << "Reading in Read Coverage List" << endl;
            ifstream inputFileCoverage(combined_file);
            string line;
            while (getline(inputFileCoverage, line, '\n'))
            {
                istringstream lin(line);
                string id;
                int len;
                float cov;
                if(cov < 1){
                	cov = 1.0;
                }
                lin >> id >> len >> cov;
                read_coverage.insert(make_pair(id, cov));
            }
        }
    }

    int mean_read_length = 0;
    map<string, vector<string> > read_indegree; // Number of times read is target
    map<string, vector<string> > read_outdegree;
    vector<string> n50_values;
    
    if(!is_file(pafFile.c_str()) && !is_dir(pafFile.c_str())){
        cerr << "Input provided is neither file nor directory" << endl;
        return 0;
    } else if(is_file(pafFile.c_str())){
        cerr << "Parsing Paf Input File" << endl;
        mean_read_length = MatchUtils::read_paf_file(all_matches, raw_matches, read_ids, read_lengths, pafFile, chimeric_reads, read_classification, false);
        cerr << read_ids.size() << " Unique Reads found in File"<< endl;
        cerr << "Average Read Length of " << mean_read_length << " base pairs" << endl;
        // Myers Transitive Reduction Alg
        cerr << "Reducing Edges" << endl;
        MatchUtils::reduce_edges(all_matches, read_ids, fuzz);
        read_indegree.clear();
        read_outdegree.clear();
        MatchUtils::compute_in_out_degree(all_matches, read_ids, read_indegree, read_outdegree);
        
        //Clean up all_matches and then the edge_lists
        cerr << "Dropping Reduced Edges" << endl;
        MatchUtils::clean_matches(all_matches);
        cerr << "Pruning Dead Ends" << endl;
        int rm_edge = 0;
        for (int i = 0; i < iterations; ++i)
        {
            set<string> de_ids;
            map<string, vector<string> > de_paths;
            read_indegree.clear();
            read_outdegree.clear();
            MatchUtils::compute_in_out_degree(all_matches, read_ids, read_indegree, read_outdegree);
            // Prune Dead End Reads
            MatchUtils::compute_dead_ends(all_matches, read_ids,read_indegree, read_outdegree, de_ids, de_paths);
            rm_edge += MatchUtils::prune_dead_paths(all_matches, read_ids, read_indegree, read_outdegree, de_paths, mean_read_length, threshold);
            MatchUtils::clean_matches(all_matches);
        }
        cerr << "\tRemoved " << rm_edge << " Edges" << endl;
        read_indegree.clear();
        read_outdegree.clear();
        MatchUtils::compute_in_out_degree(all_matches, read_ids, read_indegree, read_outdegree);

    } else if(is_dir(pafFile.c_str()) && is_binned){
        cerr << "Parsing Paf Input Directory" << endl;
        mean_read_length = MatchUtils::read_and_assemble_paf_dir_binned(all_matches, n50_values, read_ids, read_lengths, pafFile, chimeric_reads, read_classification, fuzz, iterations, threshold);
        read_indegree.clear();
        read_outdegree.clear();
        MatchUtils::compute_in_out_degree(all_matches, read_ids, read_indegree, read_outdegree);
    } else if(is_dir(pafFile.c_str())){
        cerr << "Parsing Paf Input Directory" << endl;
        mean_read_length = MatchUtils::read_and_assemble_paf_dir(all_matches, n50_values, read_ids, read_lengths, pafFile, chimeric_reads, read_classification, fuzz, iterations, threshold);
        read_indegree.clear();
        read_outdegree.clear();
        MatchUtils::compute_in_out_degree(all_matches, read_ids, read_indegree, read_outdegree);
    }
    
    map<string, string> read_names;
    for (set<string>::iterator it=read_ids.begin(); it!=read_ids.end(); ++it){
            read_names.insert(make_pair(*it, *it));
    }
    map<string, string> colours;
    for (set<string>::iterator it=read_ids.begin(); it!=read_ids.end(); ++it){
            colours.insert(make_pair(*it, "#bebebe"));// Grey in RGB Hex
    }

    if(tax){
        cerr << "Taxonomy File Detected" << endl;
        map<string, string> species_map;
        string line;
        if(col){
            ifstream inputFile_colour(colour_file);
            while (getline(inputFile_colour, line))
            {
                stringstream data(line);
                string level;
                vector<string> result;
                while(getline(data,level,'\t'))
                {
                    result.push_back(level);
                }
                if(result.size() > 1){
                    species_map[result[0]] = result[1];
                }
            }
            cerr << "Assigning Levels and Colours to Reads" << endl;
        }
        ifstream inputFile_filter(taxonomy_file);
        while (getline(inputFile_filter, line))
        {   
            istringstream lin(line);
            string id, classification, level;
            lin >> id;
            getline(lin, classification);
            stringstream  data(classification);
            vector<string> result;
            while(getline(data,level,';'))
            {
                result.push_back(level);
            }
            if(result.size() > 1){
                read_names[id] = read_names[id] + " : " + result[result.size() - 1];
            }
            read_full_taxonomy.insert(make_pair(id, classification));
            read_lowest_taxonomy.insert(make_pair(id, result[result.size() - 1]));
            if(col){
                for(map<string, string>::iterator it=species_map.begin(); it!= species_map.end(); ++it) {
                    if (result[result.size() - 1].find(it->first) != std::string::npos) {
                        colours[id] = it->second;
                        break;
                    }
                }
            }
        }
    }

    cerr << "Compute Possible Bubbles" << endl;
    //Find Paths to each node from a given start
    map<pair<string,string>, set<string> > bubble_sets;
    for (std::set<std::string>::iterator it=read_ids.begin(); it!=read_ids.end(); ++it)
    {
        //Look for nodes that have at least 2 valid neighbours
        if(read_indegree[*it].size() >= 2){
            // If we have 2 or more valid edges then check if they form a bubble
            // Do a breath first search on the edges out of this node.
            // If we get to a node we have already seen then we have a bubble, and can compute the paths between the start and this node
            //cout << "Checking Branch " << *it << " For Bubble" << endl;
            MatchUtils::find_bubble(*it, read_indegree, read_outdegree, bubble_sets, read_indegree[*it]);
        }
        if(read_outdegree[*it].size() >= 2){
            MatchUtils::find_bubble(*it, read_indegree, read_outdegree, bubble_sets, read_outdegree[*it]);
        }
    }
    // For each candidate bubble check if there is really a bubble between the two
    cerr << "Validate Possible Bubbles & Pop those that are sequencing error based" << endl;
    set<pair<string, string> > seen_bubbles;
    if(tax){
        for (map<pair<string,string>, set<string> >::iterator it=bubble_sets.begin(); it!=bubble_sets.end(); ++it)
        {
            // Need to find any bubbles that are only true bubbles when taxonomy info is present.
            // Idea is that bubbles between regions that are all from the same species or subspecies, with no ambiguity should be popped
            // Ambiguity can occur if each arm of bubble has multiple reads classified to same level but differing classifcation
            // ie. arm 1 is sub-species A and arm2 is sub-species B. Can't choose between them Vs Arm1 is subspecies A and arm2 is just the main species
            if(seen_bubbles.count(it->first) == 0 && seen_bubbles.count(std::make_pair((it->first).second, (it->first).first)) == 0){
                //if((it->first).first == "1807" || (it->first).second == "1807"){
                //    cerr << " Validate " << (it->first).first << " To " << (it->first).second << endl;
                //}
                vector<vector<string> > arms;
                MatchUtils::get_bubble_arms((it->first).first, (it->first).second, it->second, read_indegree, read_outdegree, arms);
                bool tax_only = MatchUtils::validBubbleTax(arms, read_lowest_taxonomy);
                bool true_bubble =  MatchUtils::check_bubble((it->first).first, (it->first).second, it->second, read_indegree, read_outdegree);
                if(true_bubble && !tax_only){
                    // This bubble is a true bubble but all the reads can be from the same species or subspecies
                    // So we can collapse it
                    MatchUtils::collapseBubble(arms, all_matches);
                }
                seen_bubbles.insert(it->first);
            }
        }
    }

    bubble_sets.clear();
    read_indegree.clear();
    read_outdegree.clear();
    MatchUtils::compute_in_out_degree(all_matches, read_ids, read_indegree, read_outdegree);
    for (std::set<std::string>::iterator it=read_ids.begin(); it!=read_ids.end(); ++it)
    {
        //Look for nodes that have at least 2 valid neighbours
        if(read_indegree[*it].size() >= 2){
            MatchUtils::find_bubble(*it, read_indegree, read_outdegree, bubble_sets, read_indegree[*it]);
        }
        if(read_outdegree[*it].size() >= 2){
            MatchUtils::find_bubble(*it, read_indegree, read_outdegree, bubble_sets, read_outdegree[*it]);
        }
    }
    seen_bubbles.clear();
    set<pair<string, string> > collapsed_bubbles;
    if(tax){
        for (map<pair<string,string>, set<string> >::iterator it=bubble_sets.begin(); it!=bubble_sets.end(); ++it)
        {
            // Need to find any bubbles that are only true bubbles when taxonomy info is present.
            // Idea is that bubbles between regions that are all from the same species or subspecies, with no ambiguity should be popped
            // Ambiguity can occur if each arm of bubble has multiple reads classified to same level but differing classifcation
            // ie. arm 1 is sub-species A and arm2 is sub-species B. Can't choose between them Vs Arm1 is subspecies A and arm2 is just the main species
            if(collapsed_bubbles.count(it->first) == 0 && collapsed_bubbles.count(std::make_pair((it->first).second, (it->first).first)) == 0){
                vector<vector<string> > arms;
                MatchUtils::get_bubble_arms((it->first).first, (it->first).second, it->second, read_indegree, read_outdegree, arms);
                bool tax_only = MatchUtils::validBubbleTax(arms, read_lowest_taxonomy);
                bool true_bubble =  MatchUtils::check_bubble((it->first).first, (it->first).second, it->second, read_indegree, read_outdegree);
                if(true_bubble && !tax_only){
                    // This bubble is a true bubble but all the reads can be from the same species or subspecies
                    // So we can collapse it
                    MatchUtils::collapseBubble(arms, all_matches);
                    seen_bubbles.insert(it->first);
                }
                collapsed_bubbles.insert(it->first);
            }
        }
    }
    cerr << "Collapsed Bubbles" << endl;
    // Now that simple bubbles are collapsed, remove internal bubbles
    // These are small bubbles hiding inside larger ones
    for(int i = 0; i < iterations; i++){
	    bubble_sets.clear();
	    read_indegree.clear();
	    read_outdegree.clear();
	    MatchUtils::compute_in_out_degree(all_matches, read_ids, read_indegree, read_outdegree);
	    for (std::set<std::string>::iterator it=read_ids.begin(); it!=read_ids.end(); ++it)
	    {
	        //Look for nodes that have at least 2 valid neighbours
	        if(read_indegree[*it].size() >= 2){
	            MatchUtils::find_bubble(*it, read_indegree, read_outdegree, bubble_sets, read_indegree[*it]);
	        }
	        if(read_outdegree[*it].size() >= 2){
	            MatchUtils::find_bubble(*it, read_indegree, read_outdegree, bubble_sets, read_outdegree[*it]);
	        }
	    }
	    MatchUtils::remove_internal_bubbles(bubble_sets, all_matches, read_indegree, read_outdegree);
        cerr << "Removed internal bubbles round "<< i+1 << endl;
	}

    bubble_sets.clear();
    read_indegree.clear();
    read_outdegree.clear();
    MatchUtils::compute_in_out_degree(all_matches, read_ids, read_indegree, read_outdegree);
    for (std::set<std::string>::iterator it=read_ids.begin(); it!=read_ids.end(); ++it)
    {
        //Look for nodes that have at least 2 valid neighbours
        if(read_indegree[*it].size() >= 2){
            MatchUtils::find_bubble(*it, read_indegree, read_outdegree, bubble_sets, read_indegree[*it]);
        }
        if(read_outdegree[*it].size() >= 2){
            MatchUtils::find_bubble(*it, read_indegree, read_outdegree, bubble_sets, read_outdegree[*it]);
        }
    }
    seen_bubbles.clear();

    // Check to find sets that are unique less the two ends of the pair
    // A true bubble will have reads that only appear in the bubble
    // Will also check that all reads in the set that aren't the two ends don't have any incoming or outgoing edges that are to reads not in the set
    ofstream bubbleOutput;
    bubbleOutput.open(outputFileName+"_bubble_list.txt");
    for (map<pair<string,string>, set<string> >::iterator it=bubble_sets.begin(); it!=bubble_sets.end(); ++it)
    {    
        if((it->first).first == (it->first).second){
            // Loop Bubble, don't want to consider this case
            //cerr << "Bubble has same start and end " << (it->first).first << " and " << (it->first).second << endl;
            seen_bubbles.insert(it->first);
            continue;
        }

    	// Four ways of validating.
    	// If we have taxinomic info and coverage, coverage only, taxonomy only and neither of the two
        vector<vector<string> > arms;
        MatchUtils::get_bubble_arms((it->first).first, (it->first).second, it->second, read_indegree, read_outdegree, arms);
        if(arms.size() < 2){
            //cerr << "Bubble has less than 2 arms " << (it->first).first << " and " << (it->first).second << endl;
            seen_bubbles.insert(it->first);
            continue;
        }

        if(seen_bubbles.count(it->first) != 0 || seen_bubbles.count(make_pair((it->first).second, (it->first).first)) != 0){
            //cerr << "Already Seen " << (it->first).first << " and " << (it->first).second << endl;
            continue;
        }

        for(int i = 0; i<arms.size(); i++){
            for(int j = 1; j < arms.size(); j++){
                if(i >= j){
                    continue;
                }

                vector<string> v1 = arms[i];
                vector<string> v2 = arms[j];
                vector<string> v3;
                sort(v1.begin(), v1.end());
                sort(v2.begin(), v2.end());
                set_intersection(v1.begin(),v1.end(),v2.begin(),v2.end(),back_inserter(v3));
                if(v3.size() != 2){
                    // should have just the start and end reads that are the same, otherwise the bubble are smaller than reported
                    //cerr << "Non Exclusive Arms for " << (it->first).first << " and " << (it->first).second << endl;
                    continue;
                }    
        
                /*
                //print out bubble arms for valid bubbles
                cerr << endl;
                cerr << "Arms for bubble between " << (it->first).first << " and " << (it->first).second << endl;
                for(int k = 0; k < arms[i].size(); k++){
                    cerr << arms[i][k] << " ";
                }
                cerr << endl;
                for(int k = 0; k < arms[j].size(); k++){
                    cerr << arms[j][k] << " ";
                }
                cerr << endl;
                */

                // Figure out which arm is the longer one
                // We will have the longer arm be A and our shorter arm be B
                // Don't count start and end reads
                vector<vector<string> > tmp_arms;
                if(arms[i].size() >= arms[j].size()){
                    tmp_arms.push_back(arms[i]);
                    tmp_arms.push_back(arms[j]);
                } else {
                    tmp_arms.push_back(arms[j]);
                    tmp_arms.push_back(arms[i]);
                }

                vector<float> tax_and_cov; // checks if the coverage for each arm matches the average coverage it should have based on the taxinomic classification of the reads in the arm
                bool tax_only = false; // checks to see if each arm contains at least one unique classification (ideally one read at least in each arm that has a species or subspecies that isnt in the other)
                float cov_only = 0.0; // checks each arm to see if there is a drastic difference in the coverage between them (Possible to detect small errors that cause bubbles by this method as sequencing errors should have lower coverage)
                bool true_bubble = false; // checks to see if the two arms form a true bubble, that is only the start and end nodes have edges to things not in the bubble (two clean arms)
            	if(tax && coverage && (mpa || binned)){
                    tax_and_cov = MatchUtils::validBubbleTaxCov(tmp_arms, read_coverage, per_species_coverage, read_levels, read_lengths, binned, all_matches);
            	}
                if(tax) {
                    tax_only = MatchUtils::validBubbleTax(tmp_arms, read_lowest_taxonomy);
            	}
                if (coverage){
            		cov_only = MatchUtils::validBubbleCov(tmp_arms, read_coverage, read_lengths);
            	}
        	    true_bubble =  MatchUtils::check_bubble((it->first).first, (it->first).second, it->second, read_indegree, read_outdegree);
        	    float arm_ratio = MatchUtils::getArmLengthRatio(tmp_arms, all_matches);
          
                // Score Bubbles based on values seen
                //Linear:   
                float weights[7] = {0.003137, -0.077071, 0.028428, 0.024931, 0.533754, 0.128592, 0.014294};
                //Logistic:
                //int weights[7] = {0.02297, 0.03667, 3.11104, -0.06563, 0.29392, -4.52549, 0.40195};
                float score = 0.111063;
                score += weights[0]*(tmp_arms[0].size()+tmp_arms[1].size());
                if(true_bubble){
                    score += weights[1];
                }
                if(tax_and_cov.size() == 2){
                   score += tax_and_cov[0]*weights[2];
                   score += tax_and_cov[1]*weights[3];
                }
                if(tax_only){
                    score += weights[4];
                }
                score += cov_only*weights[5];
                score += arm_ratio*weights[6];

                bubbleOutput << score << "\t";
                bubbleOutput << tmp_arms[0].size()+tmp_arms[1].size() << "\t";
                if(true_bubble){
                    bubbleOutput << 1;
                } else {
                    bubbleOutput << 0;
                }
                bubbleOutput << "\t";
                if(tax_and_cov.size() == 2){
                    bubbleOutput << tax_and_cov[0] << "\t" << tax_and_cov[1];
                } else {
                    bubbleOutput << 0 << "\t" << 0;
                }
                bubbleOutput << "\t";
                if(tax_only){
                    bubbleOutput << 1;
                } else {
                    bubbleOutput << 0;
                }
                bubbleOutput << "\t" << cov_only << "\t" << arm_ratio << "\t";
                bubbleOutput << read_names[(it->first).first] << "\t" << read_names[(it->first).second] << "\t";
                for(int k = 0; k < tmp_arms[0].size(); k++){
                    bubbleOutput << tmp_arms[0][k] << " ";
                }
                bubbleOutput << "\t";
                for(int k = 0; k < tmp_arms[1].size(); k++){
                    bubbleOutput << tmp_arms[1][k] << " ";
                }
                bubbleOutput << endl;
                
            }
        }
        seen_bubbles.insert(it->first);
    }
    bubbleOutput.close();

    cerr << "Final Clean up" << endl;
    int rm_edge = 0;
    for (int i = 0; i < iterations; ++i)
    {
        set<string> de_ids;
        map<string, vector<string> > de_paths;
        read_indegree.clear();
        read_outdegree.clear();
        MatchUtils::compute_in_out_degree(all_matches, read_ids, read_indegree, read_outdegree);
        // Prune Dead End Reads
        MatchUtils::compute_dead_ends(all_matches, read_ids,read_indegree, read_outdegree, de_ids, de_paths);
        rm_edge += MatchUtils::prune_dead_paths(all_matches, read_ids, read_indegree, read_outdegree, de_paths, mean_read_length, threshold);
        MatchUtils::clean_matches(all_matches);
    }
    cerr << "\tRemoved " << rm_edge << " Edges" << endl;
    cerr << "Output GFA: "<< outputFileName << ".gfa" << endl;
    read_indegree.clear();
    read_outdegree.clear();
    MatchUtils::compute_in_out_degree(all_matches, read_ids, read_indegree, read_outdegree);

    ofstream n50Output;
    n50Output.open(outputFileName+"_assembly_stats.txt");

    for (int k = 0; k < n50_values.size(); k++)
    {
        n50Output << n50_values[k] << "\n"; 
    }
    if(use_ng50){
        n50Output << "Overall\t" << MatchUtils::compute_ng50(all_matches, read_indegree, read_outdegree, read_ids, genome_size, 2*mean_read_length) << "\n";
    } else {
        n50Output << "Overall\t" << MatchUtils::compute_n50(all_matches, read_indegree, read_outdegree, read_ids, 2*mean_read_length) << "\n";
    }
    n50Output.close();

    MatchUtils::toGfa(all_matches,read_lengths, outputFileName+".gfa", read_indegree, read_outdegree, read_names, colours, read_coverage);

    if(collapse_contigs){
        cerr << "Collapsing overlaps to Contigs" << endl;
        MatchUtils::collapse_contigs(all_matches, read_indegree, read_outdegree, read_ids, colours, read_coverage, outputFileName+"_collapsed.gfa", 2*mean_read_length);
    }
    cerr << "Done" << endl;
  	return 0;
}

