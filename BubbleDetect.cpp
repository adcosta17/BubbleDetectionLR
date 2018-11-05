#include <iostream>
#include <unistd.h>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <regex>
#include <cstdlib>
#include <sys/stat.h>
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
    int opt, iterations = 10, fuzz = 1000, threshold = 5;
    string pafFile, outputFileName, taxonomy_file, colour_file, chimeric_read_file, coverage_file;
    while ((opt = getopt(argc,argv,"p:o:i:f:t:r:s:c:h:")) != EOF)
    {
        switch(opt)
        {
            case 'p': pafFile = optarg; break;
            case 'o': outputFileName = optarg; break;
            case 'r': taxonomy_file = optarg; tax = true; break;
            case 's': colour_file = optarg; col = true; break;
            case 'h': chimeric_read_file = optarg; chimeric = true; break;
            case 'i': iterations = atoi(optarg); break;
            case 'f': fuzz = atoi(optarg); break;
            case 't': threshold = atoi(optarg); break;
            case 'c': coverage_file = optarg; coverage = true; break;
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
        cerr << "Usage: BubbleDetect\n -p Paf_Input_File.paf or directory containing paf files\n -o Output Path and Prefix\n (optional) -i Iterations# [10]\n (optional) -f fuzz[1000]\n (optional) -t threshold[5]\n (optional) -r read classification file\n (optional) -s species to colour map\n (optional) -h chimeric read map\n (optional) -c read coverage map\n";
        return 0;
    }

    // Read in and Parse input file
	map<string, vector<Match> > all_matches;
    map<string, vector<Match> > raw_matches;
    map<string, vector<Match> > edge_lists;
    set<string> read_ids;
    map<string, int> read_lengths;
    set<string> chimeric_reads;
    map<string, Read> read_classification;
    map<string, float> read_coverage;

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
            //cout << id << endl;
        }
    }

    if(coverage){
        cerr << "Reading in Read Coverage List" << endl;
        ifstream inputFileCoverage(coverage_file);
        string line;
        while (getline(inputFileCoverage, line, '\n'))
        {
            istringstream lin(line);
            string id;
            int len;
            float cov;
            lin >> id >> len >> cov;
            read_coverage.insert(make_pair(id, cov));
            //cout << read_coverage[id] << endl;
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
        mean_read_length = MatchUtils::read_paf_file(edge_lists, all_matches, raw_matches, read_ids, read_lengths, pafFile, chimeric_reads, read_classification, false);
        cerr << read_ids.size() << " Unique Reads found in File"<< endl;
        cerr << "Average Read Length of " << mean_read_length << " base pairs" << endl;
        // Myers Transitive Reduction Alg
        cerr << "Reducing Edges" << endl;
        MatchUtils::reduce_edges(all_matches, read_ids, edge_lists, fuzz);
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
            rm_edge += MatchUtils::prune_dead_paths(all_matches, read_ids, read_indegree, read_outdegree, de_paths, mean_read_length, 2);
            MatchUtils::clean_matches(all_matches);
        }
        cerr << "\tRemoved " << rm_edge << " Edges" << endl;
        read_indegree.clear();
        read_outdegree.clear();
        MatchUtils::compute_in_out_degree(all_matches, read_ids, read_indegree, read_outdegree);

    } else if(is_dir(pafFile.c_str())){
        cerr << "Parsing Paf Input Directory" << endl;
        MatchUtils::read_and_assemble_paf_dir(all_matches, n50_values, read_ids, read_lengths, pafFile, chimeric_reads, read_classification, fuzz, iterations);
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
    map<string, string> read_full_taxonomy;
    map<string, string> read_lowest_taxonomy;
    map<string, long> classification_count;
    map<string, float> classification_avg_coverage;
    std::ofstream n50Output;
    n50Output.open(outputFileName+"_assembly_stats.txt");

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
            if(classification_count.count(classification) == 0){
                classification_count.insert(make_pair(classification,0));
            }
            classification_count[classification] += read_lengths[id];
            if(col){
                for(map<string, string>::iterator it=species_map.begin(); it!= species_map.end(); ++it) {
                    if (result[result.size() - 1].find(it->first) != std::string::npos) {
                        colours[id] = it->second;
                        break;
                    }
                }
            }
        }
        if(coverage){
            for (map<string, string>::iterator it = read_full_taxonomy.begin(); it != read_full_taxonomy.end(); ++it)
            {
                string c = it->second;
                long n = classification_count[c];
                if(classification_avg_coverage.count(c) == 0){
                    classification_avg_coverage.insert(make_pair(c, 0.0));
                }
                classification_avg_coverage[c] += (read_coverage[it->first]/n);
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
        // First check if Bubble contains 

    	// Four ways of validating.
    	// If we have taxinomic info and coverage, coverage only, taxonomy only and neither of the two
        vector<vector<string> > arms;
        if(tax || coverage){
            MatchUtils::get_bubble_arms((it->first).first, (it->first).second, it->second, read_indegree, read_outdegree, arms);
        }
        if(arms.size() != 2){
            //cerr << (it->first).first << " " << (it->first).second << " " << arms.size() << endl;
            continue;
        }

        //print out bubble arms for valid bubbles
        cerr << endl;
        cerr << "Arms for bubble between " << (it->first).first << " and " << (it->first).second << endl;
        for(int i = 0; i < arms[0].size(); i++){
            cerr << arms[0][i] << " ";
        }
        cerr << endl;
        for(int i = 0; i < arms[1].size(); i++){
            cerr << arms[1][i] << " ";
        }
        cerr << endl;

        vector<float> tax_and_cov; // checks if the coverage for each arm matches the average coverage it should have based on the taxinomic classification of the reads in the arm
        bool tax_only = false; // checks to see if each arm contains at least one unique classification (ideally one read at least in each arm that has a species or subspecies that isnt in the other)
        float cov_only = 0.0; // checks each arm to see if there is a drastic difference in the coverage between them (Possible to detect small errors that cause bubbles by this method as sequencing errors should have lower coverage)
        bool true_bubble = false; // checks to see if the two arms form a true bubble, that is only the start and end nodes have edges to things not in the bubble (two clean arms)
    	if(tax && coverage){
            tax_and_cov = MatchUtils::validBubbleTaxCov(arms, read_coverage, classification_avg_coverage, read_full_taxonomy, read_lengths);
    	}
        if(tax) {
            tax_only = MatchUtils::validBubbleTax(arms, read_lowest_taxonomy);
    	}
        if (coverage){
    		cov_only = MatchUtils::validBubbleCov(arms, read_coverage, read_lengths);
    	}
	    true_bubble =  MatchUtils::check_bubble((it->first).first, (it->first).second, it->second, read_indegree, read_outdegree);
	    float arm_ratio = MatchUtils::getArmLengthRatio(arms, all_matches);

        if(seen_bubbles.count(it->first) == 0 && seen_bubbles.count(std::make_pair((it->first).second, (it->first).first)) == 0){
            
            // Score Bubbles based on values seen
            //Linear:   
            float weights[7] = {0.004047, 0.012653, 0.452049, 0.016768, 0.068084, -0.616781, 0.058613};
            //Logistic:
            //int weights[7] = {0.02297, 0.03667, 3.11104, -0.06563, 0.29392, -4.52549, 0.40195};
            float score = weights[0]*it->second.size();
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
            bubbleOutput << read_names[(it->first).first] << "\t" << read_names[(it->first).second] << "\t" << it->second.size() << "\t";
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
            bubbleOutput << "\t";
            bubbleOutput << cov_only;
            bubbleOutput << "\t";
            bubbleOutput << arm_ratio;
            bubbleOutput << endl;
            
            seen_bubbles.insert(it->first);
        }
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
        rm_edge += MatchUtils::prune_dead_paths(all_matches, read_ids, read_indegree, read_outdegree, de_paths, mean_read_length, 2);
        MatchUtils::clean_matches(all_matches);
    }
    cerr << "\tRemoved " << rm_edge << " Edges" << endl;
    cerr << "Output GFA: "<< outputFileName << ".gfa" << endl;
    read_indegree.clear();
    read_outdegree.clear();
    MatchUtils::compute_in_out_degree(all_matches, read_ids, read_indegree, read_outdegree);
    MatchUtils::toGfa(all_matches,read_lengths, outputFileName+".gfa", read_indegree, read_outdegree, read_names, colours, read_coverage);

    for (int k = 0; k < n50_values.size(); k++)
    {
        n50Output << n50_values[k] << "\n"; 
    }
    n50Output << "Overall\t" << MatchUtils::compute_n50(all_matches, read_indegree, read_outdegree, read_ids) << "\n";
    n50Output.close();
    

  	return 0;
}

