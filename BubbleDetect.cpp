#include <iostream>
#include <unistd.h>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <regex>
#include <cstdlib>
#include "MatchUtils.hpp"

using namespace std;


int main(int argc, char** argv)
{
    bool exit = false;
    bool tax = false;
    bool col = false;
    bool chimeric = false;
    bool group = false;
    int opt, iterations, fuzz, threshold;
    string pafFile, outputFileName, taxonomy_file, colour_file, chimeric_read_file;
    while ((opt = getopt(argc,argv,"p:g:i:f:t:r:s:c:x")) != EOF)
    {
        switch(opt)
        {
            case 'p': pafFile = optarg; break;
            case 'g': outputFileName = optarg; break;
            case 'r': taxonomy_file = optarg; tax = true; break;
            case 's': colour_file = optarg; col = true; break;
            case 'c': chimeric_read_file = optarg; chimeric = true; break;
            case 'i': iterations = atoi(optarg); break;
            case 'f': fuzz = atoi(optarg); break;
            case 't': threshold = atoi(optarg); break;
            case 'x': group = true; break;
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
        cerr << "Usage: BubbleDetect\n -p Paf_Input_File.paf\n -g Gfa_Output_File.gfa\n -i Iterations#\n -f fuzz\n -t threshold\n (optional) -x group by species [False]\n (optional) -r read_classification_file\n (optional) -s species_to_colour_map\n (optional) -c chimeric_read_map\n";
        return 0;
    }

    // Read in and Parse input file
	map<string, vector<Match> > all_matches;
    map<string, vector<Match> > raw_matches;
    map<string, vector<Match> > all_matches_reversed;
    map<string, vector<Match> > edge_lists;
    set<string> read_ids;
    map<string, int> read_lengths;
    set<string> chimeric_reads;

    if(chimeric){
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

    cerr << "Parsing Paf Input File" << endl;
    int mean_read_length = MatchUtils::read_paf_file(edge_lists, all_matches, raw_matches, read_ids, read_lengths, pafFile, chimeric_reads, true);
    cerr << read_ids.size() << " Unique Reads found in File"<< endl;
    cerr << "Average Read Length of " << mean_read_length << " base pairs" << endl;
    map<string, vector<string> > read_indegree; // Number of times read is target
    map<string, vector<string> > read_outdegree;
    
    map<string, string> read_names;
    for (set<string>::iterator it=read_ids.begin(); it!=read_ids.end(); ++it){
            read_names.insert(make_pair(*it, *it));
    }
    map<string, string> colours;
    for (set<string>::iterator it=read_ids.begin(); it!=read_ids.end(); ++it){
            colours.insert(make_pair(*it, "#bebebe"));// Grey in RGB Hex
    }
    map<string, vector<string> > read_taxonomy;
    if(tax){
        cerr << "Taxonomy File Detected" << endl;
        map<string, string> species_map;
        if(col){
            ifstream inputFile_colour(colour_file);
            string line;
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
            while(getline(data,level,'|'))
            {
                result.push_back(level);
            }
            read_names[id] = read_names[id] + " : " + result[result.size() - 1];
            read_taxonomy.insert(make_pair(id, result));
            if(col){
                for(map<string, string>::iterator it=species_map.begin(); it!= species_map.end(); ++it){
                    if (result[result.size() - 1].find(it->first) != std::string::npos){
                        colours[id] = it->second;
                        break;
                    }
                }
            }
        }

        if(group){

        	// Default is to group by species, done only if we see taxonomy file
            // Read in Taxonomy file, and store each classification level for each read if possible
            // Can query struct or map for all reads with level equal to x
            // Store list of all seen species 

            // first get a list of all unique speces, and their counts,
            // Sort from least to greatest in terms of count
            // For each in this order assemble             
            map<vector<string>, int> species_count;
            for(map<string, vector<string> >::iterator it=read_taxonomy.begin(); it!= read_taxonomy.end(); ++it){
                if(it->second[it->second.size() - 1].at(0) == 's'){
                    if(species_count.count(it->second) == 0){
                        species_count.insert(make_pair(it->second, 0));
                    }
                    species_count[it->second] += 1;
                }
            }
            
        	return 0;
        }
    }

    // Myers Transitive Reduction Alg
    cerr << "Reducing Edges" << endl;
    MatchUtils::reduce_edges(all_matches, read_ids, edge_lists, fuzz);
    MatchUtils::compute_in_out_degree(all_matches, read_ids, read_indegree, read_outdegree);
    string tmp_name = outputFileName + "reduced";
    MatchUtils::toGfa(all_matches,read_lengths, tmp_name, read_indegree, read_outdegree, read_names, colours);
    
    //Clean up all_matches and then the edge_lists
    cerr << "Dropping Reduced Edges" << endl;
    MatchUtils::clean_matches(all_matches);
    cerr << "Pruning Dead Ends" << endl;
    for (int i = 0; i < iterations; ++i)
    {
        set<string> de_ids;
        map<string, vector<string> > de_paths;
        read_indegree.clear();
        read_outdegree.clear();
        MatchUtils::compute_in_out_degree(all_matches, read_ids, read_indegree, read_outdegree);
        // Prune Dead End Reads
        MatchUtils::compute_dead_ends(all_matches, read_ids,read_indegree, read_outdegree, de_ids, de_paths);
        MatchUtils::prune_dead_paths(all_matches, read_ids, read_indegree, read_outdegree, de_paths, mean_read_length, threshold);
        MatchUtils::clean_matches(all_matches);
    }
    read_indegree.clear();
    read_outdegree.clear();
    MatchUtils::compute_in_out_degree(all_matches, read_ids, read_indegree, read_outdegree);
    tmp_name = outputFileName + "pruned";
    MatchUtils::toGfa(all_matches,read_lengths, tmp_name, read_indegree, read_outdegree, read_names, colours);

    cerr << "Compute Possible Bubbles" << endl;
    //Find Paths to each node from a given start
    map<string, string>  bubble_pairs;
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
        } else if(read_outdegree[*it].size() >= 2){
            //cout << "Checking Branch " << *it << " For Bubble" << endl;
            MatchUtils::find_bubble(*it, read_indegree, read_outdegree, bubble_sets, read_outdegree[*it]);
        }
    }
    // For each candidate bubble check if there is really a bubble between the two
    cerr << "Validate Possible Bubbles" << endl;
    // Check to find sets that are unique less the two ends of the pair
    // A true bubble will have reads that only appear in the bubble
    // Will also check that all reads in the set that aren't the two ends don't have any incoming or outgoing edges that are to reads not in the set
    std::set<std::pair<std::string, std::string> > seen_bubbles;
    for (map<pair<string,string>, set<string> >::iterator it=bubble_sets.begin(); it!=bubble_sets.end(); ++it)
    {    
        if(MatchUtils::check_bubble((it->first).first, (it->first).second, it->second, read_indegree, read_outdegree)){
            // Have Valid Bubble
            if(seen_bubbles.count(it->first) == 0 && seen_bubbles.count(std::make_pair((it->first).second, (it->first).first)) == 0){
                cout << "Bubble Found Between " << read_names[(it->first).first] << " and " << read_names[(it->first).second] << " of size "<< it->second.size() <<" With Nodes: ";
                for (set<string>::iterator it2=it->second.begin(); it2!=it->second.end(); ++it2)
                {   
                    cout << read_names[*it2] << " ";
                }
                cout << endl;
                seen_bubbles.insert(it->first);
            }
        }
    }
    
  	return 0;
}

