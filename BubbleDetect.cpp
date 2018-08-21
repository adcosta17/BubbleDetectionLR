#include <iostream>
#include <unistd.h>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <regex>
#include <cstdlib>
#include "MatchUtils.hpp"

bool sortSecond(const std::pair<std::string, int>& a, const std::pair<std::string, int>& b){
	return (a.second > b.second);
}

using namespace std;

int main(int argc, char** argv)
{
    bool exit = false;
    bool tax = false;
    bool col = false;
    bool chimeric = false;
    bool group = false;
    bool coverage = false;
    int opt, iterations, fuzz, threshold;
    char group_level = 's';
    string pafFile, outputFileName, taxonomy_file, colour_file, chimeric_read_file, mpa_file, coverage_file;
    while ((opt = getopt(argc,argv,"p:o:i:f:t:r:s:c:x:m:h:")) != EOF)
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
            case 'm': mpa_file = optarg; group = true; break;
            case 'x': group_level = optarg[0]; break;
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
        cerr << "Usage: BubbleDetect\n -p Paf_Input_File.paf\n -o Output Path and Prefix\n -i Iterations#\n -f fuzz\n -t threshold\n (optional) -m mpa taxonomy file\n (optional) -x group by level [s: species (default), g: genus, f: family, o: order, c: class, p: phylum, d: domain]\n (optional) -r read classification file\n (optional) -s species to colour map\n (optional) -h chimeric read map\n (optional) -c read coverage map\n";
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
        ifstream inputFileCoverage(coverage_file);
        string line;
        while (getline(inputFileCoverage, line))
        {
            istringstream lin(line);
            string id;
            float cov;
            lin >> id >> cov;
            read_coverage.insert(make_pair(id, cov));
        }
    }

    cerr << "Parsing Paf Input File" << endl;
    int mean_read_length = MatchUtils::read_paf_file(edge_lists, all_matches, raw_matches, read_ids, read_lengths, pafFile, chimeric_reads, read_classification, true);
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
    map<string, string> read_full_taxonomy;
    map<string, string> read_lowest_taxonomy;
    map<string, int> classification_count;
    map<string, float> classification_avg_coverage;
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
            read_names[id] = read_names[id] + " : " + result[result.size() - 1];
            read_full_taxonomy.insert(make_pair(id, classification));
            read_lowest_taxonomy.insert(make_pair(id, result[result.size() - 1]));
            if(classification_count.count(classification) == 0){
                classification_count.insert(make_pair(classification,0));
            }
            classification_count[classification] += 1;
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
                int n = classification_count[c];
                if(classification_avg_coverage.count(c) == 0){
                    classification_avg_coverage.insert(make_pair(c, 0.0));
                }
                classification_avg_coverage[c] += (read_coverage[it->first]/n);
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

            // Need to read in mpa file to get 
            ifstream inputFile_mpa(mpa_file);
            map<string, int> species_count;
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
                if(read_classification.count(id) != 0){
                    read_classification.at(id).setTaxonomy(result);
                    tmp = read_classification.at(id).getClassification(group_level);
                }
                if(tmp != ""){
                	if(species_count.count(tmp) == 0){
                		species_count.insert(make_pair(tmp, 0));
                	}
                	species_count[tmp] += 1;
                }
            }
            // Convert species_count to vector of pairs, sort by second pair value 
            // Then in decreasing order of reads classified assemble 

          	vector<pair<string,int> > sorted_species_counts;
            for(map<string, int>::iterator it=species_count.begin(); it!= species_count.end(); ++it){
                pair<string, int> tmp;
                tmp.first = it->first;
                tmp.second = it->second;
                sorted_species_counts.push_back(tmp);
            }
            species_count.clear();
           	sort(sorted_species_counts.begin(), sorted_species_counts.end(), sortSecond);


           	//Iterate over the now sorted list of groups, filter reads for each group, assemble them
           	//Prune out any reads that are in a connected component after cleanup, and return rest to pool
            map<string, vector<Match> > graph_edges;
           	set<string> available_reads = read_ids;
            vector<string> n50_values;
           	for(int i=0; i < sorted_species_counts.size(); i++)
            {
           		string level = sorted_species_counts[i].first;
                cerr << "Assembling: " << level << endl;
                set<string> ids_for_group;
                //compute all reads that have this level or also fall into the bins above it
                for(map<string, Read>::iterator it=read_classification.begin(); it!= read_classification.end(); ++it){
                    if(level == it->second.getClassification(group_level)){
                        ids_for_group.insert(it->first);
                    } else if(it->second.parentLevel(level)){
                        // Here Read is subset of level or unclassifed
                        ids_for_group.insert(it->first);
                    }
                }
                // Compute difference between all valid ids and ids that are available
                set<string> ids_to_use;
                set_intersection(available_reads.begin(), available_reads.end(), ids_for_group.begin(), ids_for_group.end(), inserter(ids_to_use, ids_to_use.begin()));
            
                // Now that we have set of ids to use, prune all overlaps to get subset of valid overlaps to use for assembly
                map<string, vector<Match> > species_matches;
                map<string, vector<Match> > species_edge_lists;
                MatchUtils::subset_matches(all_matches, edge_lists, species_matches, species_edge_lists, ids_to_use);
                cerr << "Using " << ids_to_use.size() << " reads" << endl;
                // Myers Transitive Reduction Alg
                cerr << "\tReducing Edges" << endl;
                MatchUtils::reduce_edges(species_matches, ids_to_use, species_edge_lists, fuzz);
                read_indegree.clear();
                read_outdegree.clear();
                MatchUtils::compute_in_out_degree(species_matches, ids_to_use, read_indegree, read_outdegree);
                cerr << "\tDropping Reduced Edges" << endl;
                MatchUtils::clean_matches(species_matches);
                cerr << "\tPruning Dead Ends" << endl;
                for (int j = 0; j < iterations; ++j)
                {
                    set<string> de_ids;
                    map<string, vector<string> > de_paths;
                    read_indegree.clear();
                    read_outdegree.clear();
                    MatchUtils::compute_in_out_degree(species_matches, ids_to_use, read_indegree, read_outdegree);
                    // Prune Dead End Reads
                    MatchUtils::compute_dead_ends(species_matches, ids_to_use,read_indegree, read_outdegree, de_ids, de_paths);
                    MatchUtils::prune_dead_paths(species_matches, ids_to_use, read_indegree, read_outdegree, de_paths, mean_read_length, threshold);
                    MatchUtils::clean_matches(species_matches);
                }
                read_indegree.clear();
                read_outdegree.clear();
                MatchUtils::compute_in_out_degree(species_matches, ids_to_use, read_indegree, read_outdegree);
                string n50_val = MatchUtils::compute_n50(species_matches, read_indegree, read_outdegree, ids_to_use);
                if(n50_val != ""){
                    n50_values.push_back(level + "\n" + n50_val);
                }
                // Need to compute set of reads that was actually used and then remove them from the valid read set
                set<string> reads_in_graph;
                for (set<string>::iterator it=ids_to_use.begin(); it!=ids_to_use.end(); ++it){
                    if(read_outdegree[*it].size() > 0 || read_indegree[*it].size() > 0){
                        // Read is in connected componenet, add to set to drop
                        reads_in_graph.insert(*it);
                    }
                }
                set<string> diff;
                set_difference(available_reads.begin(), available_reads.end(), reads_in_graph.begin(), reads_in_graph.end(), inserter(diff, diff.begin()));
                available_reads = diff;

                // Then we also need to take the overlaps for this species and add it to the master set of overalps to look at that are valid
                // This valid set is what we will call bubbles on
                for (map<string, vector<Match> >::iterator it = species_matches.begin(); it != species_matches.end(); ++it)
                {
                    if(graph_edges.count(it->first) == 0){
                        vector<Match> tmp;
                        graph_edges.insert(pair<string, vector<Match> >(it->first,tmp));
                    }
                    for (int j = 0; j < it->second.size(); ++j)
                    {
                        if(!it->second[j].reduce){
                            graph_edges[it->first].push_back(it->second[j]);
                        }
                    }
                }
           	}
            read_indegree.clear();
            read_outdegree.clear();
            MatchUtils::compute_in_out_degree(graph_edges, read_ids, read_indegree, read_outdegree);
            MatchUtils::toGfa(graph_edges,read_lengths, outputFileName+".gfa", read_indegree, read_outdegree, read_names, colours);

            std::ofstream n50Output;
            n50Output.open(outputFileName+"_assembly_stats.txt");
            for (int k = 0; k < n50_values.size(); k++)
            {
                n50Output << n50_values[k] << "\n"; 
            }
            n50Output << "Overall: " << MatchUtils::compute_n50(all_matches, read_indegree, read_outdegree, read_ids) << "\n";
            n50Output.close();
            return 0;
        }
    }

    if(!group) {

        // Myers Transitive Reduction Alg
        cerr << "Reducing Edges" << endl;
        MatchUtils::reduce_edges(all_matches, read_ids, edge_lists, fuzz);
        MatchUtils::compute_in_out_degree(all_matches, read_ids, read_indegree, read_outdegree);
        
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
        MatchUtils::toGfa(all_matches,read_lengths, outputFileName+".gfa", read_indegree, read_outdegree, read_names, colours);

        std::ofstream n50Output;
        n50Output.open(outputFileName+"_assembly_stats.txt");
        n50Output << "Overall: " << MatchUtils::compute_n50(all_matches, read_indegree, read_outdegree, read_ids) << "\n";
        n50Output.close();
    }

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
    set<pair<string, string> > seen_bubbles;
    for (map<pair<string,string>, set<string> >::iterator it=bubble_sets.begin(); it!=bubble_sets.end(); ++it)
    {    
    	// Four ways of validating.
    	// If we have taxinomic info and coverage, coverage only, taxonomy only and neither of the two
        vector<vector<string> > arms;
        if(tax || coverage){
            MatchUtils::get_bubble_arms((it->first).first, (it->first).second, it->second, read_indegree, read_outdegree, arms);
        }
    	if(tax && coverage){
    		// Can use both the coverage info we have and the taxonomy to get per species/subspecies coverage
    		// And then see if it matches what coverage the arms have
    		bool valid = true;
    		// We now have each arm of the bubble. Now can compute the coverage for each arm, as the average of the coverage for each read in the arm
    		if(arms.size() < 2){
    			continue;
    		}
    		vector<float> avg_cov;
    		for (int i = 0; i < arms.size(); ++i)
    		{
    			float tmp = 0.0;
    			for (int j = 0; j < arms[i].size(); ++j)
    			{
    				tmp += read_coverage[arms[i][j]];
    			}
    			avg_cov.push_back(tmp/arms[i].size());
    		}
    		// Now need to compute the average coverage for each species/supspecies found in the arm
    		// First get the species/subspecies in each arm using read_full_taxonomy
    		// Want to also get the length of each arm as a contig at this spot at well
            vector<float> arm_tax_cov;
    		for (int i = 0; i < arms.size(); ++i)
    		{
                float tmp = 0.0;
    			for (int j = 0; j < arms[i].size(); ++j)
    			{
    				if(j == 0 || j == arms[i].size()){
    					//ignore the first read, this is the start, and ignore the last one, the end. These are shared so they shouldn't be used to distinguish arms
                        continue;
    				}
                    // Otherwise get the edge lengths of each overlap, and then weight each read's coverage. Each read should have a classification. 
                    // If it doesn't then use the read coverage
                    tmp += classification_avg_coverage[read_full_taxonomy[arms[i][j]]];
    			}
                arm_tax_cov.push_back(tmp/arms[i].size());
    		}
            // Now that we have the coverage for each arm, and the average coverage based on the taxonomic classifications of the reads in the arm, we can compare and see if they match what we should be seeing
            // Assumption is that each arm should have similar coverage as the species it comes from. If it does we call it a bubble
            for (int i = 0; i < avg_cov.size(); ++i)
            {
                // Compare the average coverages between the arm itself and the species
                float ratio;
                if(avg_cov[i] > arm_tax_cov[i]){
                    ratio = arm_tax_cov[i]/avg_cov[i];
                } else {
                    ratio = avg_cov[i]/arm_tax_cov[i];
                }
                if(ratio < 0.75){
                    valid = false;
                }
            }
            if(valid){
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
    	} else if(tax) {
    		// Can use the taxonomy information to see if the arms have distinct species or subspecies
            // Don't want to have shared lowest level of taxinomic classification between arms
            // I.E. if each arm has a subspecies classification and share reads that have the same species classification, ok
            // But if they share a subspecies classification, then need to see what percentage of each arm is what subspecies
            // Base it being a bubble on if it has a majority of one subspecies in each arm (> 75%)

            bool valid = false;
            if(arms.size() < 2){
                continue;
            }
            vector<set<string> > arm_classifcation;
            for (int i = 0; i < arms.size(); ++i)
            {
                set<string> tmp;
                for (int j = 0; j < arms[i].size(); ++j)
                {
                    tmp.insert(read_lowest_taxonomy[arms[i][j]]);
                }
                arm_classifcation.push_back(tmp);
            }
            // Now compare the species in each arm to see if there are any that are shared. We should have at least one unique classification in each, ie. distinct sets
            for(int i = 0; i < arm_classifcation.size(); i++){
                for (int j = 0; j < arm_classifcation.size(); ++j)
                {
                    if(i == j){
                        continue;
                    }
                    set<string> res1;
                    set<string> res2;
                    set_difference( arm_classifcation[i].begin(), arm_classifcation[i].end(), arm_classifcation[j].begin(), arm_classifcation[j].end(), inserter(res1, res1.begin()));
                    set_difference( arm_classifcation[j].begin(), arm_classifcation[j].end(), arm_classifcation[i].begin(), arm_classifcation[i].end(), inserter(res2, res2.begin()));
                    if(res1.size() > 0 && res2.size() > 0){
                        valid = true;
                    }
                }
            }
            if(valid){
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
    	} else if (coverage){
    		// Can see if there is drastic differences in coverage between the two arms, Assuming that short bubbles can be caused by a sequencing error
    		// Smaller sequencing errors should have a lower coverage than true variation
    		bool valid = true;
    		// We now have each arm of the bubble. Now can compute the coverage for each arm, as the average of the coverage for each read in the arm
    		if(arms.size() < 2){
    			continue;
    		}
    		vector<float> avg_cov;
    		for (int i = 0; i < arms.size(); ++i)
    		{
    			float tmp = 0.0;
    			for (int j = 0; j < arms[i].size(); ++j)
    			{
    				tmp += read_coverage[arms[i][j]];
    			}
    			avg_cov.push_back(tmp);
    		}
    		for (int i = 0; i < avg_cov.size(); ++i)
    		{
    			for (int j = 0; j < avg_cov.size(); ++j)
    			{
    				// Compare the average coverages between arms, ensure that one isn't super low compared to the other
    				if(i == j) {
    					continue;
    				} 
    				float ratio;
    				if(avg_cov[i] > avg_cov[j]){
    					ratio = avg_cov[j]/avg_cov[i];
    				} else {
    					ratio = avg_cov[i]/avg_cov[j];
    				}
    				if(ratio < 0.25){
    					valid = false;
    				}
    			}
    		}
    		if(valid){
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
    	} else {
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
    }
           
  	return 0;
}

