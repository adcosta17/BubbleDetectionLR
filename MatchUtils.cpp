#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include <fstream>
#include <algorithm>
#include <list>
#include <utility>
#include <cmath>

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>    
namespace io = boost::iostreams;

#include <dirent.h>

#include "MatchUtils.hpp"

void get_contained_and_chimeric_reads(std::set<std::string>& to_drop, std::set<std::string>& chimeric_reads, std::set<std::string>& read_ids, std::string file_name){
	using namespace std;
    io::filtering_istream in_filter;
    in_filter.push(io::gzip_decompressor());
    in_filter.push(io::file_source(file_name));

	//ifstream inputFile_filter(file_name);
	string line;
	while (getline(in_filter, line, '\n'))
	{
		istringstream lin(line);
    	string c1, c6, meta, cg;
    	char c5;
    	int c2, c3, c4, c7, c8, c9, c10, c11;
    	lin >> c1 >> c2 >> c3 >> c4 >> c5 >> c6 >> c7 >> c8 >> c9 >> c10 >> c11;
        //getline(lin, meta);
        read_ids.insert(c1);
        read_ids.insert(c6);
        if(to_drop.count(c1) > 0 || to_drop.count(c6) > 0 ){
        	continue;
        }
        if(chimeric_reads.count(c1) != 0){
            to_drop.insert(c1);
        }
        if(chimeric_reads.count(c6) != 0){
            to_drop.insert(c6);
            continue;
        }
        Match tmpLine(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,0,0,"");
        if(tmpLine.internal_edge()){
        	continue;
        }
        int contained = tmpLine.check_match_contained();
	    if(contained == 1 && to_drop.count(tmpLine.target_read_id) == 0){
	        to_drop.insert(tmpLine.target_read_id);
            //cout << c6 << "\t" << c7 << "\t" << c8 << "\t" << c9 << "\t" << c5 << "\t" << c1 << "\t" << c2 << "\t" << c3 << "\t" << c4 << endl; 
	    } else if (contained == -1 && to_drop.count(tmpLine.query_read_id) == 0) {
	        to_drop.insert(tmpLine.query_read_id);
            //cout << c1 << "\t" << c2 << "\t" << c3 << "\t" << c4 << "\t" << c5 << "\t" << c6 << "\t" << c7 << "\t" << c8 << "\t" << c9 << endl;
	    }
    }
}

std::vector<int> get_all_matches_for_file(std::map<std::string, std::vector<Match> >& edge_lists, std::map<std::string, std::vector<Match> >& all_matches, std::map<std::string, std::vector<Match> >& raw_matches, std::set<std::string>& read_ids, std::map<std::string, int>& read_lengths, std::string file_name, std::map<std::string, Read>& read_classification, std::set<std::string>& to_drop){
	using namespace std;

	io::filtering_istream in;
    in.push(io::gzip_decompressor());
    in.push(io::file_source(file_name));
    vector<int> sizes;
    string line;
	while (getline(in, line, '\n'))
	{
        // Each line of input is split into columns
        // Each line coresponds to a called overlap
		istringstream lin(line);
    	string c1, c6, meta, cg;
    	char c5;
    	int c2, c3, c4, c7, c8, c9, c10, c11;
    	lin >> c1 >> c2 >> c3 >> c4 >> c5 >> c6 >> c7 >> c8 >> c9 >> c10 >> c11;
        //getline(lin, meta);
        //size_t idx = meta.find("cg:");
        cg = "";
        //if(idx != string::npos){
        //    cg = meta.substr(idx+5);
        //}
        //cerr << cg << endl;
        // Check for self alignments && contained reads
        Match tmpLine(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,0,0,cg);
        if(raw_matches.count(c1) == 0){
            vector<Match> tmp;
            raw_matches.insert(pair<string, vector<Match> >(c1,tmp));
        }
        if(raw_matches.count(c6) == 0){
            vector<Match> tmp;
            raw_matches.insert(pair<string, vector<Match> >(c6,tmp));
        }
        if(c1 < c6) {
            raw_matches[c1].push_back(tmpLine);
        } else {
            raw_matches[c6].push_back(tmpLine);
        }
        if(c1 == c6 || to_drop.count(c1) > 0 || to_drop.count(c6) > 0) {
            continue;
        }
        // Match is not long enough
        if((c4 - c3) < 2000 || (c9 - c8) < 2000 || c10 < 100){
        	//cerr << c4 - c3 << " " << c8 - c7  << endl;
        	continue;
        }
        // Check that the alignments are proper, and at the ends of reads, if not edge reduction will fail    	
        // Generate Unflitered Matches for the raw gfa file
        if(!tmpLine.internal_edge() && tmpLine.check_match_contained() == 0){
        	if(all_matches.count(c1) == 0){
        		vector<Match> tmp;
        		all_matches.insert(pair<string, vector<Match> >(c1,tmp));
        	}
        	if(all_matches.count(c6) == 0){
        		vector<Match> tmp;
        		all_matches.insert(pair<string, vector<Match> >(c6,tmp));
        	}
	        // Determine which list to store under based on lexographic comparison
	        // Should never have equality here, check done above already 
	        if(c1 < c6) {
	        	all_matches[c1].push_back(tmpLine);
	        } else {
	        	all_matches[c6].push_back(tmpLine);
	        }
	        if(edge_lists.count(c1) == 0){
	        	vector<Match> tmp;
	        	edge_lists.insert(pair<string, vector<Match> >(c1,tmp));
	        }
	        if(edge_lists.count(c6) == 0){
	        	vector<Match> tmp;
	        	edge_lists.insert(pair<string, vector<Match> >(c6,tmp));
	        }
	        tmpLine.cigar = "";
	        edge_lists[c1].push_back(tmpLine);
	        Match tmpLine2(c6,c7,c8,c9,c5,c1,c2,c3,c4,c10,c11,0,0,"");
	        tmpLine2.check_match_contained();
	        edge_lists[c6].push_back(tmpLine2);

	        read_ids.insert(c1);
	        read_ids.insert(c6);
            if(read_classification.count(c1) == 0){
                Read tmp = Read(c1, c2);
                read_classification.insert(make_pair(c1, tmp));
            }
            if(read_classification.count(c6) == 0){
                Read tmp = Read(c6, c7);
                read_classification.insert(make_pair(c6, tmp));
            }
	        if(read_lengths.count(c1) == 0){
	            read_lengths.insert(pair<string, int>(c1, c2));
                sizes.push_back(c2);
	        }
	        if(read_lengths.count(c6) == 0){
	            read_lengths.insert(pair<string, int>(c6, c7));
                sizes.push_back(c7);
	        }
    	}
    }
    return sizes;
}

int MatchUtils::read_paf_dir(std::map<std::string, std::vector<Match> >& edge_lists, std::map<std::string, std::vector<Match> >& all_matches, std::map<std::string, std::vector<Match> >& raw_matches, std::set<std::string>& read_ids, std::map<std::string, int>& read_lengths, std::string file_name, std::set<std::string>& chimeric_reads, std::map<std::string, Read>& read_classification){
	using namespace std;
	set<string> paf_files;
	DIR *dir;
	struct dirent *ent;
	string suffix = ".ava.paf.gz";
	if ((dir = opendir (file_name.c_str())) != NULL) {
	  /* print all the files and directories within directory */
	  while ((ent = readdir (dir)) != NULL) {
	  	string name(ent->d_name);
	  	//cerr << name << endl;
	  	if (name.find(suffix) != string::npos){
	    	paf_files.insert(file_name + "/" + name);
	    	//cerr << file_name + "/" + name << endl;
        }
	  }
	  closedir (dir);
	} else {
	  /* could not open directory */
	  cerr << "ERROR: Could not open " << file_name << endl; 
	  return 0;
	}
	cerr << "Found " << paf_files.size() << " Paf Files" << endl;

	// Take the list of paf_files and then for each of them read in the file
	set<string> to_drop;
	cerr << "Caculating Contained Reads" << endl;
	int count = 0;
	for (set<string>::iterator it = paf_files.begin(); it != paf_files.end(); ++it) {
		get_contained_and_chimeric_reads(to_drop, chimeric_reads, read_ids, *it);
		count++;
		if(count % 20 == 0){
			cerr << "Processed " << count << " Pafs" << endl;
		}
	}
	count = 0;
	read_ids.clear();
	// Now that we have the contained read, read in the files again and then merge the contents after every pass into our complete set
	cerr << "Reading in Valid Matches and merging overlap sets" << endl; 
	for (set<string>::iterator it = paf_files.begin(); it != paf_files.end(); ++it) {
        map<string, vector<Match>> tmp_edge_lists;
        map<string, vector<Match>> tmp_all_matches;
        map<string, vector<Match>> tmp_raw_matches;
        set<string> tmp_read_ids;
        count++;
		if(count % 20 == 0){
			cerr << "Processed " << count << " Pafs" << endl;
		}
		get_all_matches_for_file(tmp_edge_lists, tmp_all_matches, tmp_raw_matches, tmp_read_ids, read_lengths, *it, read_classification, to_drop);
		// need to combine out tmp sets with the actual sets, ensuring to not add duplicate matches

        for(set<string>::iterator it = tmp_read_ids.begin(); it != tmp_read_ids.end(); ++it){
            if(read_ids.count(*it) == 0){
                read_ids.insert(*it);
            }
        }

		for (map<string, vector<Match>>::iterator it = tmp_edge_lists.begin(); it != tmp_edge_lists.end(); ++it)
		{
			if(edge_lists.count(it->first) == 0){
				edge_lists.insert(make_pair(it->first, it->second));
				continue;
			}
			for (int i = 0; i < it->second.size(); ++i)
			{
				// Check the edge list for just the query read
                // Because of the way we add edges we know that if the query reads vector has the edge then the target will have it also and vice versa
                string q = it->second[i].query_read_id;
                string t = it->second[i].target_read_id;
                // Its possible that we haven't seen this read before
                // Since we could have an edge that is between something we have seen and something we haven't
                //  q might be the thing we haven't seen, which means that we will have HAVE to have seen t already
                //  based on the if statement check above
                bool found = false;
                if(edge_lists.count(q) > 0){
                    for(int j = 0; j < edge_lists[q].size(); j++){
                        if((edge_lists[q][j].query_read_id == q && edge_lists[q][j].target_read_id == t) ||
                            (edge_lists[q][j].query_read_id == t && edge_lists[q][j].target_read_id == q)) {
                            found = true;
                            break;
                        }
                    }
                    if(!found){
                        // Didn't find an entry so we should add ours to the respective lists
                        edge_lists[q].push_back(it->second[i]);
                        if(edge_lists.count(t) == 0){
                            vector<Match> tmp;
                            edge_lists.insert(make_pair(t, tmp));
                        }
                        edge_lists[t].push_back(it->second[i]);
                    }
                } else if(edge_lists.count(t) > 0){
                    for(int j = 0; j < edge_lists[t].size(); j++){
                        if((edge_lists[t][j].query_read_id == q && edge_lists[t][j].target_read_id == t) ||
                            (edge_lists[t][j].query_read_id == t && edge_lists[t][j].target_read_id == q)) {
                            found = true;
                            break;
                        }
                    }
                    if(!found){
                        // Didn't find an entry so we should add ours to the respective lists
                        edge_lists[t].push_back(it->second[i]);
                        if(edge_lists.count(q) == 0){
                            vector<Match> tmp;
                            edge_lists.insert(make_pair(q, tmp));
                        }
                        edge_lists[q].push_back(it->second[i]);
                    }
                } else {
                    // We shouldn't get here
                }
			}
		}

        for (map<string, vector<Match>>::iterator it = tmp_raw_matches.begin(); it != tmp_raw_matches.end(); ++it)
        {
            if(raw_matches.count(it->first) == 0){
                raw_matches.insert(make_pair(it->first, it->second));
                continue;
            }
            for (int i = 0; i < it->second.size(); ++i)
            {
                // Check the edge list for just the query read
                // Because of the way we add edges we know that if the query reads vector has the edge then the target will have it also and vice versa
                string q = it->second[i].query_read_id;
                string t = it->second[i].target_read_id;
                bool found = false;
                if(q < t){
                    // Edge should be indexed via q
                    if(raw_matches.count(q) > 0){
                        for (int j = 0; j < raw_matches[q].size(); ++j)
                        {
                            if((raw_matches[q][j].query_read_id == q && raw_matches[q][j].target_read_id == t) ||
                                (raw_matches[q][j].query_read_id == t && raw_matches[q][j].target_read_id == q)) {
                                found = true;
                                break;
                            }
                        }
                    } else {
                        vector<Match> tmp;
                        raw_matches.insert(make_pair(q, tmp));
                    }
                    if(!found){
                        raw_matches[q].push_back(it->second[i]);
                    }
                } else {
                    // Edge should be indexed via t
                    if(raw_matches.count(t) > 0){
                        for (int j = 0; j < raw_matches[t].size(); ++j)
                        {
                            if((raw_matches[t][j].query_read_id == q && raw_matches[t][j].target_read_id == t) ||
                                (raw_matches[t][j].query_read_id == t && raw_matches[t][j].target_read_id == q)) {
                                found = true;
                                break;
                            }
                        }
                    } else {
                        vector<Match> tmp;
                        raw_matches.insert(make_pair(t, tmp));
                    }
                    if(!found){
                        raw_matches[t].push_back(it->second[i]);
                    }
                }
            }
        }

        for (map<string, vector<Match>>::iterator it = tmp_all_matches.begin(); it != tmp_all_matches.end(); ++it)
        {
            if(all_matches.count(it->first) == 0){
                all_matches.insert(make_pair(it->first, it->second));
                continue;
            }
            for (int i = 0; i < it->second.size(); ++i)
            {
                // Check the edge list for just the query read
                // Because of the way we add edges we know that if the query reads vector has the edge then the target will have it also and vice versa
                string q = it->second[i].query_read_id;
                string t = it->second[i].target_read_id;
                bool found = false;
                if(q < t){
                    // Edge should be indexed via q
                    if(all_matches.count(q) > 0){
                        for (int j = 0; j < all_matches[q].size(); ++j)
                        {
                            if((all_matches[q][j].query_read_id == q && all_matches[q][j].target_read_id == t) ||
                                (all_matches[q][j].query_read_id == t && all_matches[q][j].target_read_id == q)) {
                                found = true;
                                break;
                            }
                        }
                    } else {
                        vector<Match> tmp;
                        all_matches.insert(make_pair(q, tmp));
                    }
                    if(!found){
                        all_matches[q].push_back(it->second[i]);
                    }
                } else {
                    // Edge should be indexed via t
                    if(all_matches.count(t) > 0){
                        for (int j = 0; j < all_matches[t].size(); ++j)
                        {
                            if((all_matches[t][j].query_read_id == q && all_matches[t][j].target_read_id == t) ||
                                (all_matches[t][j].query_read_id == t && all_matches[t][j].target_read_id == q)) {
                                found = true;
                                break;
                            }
                        }
                    } else {
                        vector<Match> tmp;
                        all_matches.insert(make_pair(t, tmp));
                    }
                    if(!found){
                        all_matches[t].push_back(it->second[i]);
                    }
                }
            }
        }
	}


}

int MatchUtils::read_paf_file(std::map<std::string, std::vector<Match> >& edge_lists, std::map<std::string, std::vector<Match> >& all_matches, std::map<std::string, std::vector<Match> >& raw_matches, std::set<std::string>& read_ids, std::map<std::string, int>& read_lengths, std::string file_name, std::set<std::string>& chimeric_reads, std::map<std::string, Read>& read_classification, bool gfa)
{
	using namespace std;
    io::filtering_istream in_filter;
    in_filter.push(io::gzip_decompressor());
    in_filter.push(io::file_source(file_name));

	//ifstream inputFile_filter(file_name);
	string line;
    set<string> to_drop;
    int count = 0;
    cerr << "Checking for Contained Reads" << endl;
    get_contained_and_chimeric_reads(to_drop, chimeric_reads, read_ids, file_name);
	cerr << "Found: " << to_drop.size() << " Contained Reads and "<< read_ids.size() << " total reads" << endl;
	cerr << "Reading in all valid Matches" << endl;

	read_ids.clear();
	vector<int> sizes = get_all_matches_for_file(edge_lists, all_matches, raw_matches, read_ids, read_lengths, file_name, read_classification, to_drop);

    return accumulate( sizes.begin(), sizes.end(), 0)/sizes.size();
}

void MatchUtils::collapseBubble(std::vector<std::vector<std::string> >& arms, std::map<std::string, std::vector<Match> >& all_matches){
    // Keep the first arm and remove the rest
    if(arms.size() < 2){
    	//std::cout << "Not enough arms" << std::endl;
        return;
    }
    for (int i = 1; i < arms.size(); ++i)
    {
        //Remove this arm
        for (int j = 0; j < arms[i].size()-1; ++j)
        {
            // Edges are ordered so we should be able to find the respective edge in all matches and remove it
            std::string first = arms[i][j];
            std::string second = arms[i][j+1];
            if(first < second){
                for(int k = 0; k < all_matches[first].size(); k++){
                    if(all_matches[first][k].query_read_id == second || all_matches[first][k].target_read_id == second){
                        all_matches[first][k].reduce = true;
                        //std::cout << "Reduced: "<< all_matches[first][k].query_read_id << " to " << all_matches[first][k].target_read_id << std::endl;
                        break;
                    }
                }
            } else {
                for(int k = 0; k < all_matches[second].size(); k++){
                    if(all_matches[second][k].query_read_id == first || all_matches[second][k].target_read_id == first){
                        all_matches[second][k].reduce = true;
                        //std::cout << "Reduced: "<< all_matches[second][k].query_read_id << " to " << all_matches[second][k].target_read_id << std::endl;
                        break;
                    }
                }
            }
        }
    }
}

bool MatchUtils::validBubbleCov(std::vector<std::vector<std::string> >& arms, std::map<std::string, float>& read_coverage){
    // Can see if there is drastic differences in coverage between the two arms, Assuming that short bubbles can be caused by a sequencing error
    // Smaller sequencing errors should have a lower coverage than true variation
    bool valid = true;
    // We now have each arm of the bubble. Now can compute the coverage for each arm, as the average of the coverage for each read in the arm
    std::vector<float> avg_cov;
    for (int i = 0; i < arms.size(); ++i)
    {
        float tmp = 0.0;
        for (int j = 0; j < arms[i].size(); ++j)
        {
            tmp += read_coverage[arms[i][j]];
        }
        avg_cov.push_back(tmp/arms[i].size());
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
    return valid;
}

bool MatchUtils::validBubbleTax(std::vector<std::vector<std::string> >& arms, std::map<std::string, std::string>& read_lowest_taxonomy){
    // Can use the taxonomy information to see if the arms have distinct species or subspecies
    // Don't want to have shared lowest level of taxinomic classification between arms
    // I.E. if each arm has a subspecies classification and share reads that have the same species classification, ok
    // But if they share a subspecies classification, then need to see what percentage of each arm is what subspecies
    // Base it being a bubble on if it has a majority of one subspecies in each arm (> 75%)
    using namespace std;

    bool valid = false;
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
    return valid;
}

bool MatchUtils::validBubbleTaxCov(std::vector<std::vector<std::string> >& arms, std::map<std::string, float>& read_coverage, std::map<std::string, float>& classification_avg_coverage, std::map<std::string, std::string>& read_full_taxonomy){
    // Can use both the coverage info we have and the taxonomy to get per species/subspecies coverage
    // And then see if it matches what coverage the arms have
    bool valid = true;
    // We now have each arm of the bubble. Now can compute the coverage for each arm, as the average of the coverage for each read in the arm

    std::vector<float> avg_cov;
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
    std::vector<float> arm_tax_cov;
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
    return valid;
}

std::string MatchUtils::compute_n50(std::map<std::string, std::vector<Match> >& all_matches, std::map<std::string,std::vector<std::string> >& read_indegree, std::map<std::string,std::vector<std::string> >& read_outdegree, std::set<std::string>& read_ids){
    // Need to evaluate the assembly stats (N50, Overall length, and total number of contigs.)
    // Don't want to fully visualize as we need to see where each read comes from
    // Last step of Layout phase in OLC
    // Look for nodes that are branches. Iterate along each branch and add overlaps to list of edges against respective contig
    // Keep going until we get to a node that is a branch, in either direction.

    using namespace std;

    map<int, vector<Match> > contigs_sets;
    int contig_num = 0;
    for (set<string>::iterator it=read_ids.begin(); it!=read_ids.end(); ++it)
    {
        // Look for nodes that have at least 2 valid neighbours, or are dead ends
        // For each one we need to iterate over each neighbour and find path to next branch node
        if(read_indegree[*it].size() >= 2 || read_outdegree[*it].size() >= 2 || (read_indegree[*it].size() == 0 && read_outdegree[*it].size() >= 1) || (read_outdegree[*it].size() == 0 && read_indegree[*it].size() >= 1)) {
            contig_num = MatchUtils::compute_contigs(*it, all_matches, read_indegree, read_outdegree, contigs_sets, contig_num);
        }
    }
    if(contig_num == 0){
        // No contigs found
        return "";
    }
    vector<int> contig_lengths;
    int total_length = 0;
    // Now compute the length of each contig
    for (map<int, vector<Match> >::iterator it = contigs_sets.begin(); it != contigs_sets.end(); ++it)
    {
        if(it->second.size() == 1){
            if(it->second[0].length_to_use == 0){
                continue;
            }
        }
        int len = 0;
        for (int i = 0; i < it->second.size(); ++i)
        {
            if(i == 0){
                // Know that we must have at least two edges, otherwise they wouldn't be included here
                if(it->second[i].length_to_use > 0){
                    len = it->second[i].length_to_use;
                    break;
                }
                else if(it->second[i].query_read_id == it->second[i+1].query_read_id || it->second[i].query_read_id == it->second[i+1].target_read_id){
                    len = it->second[i].target_read_length;
                } else {
                    len = it->second[i].query_read_length;
                }
            }
            len += it->second[i].length;
        }
        if(it->second.size() > 1){
            //cerr << "Contig Between: " << len << " " << it->second[0].query_read_id << " : " << it->second[0].target_read_id << " and " << it->second[it->second.size() -1].query_read_id << " : " << it->second[it->second.size() -1].target_read_id << endl;
        }
        contig_lengths.push_back(len);
        total_length += len;
    }

    // Compute N50 and report number of contigs and total length
    sort(contig_lengths.begin(), contig_lengths.end(), greater<int>());
    int n50 = 0;
    int sum_len = 0;
    for (int i = 0; i < contig_lengths.size(); ++i)
    {
        sum_len += contig_lengths[i];
        //cout << sum_len << " " << contig_lengths[i] << " " << total_length << " " << sum_len/float(total_length) << endl;
        if(sum_len/float(total_length) > 0.5){
            n50 = contig_lengths[i];
            break;
        }
    }

    return string("") + to_string(n50) + "\t" + to_string(total_length) + "\t" + to_string(contig_lengths.size()) + "\t" + to_string(contig_lengths[0]);
}


void MatchUtils::subset_matches(std::map<std::string, std::vector<Match> >& all_matches, std::map<std::string, std::vector<Match> >& edge_lists, std::map<std::string, std::vector<Match> >& species_matches, std::map<std::string, std::vector<Match> >&  species_edge_lists, std::set<std::string> ids_to_use){
    for (std::map<std::string, std::vector<Match> >::iterator it=all_matches.begin(); it!=all_matches.end(); ++it)
    {
        for (int i = 0; i < it->second.size(); ++i)
        {
            if(ids_to_use.count(it->second[i].query_read_id) > 0 && ids_to_use.count(it->second[i].target_read_id) > 0){
                if(species_matches.count(it->first) == 0){
                    std::vector<Match> tmp;
                    species_matches.insert(std::pair<std::string, std::vector<Match> >(it->first,tmp));
                }
                species_matches[it->first].push_back(it->second[i]);
            }
        }
    }

    for (std::map<std::string, std::vector<Match> >::iterator it=edge_lists.begin(); it!=edge_lists.end(); ++it)
    {
        for (int i = 0; i < it->second.size(); ++i)
        {
            if(ids_to_use.count(it->second[i].query_read_id) > 0 && ids_to_use.count(it->second[i].target_read_id) > 0){
                if(species_edge_lists.count(it->first) == 0){
                    std::vector<Match> tmp;
                    species_edge_lists.insert(std::pair<std::string, std::vector<Match> >(it->first,tmp));
                }
                species_edge_lists[it->first].push_back(it->second[i]);
            }
        }
    }
}


void MatchUtils::find_bubble(std::string start, std::map<std::string,std::vector<std::string> >& read_indegree, std::map<std::string,std::vector<std::string> >& read_outdegree, std::map<std::pair<std::string,std::string>, std::set<std::string> >& bubble_sets, std::vector<std::string>& start_ids){

    std::set<std::string> visited;
    std::list<std::pair<std::string,std::string> > q;
    std::pair<std::string,std::string> end = std::make_pair("","");
    visited.insert(start);
    bool print = false;
    //if(start == "60275"){
    //    print = true;
    //}
    for (int i = 0; i < start_ids.size(); ++i)
    {   
        visited.insert(start_ids[i]);
        q.push_back(std::make_pair(start_ids[i], start));
    }
    while(!q.empty())
    {
        std::pair<std::string,std::string> s = q.front();
        if(print){
            std::cout << "At: " << s.first << " From " << s.second << std::endl;
        }
        q.pop_front();
        bool check_outdegree = false;
        if(std::find(read_indegree[s.first].begin(), read_indegree[s.first].end(), s.second) != read_indegree[s.first].end()){
            check_outdegree = true;
        }
        if(check_outdegree){
            for (int i = 0; i < read_outdegree[s.first].size(); ++i)
            {
                if(visited.count(read_outdegree[s.first][i]) == 0){
                    visited.insert(read_outdegree[s.first][i]);
                    q.push_back(std::make_pair(read_outdegree[s.first][i], s.first));
                } else {
                    // Found a node we have already seen, possible bubble
                    end = std::make_pair(read_outdegree[s.first][i],s.first);
                    q.clear();
                    break;
                }
            }
        } else {
            for (int i = 0; i < read_indegree[s.first].size(); ++i)
            {
                if(visited.count(read_indegree[s.first][i]) == 0){
                    visited.insert(read_indegree[s.first][i]);
                    q.push_back(std::make_pair(read_indegree[s.first][i], s.first));
                } else {
                    // Found a node we have already seen, possible bubble
                    end = std::make_pair(read_indegree[s.first][i],s.first);
                    q.clear();
                    break;
                }
            }
        }
    }

    // Once we have a node that we have already seen, backtrack and compute the set of visited nodes via a BFS from node we've already seen to the start node
    // Know that 2 paths exist, as proven in the first half so we should be able to find them.
    // Can take the two visited sets and do a intersection, giving us the set of nodes that is on the path between the two exactly

    if(end.first != ""){
        if(print){
            std::cout << "Possible Bubble between: " << start << " and " << end.first << std::endl;
        }
        q.clear();
        std::set<std::string> visited_back;
        visited_back.insert(end.first);
        // Need to check on the direction of our end node, see what the predacesor was and if it is in indegreee or outdegree
        if(std::find(read_indegree[end.first].begin(), read_indegree[end.first].end(), end.second) != read_indegree[end.first].end()){
            for (int i = 0; i < read_indegree[end.first].size(); ++i)
            {   
                visited_back.insert(read_indegree[end.first][i]);
                q.push_back(std::make_pair(read_indegree[end.first][i], end.first));
            }
        } else {
            for (int i = 0; i < read_outdegree[end.first].size(); ++i)
            {   
                visited_back.insert(read_outdegree[end.first][i]);
                q.push_back(std::make_pair(read_outdegree[end.first][i], end.first));
            }
        }
        std::string bubble_end = start;
        while(!q.empty())
        {
            std::pair<std::string,std::string> e = q.front();
            if(print){
                std::cout << "At: " << e.first << " From " << e.second << std::endl;
            }
            q.pop_front();
            bool check_outdegree = false;
            if(std::find(read_indegree[e.first].begin(), read_indegree[e.first].end(), e.second) != read_indegree[e.first].end()){
                check_outdegree = true;
            }
            if(check_outdegree){
                for (int i = 0; i < read_outdegree[e.first].size(); ++i)
                {
                    if(visited_back.count(read_outdegree[e.first][i]) == 0){
                        visited_back.insert(read_outdegree[e.first][i]);
                        q.push_back(std::make_pair(read_outdegree[e.first][i], e.first));
                    } else {
                        // Found a node we have already seen, possible bubble
                        bubble_end = read_outdegree[e.first][i];
                        q.clear();
                        break;
                    }
                }
            } else {
                for (int i = 0; i < read_indegree[e.first].size(); ++i)
                {
                    if(visited_back.count(read_indegree[e.first][i]) == 0){
                        visited_back.insert(read_indegree[e.first][i]);
                        q.push_back(std::make_pair(read_indegree[e.first][i], e.first));
                    } else {
                        // Found a node we have already seen, possible bubble
                        bubble_end = read_indegree[e.first][i];
                        q.clear();
                        break;
                    }
                }
            }
        }
        std::set<std::string> combined;
        std::set_intersection(visited.begin(), visited.end(), visited_back.begin(), visited_back.end(), std::inserter(combined, combined.begin()));
        bubble_sets.insert(std::pair<std::pair<std::string,std::string>, std::set<std::string>>(make_pair(bubble_end, end.first), combined));
    } 
}

void recurse_bubble_arm(std::string id, std::string end, std::set<std::string> reads, std::map<std::string,std::vector<std::string> >& read_indegree, std::map<std::string,std::vector<std::string> >& read_outdegree, std::vector<std::string>& tmp){
	if(std::find(tmp.begin(), tmp.end(), id) == tmp.end()){
		tmp.push_back(id);
	}
	// Next check each neighbour of the id
	// There should only be one neighbour that is in the set of reads
	// If this read is our end, then return this arm, otherwise recurse to it
	std::string id_to_take = "";
	for(int i = 0; i < read_indegree[id].size(); i++){
        if(reads.count(read_indegree[id][i]) > 0){
            if(std::find(tmp.begin(), tmp.end(), read_indegree[id][i]) == tmp.end()){
				// Id is in our read set and we haven't taken it
				id_to_take = read_indegree[id][i];
				break;
			}           
        }
    }
    if(id_to_take == ""){
	    for(int i = 0; i < read_outdegree[id].size(); i++){
	        if(reads.count(read_outdegree[id][i]) > 0){
	            if(std::find(tmp.begin(), tmp.end(), read_outdegree[id][i]) == tmp.end()){
					// Id is in our read set and we haven't taken it
					id_to_take = read_outdegree[id][i];
				}
	        }
	    }
	}
	if(id_to_take == end){
		tmp.push_back(id_to_take);
		return;
	}
	if(id_to_take != ""){
		recurse_bubble_arm(id_to_take, end, reads, read_indegree, read_outdegree, tmp);
	}

}

void MatchUtils::get_bubble_arms(std::string start, std::string end, std::set<std::string> reads, std::map<std::string,std::vector<std::string> >& read_indegree, std::map<std::string,std::vector<std::string> >& read_outdegree, std::vector<std::vector<std::string> >& arms){
    // Know that if we get here there is a path between start and end that is seperate from one another.
    // Check indegree and then outdgree. If we have an edge from start to a read in our read set, then recurse over it, until we get to end
    for(int i = 0; i < read_indegree[start].size(); i++){
        if(reads.count(read_indegree[start][i]) > 0){
            // Edge is valid, take it
            std::vector<std::string> tmp;
            tmp.push_back(start);
            recurse_bubble_arm(read_indegree[start][i], end, reads, read_indegree, read_outdegree, tmp);
            arms.push_back(tmp);
        }
    }
    for(int i = 0; i < read_outdegree[start].size(); i++){
        if(reads.count(read_outdegree[start][i]) > 0){
            // Edge is valid, take it
            std::vector<std::string> tmp;
            tmp.push_back(start);
            recurse_bubble_arm(read_outdegree[start][i], end, reads, read_indegree, read_outdegree, tmp);
            arms.push_back(tmp);
        }
    }

}

bool MatchUtils::check_bubble(std::string start, std::string end, std::set<std::string> reads, std::map<std::string,std::vector<std::string> >& read_indegree, std::map<std::string,std::vector<std::string> >& read_outdegree)
{
    // other than start and end reads any reads in set should only have in/out edges to other reads in set
    for (std::set<std::string>::iterator it=reads.begin(); it!=reads.end(); ++it)
    {
        if(start == *it || end == *it){
            // These are the exception nodes, they should be connected to other things in our graph
            continue;
        }
        for (int i = 0; i < read_indegree[*it].size(); ++i)
        {
            if(reads.count(read_indegree[*it][i]) == 0){
                // read at *it has an edge to something that isn't in our read set
                // Thus not a bubble
                return false;
            }
        }
        for (int i = 0; i < read_outdegree[*it].size(); ++i)
        {
            if(reads.count(read_outdegree[*it][i]) == 0){
                // read at *it has an edge to something that isn't in our read set
                // Thus not a bubble
                return false;
            }
        }
    }
    return true;
}


void MatchUtils::compute_sets(std::map<std::string, std::vector<Match> >& all_matches, std::string start, std::string current, std::string end, std::map<std::pair<std::string,std::string>, std::set<std::string> >& bubble_sets, std::set<std::string> seen){
    if(current == end){
        // Reached the node we care about, add the nodes we've seen on the path to here
        bubble_sets[make_pair(start,end)].insert(seen.begin(), seen.end());
        return;
    }
    if(seen.count(current) > 0){
        // Cycle detected
        return;
    }
    seen.insert(current);
    for (int i = 0; i < all_matches[current].size(); ++i)
    {
        if(!all_matches[current][i].reduce){
            MatchUtils::compute_sets(all_matches, start, all_matches[current][i].target_read_id, end, bubble_sets, seen);
        }
    }
}

void MatchUtils::compute_in_out_degree(std::map<std::string, std::vector<Match> >& all_matches, std::set<std::string>& read_ids, std::map<std::string, std::vector<std::string> >& read_indegree, std::map<std::string, std::vector<std::string> >& read_outdegree)
{
    for (std::set<std::string>::iterator it=read_ids.begin(); it!=read_ids.end(); ++it)
    {
    	std::vector<std::string> tmp;
        read_indegree.insert(std::pair<std::string, std::vector<std::string> >(*it,tmp));
        read_outdegree.insert(std::pair<std::string, std::vector<std::string> >(*it,tmp));
    }
    for (std::map<std::string, std::vector<Match> >::iterator it=all_matches.begin(); it!=all_matches.end(); ++it)
    {
        for (int i = 0; i < it->second.size(); ++i)
        {
            if(!it->second[i].reduce){
                if(it->second[i].strand == '+'){
					if(it->second[i].query_read_start > it->second[i].target_read_start){
						read_indegree[it->second[i].target_read_id].push_back(it->second[i].query_read_id);
                		read_outdegree[it->second[i].query_read_id].push_back(it->second[i].target_read_id);
					} else {
						read_outdegree[it->second[i].target_read_id].push_back(it->second[i].query_read_id);
                		read_indegree[it->second[i].query_read_id].push_back(it->second[i].target_read_id);
					}
				} else {
					if(it->second[i].query_read_start > it->second[i].query_read_length - it->second[i].query_read_end && it->second[i].target_read_start > it->second[i].target_read_length - it->second[i].target_read_end){
						read_outdegree[it->second[i].target_read_id].push_back(it->second[i].query_read_id);
                		read_outdegree[it->second[i].query_read_id].push_back(it->second[i].target_read_id);
					} else {
						read_indegree[it->second[i].target_read_id].push_back(it->second[i].query_read_id);
                		read_indegree[it->second[i].query_read_id].push_back(it->second[i].target_read_id);
					}
				}
            }
        }
    }
}

void get_path_to_branch(std::map<std::string, std::vector<Match> >& all_matches, std::vector<std::string>& path, std::string id, std::string prev, std::map<std::string, std::vector<std::string> >& read_indegree, std::map<std::string, std::vector<std::string> >& read_outdegree){
    // First check if the node has any other incoming valid targets, If so then this is branch
    path.push_back(id);
    // need to check if prev is in id's in set or outset
    bool check_outdegree = false;
    if(std::find(read_indegree[id].begin(), read_indegree[id].end(), prev) != read_indegree[id].end()){
    	check_outdegree = true;
    }

    if(check_outdegree){
    	// we came in on indegree, so need to look to see if this is a branch on indegree
    	if(read_indegree[id].size() > 1){
    		// We've reached a branch here
    		return;
    	}

	    if(read_outdegree[id].size() == 1){
	        // If query is a dead end, then it can only have targets, if only one target we can repeat process on target
	        get_path_to_branch(all_matches, path, read_outdegree[id][0] ,id, read_indegree, read_outdegree);
	        return;
	    }
	    std::vector<std::string> tmp;
	    for (int i = 0; i < read_outdegree[id].size(); ++i)
	    {
            if(tmp.size() == 0){
                get_path_to_branch(all_matches, tmp, read_outdegree[id][i], id, read_indegree, read_outdegree);
            } else {
                std::vector<std::string> tmp2;
                get_path_to_branch(all_matches, tmp2, read_outdegree[id][i], id, read_indegree, read_outdegree);
                if(tmp2.size() < tmp.size()){
                    tmp = tmp2;
                }
            }
	    }
	    path.insert(path.end(),tmp.begin(),tmp.end());
	} else {
		if(read_outdegree[id].size() > 1){
    		// We've reached a branch here
    		return;
    	}

	    if(read_indegree[id].size() == 1){
	        // If query is a dead end, then it can only have targets, if only one target we can repeat process on target
	        get_path_to_branch(all_matches, path, read_indegree[id][0] ,id, read_indegree, read_outdegree);
	        return;
	    }
	    std::vector<std::string> tmp;
	    for (int i = 0; i < read_indegree[id].size(); ++i)
	    {
            if(tmp.size() == 0){
                get_path_to_branch(all_matches, tmp, read_indegree[id][i], id, read_indegree, read_outdegree);
            } else {
                std::vector<std::string> tmp2;
                get_path_to_branch(all_matches, tmp2, read_indegree[id][i], id, read_indegree, read_outdegree);
                if(tmp2.size() < tmp.size()){
                    tmp = tmp2;
                }
            }
	    }
	    path.insert(path.end(),tmp.begin(),tmp.end());
	}
}

void MatchUtils::compute_dead_ends(std::map<std::string, std::vector<Match> >& all_matches, std::set<std::string>& read_ids, std::map<std::string, std::vector<std::string> >& read_indegree, std::map<std::string, std::vector<std::string> >& read_outdegree, std::set<std::string>& de_ids, std::map<std::string, std::vector<std::string> >& de_paths)
{
    for (std::set<std::string>::iterator it=read_ids.begin(); it!=read_ids.end(); ++it)
    {
        if(read_indegree[*it].size() == 0 && read_outdegree[*it].size() == 0)
        {
            // Read isn't connected to anything, so It cant be a dead end
            continue;
        } else if(read_indegree[*it].size() == 0 || read_outdegree[*it].size() == 0) 
        {
            //Read must be a dead end, either on in our out direction
            de_ids.insert(*it);
        }
    }

    for (std::set<std::string>::iterator it=de_ids.begin(); it!=de_ids.end(); ++it)
    {
        // Compute the path from each dead end to the branch it came from.
        // If Dead end node starts path then it is query node, but if it ends path then it is the target node.
        if(read_indegree[*it].size() == 0)
        {
            // Read is a valid query only, meaning path starts here and goes to a node that could branch
            std::vector<std::string> tmp;
            tmp.push_back(*it);
            // Compute number of valid targets
            if(read_outdegree[*it].size() == 1){
                // If query is a dead end, then it can only have targets, if only one target we can repeat process on target
                get_path_to_branch(all_matches, tmp, read_outdegree[*it][0], *it, read_indegree, read_outdegree);
            }
            // If query has more than one valid target it is the branch node itself, so no need to find path
            de_paths.insert(std::pair<std::string, std::vector<std::string> >(*it, tmp));
        }
        if(read_outdegree[*it].size() == 0)
        {
            std::vector<std::string> tmp;
            tmp.push_back(*it);
            // Need to find all predacessors of node, if more than one node is the branch, else we have only one and we back track and repeat
            if(read_indegree[*it].size() == 1){
                get_path_to_branch(all_matches, tmp, read_indegree[*it][0], *it, read_indegree, read_outdegree);
            }
            // If query has more than one valid target it is the branch node itself, so no need to find path
            de_paths.insert(std::pair<std::string, std::vector<std::string> >(*it, tmp));
        }
    }
}


void MatchUtils::prune_dead_paths(std::map<std::string, std::vector<Match> >& all_matches, std::set<std::string>& read_ids, std::map<std::string, std::vector<std::string> >& read_indegree, std::map<std::string, std::vector<std::string> >& read_outdegree, std::map<std::string, std::vector<std::string> >& de_paths, int mean_read_length, int threshold)
{
    std::map<std::string, std::vector<int> > branch_count;
    int count = 0;
    for (std::map<std::string, std::vector<std::string> >::iterator it=de_paths.begin(); it!=de_paths.end(); ++it)
    {
    	if(branch_count.count(it->second[it->second.size() - 1]) == 0){
            std::vector<int> tmp;
    		branch_count.insert(make_pair(it->second[it->second.size() - 1], tmp));
    	}
    	branch_count[it->second[it->second.size() - 1]].push_back(it->second.size());
    }
    // Need to check two things, First if a branch has more than one dead end coming back to it, remove the shorter paths, only keep the longest one
    // If there is only one path back to the branch, see how long the path is, if it is short(ie caused by seqencing error ) If sum of suffixes in path is less than 2x mean read length, remove it
    for (std::map<std::string, std::vector<std::string> >::iterator it=de_paths.begin(); it!=de_paths.end(); ++it)
    {
        if(it->second.size() == 1){
            // Branch is a dead end, look at all the neighbors is has, if any of them are dead ends too, remove the edge
            //std::cout << it->first << std::endl;
            for (int i = 0; i < read_indegree[it->first].size(); ++i)
            {
                // for each branch check to see if its neighbors are also dead ends 
                if(read_indegree[read_indegree[it->first][i]].size() == 0 || read_outdegree[read_indegree[it->first][i]].size() == 0){
                    if(it->first < read_indegree[it->first][i]){
                        for (int j = 0; j < all_matches[it->first].size(); ++j)
                        {
                            if(all_matches[it->first][j].target_read_id == read_indegree[it->first][i] || all_matches[it->first][j].query_read_id == read_indegree[it->first][i]){
                                all_matches[it->first][j].reduce = true;
                                count++;
                            }
                        }
                    } else {
                        for (int j = 0; j < all_matches[read_indegree[it->first][i]].size(); ++j)
                        {
                            if(all_matches[read_indegree[it->first][i]][j].target_read_id == it->first || all_matches[read_indegree[it->first][i]][j].query_read_id == it->first){
                                all_matches[read_indegree[it->first][i]][j].reduce = true;
                                count++;
                            }
                        }
                    }
                }
            }
            for (int i = 0; i < read_outdegree[it->first].size(); ++i)
            {
                if(read_outdegree[read_outdegree[it->first][i]].size() == 0 || read_indegree[read_outdegree[it->first][i]].size() == 0){
                    
                    if(it->first < read_outdegree[it->first][i]){
                        for (int j = 0; j < all_matches[it->first].size(); ++j)
                        {
                            if(all_matches[it->first][j].target_read_id == read_outdegree[it->first][i] || all_matches[it->first][j].query_read_id == read_outdegree[it->first][i]){
                                all_matches[it->first][j].reduce = true;
                                count++;
                            }
                        }
                    } else {
                        for (int j = 0; j < all_matches[read_outdegree[it->first][i]].size(); ++j)
                        {
                            if(all_matches[read_outdegree[it->first][i]][j].target_read_id == it->first || all_matches[read_outdegree[it->first][i]][j].query_read_id == it->first){
                                all_matches[read_outdegree[it->first][i]][j].reduce = true;
                                count++;
                            }
                        }
                    }

                }
            }
            continue;
        } else if(branch_count[it->second[it->second.size() - 1]].size() > 1){
            // Have more than one dead end at this branch, check if it is the largest
            if(it->second.size() < *std::max_element( branch_count[it->second[it->second.size() - 1]].begin(), branch_count[it->second[it->second.size() - 1]].end() )){
                // Not the largest path to a dead end out of this branch, so remove it
                for (int i = 0; i < it->second.size()-1; ++i)
                {
                    if(it->second[i] < it->second[i+1]){
                        for (int j = 0; j < all_matches[it->second[i]].size(); ++j)
                        {
                            if(all_matches[it->second[i]][j].target_read_id == it->second[i+1] || all_matches[it->second[i]][j].query_read_id == it->second[i+1]){
                                all_matches[it->second[i]][j].reduce = true;
                                count++;
                            }
                        }
                    } else {
                        for (int j = 0; j < all_matches[it->second[i+1]].size(); ++j)
                        {
                            if(all_matches[it->second[i+1]][j].target_read_id == it->second[i] || all_matches[it->second[i+1]][j].query_read_id == it->second[i]){
                                all_matches[it->second[i+1]][j].reduce = true;
                                count++;
                            }
                        }
                    }
                }
                continue;
            }
        }
        // check the length of the branch
        int size = 0;
        bool larger = false;
        for (int i = 0; i < it->second.size()-1; ++i)
        {   
            if(size > threshold*mean_read_length){
                larger = true;
                break;
            }
            if(it->second[i] < it->second[i+1]){
                for (int j = 0; j < all_matches[it->second[i]].size(); ++j)
                {
                    if(all_matches[it->second[i]][j].target_read_id == it->second[i+1]){
                        size += all_matches[it->second[i]][j].suffix_length;
                    } else if(all_matches[it->second[i]][j].query_read_id == it->second[i+1]){
                        size += all_matches[it->second[i]][j].prefix_length;
                    }
                }
            } else {
                for (int j = 0; j < all_matches[it->second[i+1]].size(); ++j)
                {
                    if(all_matches[it->second[i+1]][j].target_read_id == it->second[i]){
                        size += all_matches[it->second[i+1]][j].suffix_length;
                    } else if(all_matches[it->second[i+1]][j].query_read_id == it->second[i]){
                        size += all_matches[it->second[i+1]][j].prefix_length;
                    }
                }
            }
        }
        if(!larger){
            for (int i = 0; i < it->second.size()-1; ++i)
            {
                if(it->second[i] < it->second[i+1]){
                    for (int j = 0; j < all_matches[it->second[i]].size(); ++j)
                    {
                        if(all_matches[it->second[i]][j].target_read_id == it->second[i+1] || all_matches[it->second[i]][j].query_read_id == it->second[i+1]){
                            all_matches[it->second[i]][j].reduce = true;
                            count++;
                        }
                    }
                } else {
                    for (int j = 0; j < all_matches[it->second[i+1]].size(); ++j)
                    {
                        if(all_matches[it->second[i+1]][j].target_read_id == it->second[i] || all_matches[it->second[i+1]][j].query_read_id == it->second[i]){
                            all_matches[it->second[i+1]][j].reduce = true;
                            count++;
                        }
                    }
                }
            }
        }
    }
    std::cerr << "\tRemoved Edges " << count << std::endl;
}


void MatchUtils::clean_matches(std::map<std::string, std::vector<Match> >& all_matches){
    std::map<std::string, std::vector<Match> > new_matches;
    int count = 0;
    for (std::map<std::string, std::vector<Match> >::iterator it=all_matches.begin(); it!=all_matches.end(); ++it){
        for (int i = 0; i < it->second.size(); ++i){
            Match val = it->second[i];
            if(!val.reduce){
	            if(new_matches.count(val.query_read_id) == 0){
	                std::vector<Match> tmp;
	                new_matches.insert(std::pair<std::string, std::vector<Match> >(val.query_read_id, tmp));
	            }
	            if(new_matches.count(val.target_read_id) == 0){
	                std::vector<Match> tmp;
	                new_matches.insert(std::pair<std::string, std::vector<Match> >(val.target_read_id, tmp));
	            }
	            // need to check if val is already in new matches
	        	if(val.query_read_id > val.target_read_id){
	        		new_matches[val.target_read_id].push_back(val);
	        	} else {
	        		new_matches[val.query_read_id].push_back(val);
	        	}
        	}
        }
    }
    all_matches.clear();
    all_matches = new_matches;
}


int MatchUtils::mark_matches_for_node(std::map<std::string, std::vector<Match> >& all_matches, std::string id, std::map<std::string, int>& mark){
	int count = 0;
	for (std::map<std::string, int>::iterator it=mark.begin(); it!=mark.end(); ++it){
		if(it->second != -1){
			// only care about values that have to be reduced
			continue;
		}
		// compare v and w. If v < w then deal with it in the second section
		// if v > w then we iterate through w's edges and find v and reduce it
		if(id > it->first){
			for (int i = 0; i < all_matches[it->first].size(); ++i){
				if(all_matches[it->first][i].query_read_id == id || all_matches[it->first][i].target_read_id == id){
					if(!all_matches[it->first][i].reduce){
						count++;
						all_matches[it->first][i].reduce = true;
					}
					break;
				}
			}
		}
	}
	for (int i = 0; i < all_matches[id].size(); ++i){
		if(all_matches[id][i].query_read_id == id){
			if(mark[all_matches[id][i].target_read_id] == -1){
				if(!all_matches[id][i].reduce){
					count++;
					all_matches[id][i].reduce = true;
				}
			}
		} else if(all_matches[id][i].target_read_id == id){
			if(mark[all_matches[id][i].query_read_id] == -1){
				if(!all_matches[id][i].reduce){
					count++;
					all_matches[id][i].reduce = true;
				}
			}
		}
	}
	return count;
}


bool check_contigs_for_match(std::string& id, std::string& id2, std::map<int, std::vector<Match> >& contig_map){
    for (std::map<int,std::vector<Match> >::iterator it = contig_map.begin(); it != contig_map.end(); ++it)
    {
        for(int i = 0; i < it->second.size(); i++){
            if(id == it->second[i].target_read_id && id2 == it->second[i].query_read_id){
                return true;
            }
            if(id2 == it->second[i].target_read_id && id == it->second[i].query_read_id){
                return true;
            }
        }
    }
    return false;
}

void get_matches_for_contig(std::string& id, std::string& id2, std::map<std::string, std::vector<Match> >& all_matches, std::map<std::string,std::vector<std::string> >& read_indegree, std::map<std::string,std::vector<std::string> >& read_outdegree, std::vector<Match>& tmp){
    // Always going from id to id2
    // First check that we haven't already visted this edge before for the same contig (Circular contigs specifically)
    for(int i = 0; i < tmp.size(); i++){
        if ((id == tmp[i].target_read_id && id2 == tmp[i].query_read_id) || (id2 == tmp[i].target_read_id && id == tmp[i].query_read_id)){
            return;
        }
    }
    // Second locate the edge from id to id2, and add it to the vector
    // Sorted by lexographically smaller value
    if(id < id2){
        for(int i = 0; i < all_matches[id].size(); i++){
            if ((id == all_matches[id][i].target_read_id && id2 == all_matches[id][i].query_read_id) || (id2 == all_matches[id][i].target_read_id && id == all_matches[id][i].query_read_id)){
                tmp.push_back(all_matches[id][i]);
            }
        }
    } else {
        for(int i = 0; i < all_matches[id2].size(); i++){
            if ((id == all_matches[id2][i].target_read_id && id2 == all_matches[id2][i].query_read_id) || (id2 == all_matches[id2][i].target_read_id && id == all_matches[id2][i].query_read_id)){
                tmp.push_back(all_matches[id2][i]);
            }
        }
    }
    // Next see if id2 is a branch (Either indegree or outdegree must be 2)
    if(read_indegree[id2].size() >= 2 || read_outdegree[id2].size() >= 2){
        return;
    }
    // Check to see if this is a dead end, if id2 is at a dead end then we can end the contig
    // Know that id must be on the otherside to get us here
    if(read_indegree[id2].size() == 0 || read_outdegree[id2].size() == 0){
        return;
    }
    // If not we need to see if we got to id2 from indegree or outdegree, and then continue to the other
    // Know that there must be either 1  in either.
    if(read_indegree[id2][0] == id){
        get_matches_for_contig(id2, read_outdegree[id2][0], all_matches, read_indegree, read_outdegree, tmp);
    } else {
        get_matches_for_contig(id2, read_indegree[id2][0], all_matches, read_indegree, read_outdegree, tmp);
    }
}

int MatchUtils::compute_contigs(std::string id, std::map<std::string, std::vector<Match> >& all_matches, std::map<std::string,std::vector<std::string> >& read_indegree, std::map<std::string,std::vector<std::string> >& read_outdegree, std::map<int, std::vector<Match> >& contig_map, int contig_number){
    // Look at indegree and outdegree of id
    // If it is a branch or a dead end, iterate until we get to a branch
    if((read_indegree[id].size() >= 2) || (read_outdegree[id].size() == 0 && read_indegree[id].size() >=1)){
        // Iterate over each branch
        int count_already_matched = 0;
        for(int i = 0; i < read_indegree[id].size(); i++){
            // Check if branch has already been taken from other direction
            if(!check_contigs_for_match(id, read_indegree[id][i], contig_map)){
                // Get vector of overlaps and save as a contig
                std::vector<Match> tmp;
                get_matches_for_contig(id, read_indegree[id][i], all_matches, read_indegree, read_outdegree, tmp);
                contig_map.insert(make_pair(contig_number, tmp));
                contig_number++;
            } else {
                count_already_matched++;
            }
        }
        // Have the case where the outdegree can be > 0.
        if(read_outdegree[id].size() >= 0){
            for(int i = 0; i < read_outdegree[id].size(); i++){
                // Check if branch has already been taken from other direction
                if(!check_contigs_for_match(id, read_outdegree[id][i], contig_map)){
                    // Get vector of overlaps and save as a contig
                    std::vector<Match> tmp;
                    get_matches_for_contig(id, read_outdegree[id][i], all_matches, read_indegree, read_outdegree, tmp);
                    contig_map.insert(make_pair(contig_number, tmp));
                    contig_number++;
                }
            }
        }
        // Here we have the case where a read is a branch and dead end. ie the end of a bubble. But it needs to be its own contig
        if((read_outdegree[id].size() == 0 && read_indegree[id].size() >=2) && count_already_matched == read_indegree[id].size()){
            std::string id2 = read_indegree[id][0];
            std::vector<Match> tmp;
            if(id < id2){
                for(int i = 0; i < all_matches[id].size(); i++){
                    if ((id == all_matches[id][i].target_read_id && id2 == all_matches[id][i].query_read_id) || (id2 == all_matches[id][i].target_read_id && id == all_matches[id][i].query_read_id)){
                        Match tmp_match = all_matches[id][i];
                        // Want to make sure we use the correct length
                        if(id == tmp_match.target_read_id){
                            tmp_match.length_to_use = tmp_match.target_read_length;
                        } else {
                            tmp_match.length_to_use = tmp_match.query_read_length;
                        }
                        tmp.push_back(tmp_match);
                    }
                }
            } else {
                for(int i = 0; i < all_matches[id2].size(); i++){
                    if ((id == all_matches[id2][i].target_read_id && id2 == all_matches[id2][i].query_read_id) || (id2 == all_matches[id2][i].target_read_id && id == all_matches[id2][i].query_read_id)){
                        Match tmp_match = all_matches[id2][i];
                        // Want to make sure we use the correct length
                        if(id == tmp_match.target_read_id){
                            tmp_match.length_to_use = tmp_match.target_read_length;
                        } else {
                            tmp_match.length_to_use = tmp_match.query_read_length;
                        }
                        tmp.push_back(tmp_match);
                    }
                }
            }
            contig_map.insert(make_pair(contig_number, tmp));
            contig_number++;
        }
    }
    if((read_outdegree[id].size() >= 2) || (read_indegree[id].size() == 0 && read_outdegree[id].size() >=1)){
        int count_already_matched = 0;
        for(int i = 0; i < read_outdegree[id].size(); i++){
            // Check if branch has already been taken from other direction
            if(!check_contigs_for_match(id, read_outdegree[id][i], contig_map)){
                // Get vector of overlaps and save as a contig
                std::vector<Match> tmp;
                get_matches_for_contig(id, read_outdegree[id][i], all_matches, read_indegree, read_outdegree, tmp);
                contig_map.insert(make_pair(contig_number, tmp));
                contig_number++;
            } else {
                count_already_matched++;
            }
        }
        // Have the case where the outdegree can be > 0.
        if(read_indegree[id].size() >= 0){
            for(int i = 0; i < read_indegree[id].size(); i++){
                // Check if branch has already been taken from other direction
                if(!check_contigs_for_match(id, read_indegree[id][i], contig_map)){
                    // Get vector of overlaps and save as a contig
                    std::vector<Match> tmp;
                    get_matches_for_contig(id, read_indegree[id][i], all_matches, read_indegree, read_outdegree, tmp);
                    contig_map.insert(make_pair(contig_number, tmp));
                    contig_number++;
                }
            }
        }
        if((read_indegree[id].size() == 0 && read_outdegree[id].size() >=2) && count_already_matched == read_outdegree[id].size()){
            std::string id2 = read_outdegree[id][0];
            std::vector<Match> tmp;
            if(id < id2){
                for(int i = 0; i < all_matches[id].size(); i++){
                    if ((id == all_matches[id][i].target_read_id && id2 == all_matches[id][i].query_read_id) || (id2 == all_matches[id][i].target_read_id && id == all_matches[id][i].query_read_id)){
                        Match tmp_match = all_matches[id][i];
                        // Want to make sure we use the correct length
                        if(id == tmp_match.target_read_id){
                            tmp_match.length_to_use = tmp_match.target_read_length;
                        } else {
                            tmp_match.length_to_use = tmp_match.query_read_length;
                        }
                        tmp.push_back(tmp_match);
                    }
                }
            } else {
                for(int i = 0; i < all_matches[id2].size(); i++){
                    if ((id == all_matches[id2][i].target_read_id && id2 == all_matches[id2][i].query_read_id) || (id2 == all_matches[id2][i].target_read_id && id == all_matches[id2][i].query_read_id)){
                        Match tmp_match = all_matches[id2][i];
                        // Want to make sure we use the correct length
                        if(id == tmp_match.target_read_id){
                            tmp_match.length_to_use = tmp_match.target_read_length;
                        } else {
                            tmp_match.length_to_use = tmp_match.query_read_length;
                        }
                        tmp.push_back(tmp_match);
                    }
                }
            }
            contig_map.insert(make_pair(contig_number, tmp));
            contig_number++;
        }
    }
    return contig_number;
}


void MatchUtils::reduce_edges(std::map<std::string, std::vector<Match> >& all_matches, std::set<std::string>& read_ids, std::map<std::string, std::vector<Match> >& edge_lists, int fuzz)
{
	//int fuzz = 3500;
	int count = 0;
	std::map<std::string, int> mark;
    // Myers Transitive Reduction Alg
    //Mark each node as vacant (0) & reduce = false
    for (std::set<std::string>::iterator it=read_ids.begin(); it!=read_ids.end(); ++it){
        mark.insert(std::pair<std::string,int>(*it, 0));
    }

    for (std::set<std::string>::iterator it=read_ids.begin(); it!=read_ids.end(); ++it){
        std::string v = *it;
        //std::cerr << v << " " << edge_lists[v].size() << std::endl;
        // Get a vector that has all the edges that come in/out of v (bidirected edges)
        std::vector<Match> v_to_w = edge_lists[v];
        int flongest = 0;
        int rlongest = 0;
        for (int i = 0; i < v_to_w.size(); ++i)
        {	
        	// First compute the suffix/edge length and store it in the match, relative to a specific read
        	// Needed so we can sort by increasing distance
        	//v_to_w[i].set_length(v);
        	//v_to_w[i].set_orientation(v);
            // Set all w for each edge v->w as in play
            if(v_to_w[i].reduce){
            	std::cerr << "Found Already reduced edge from: " << v_to_w[i].query_read_id << " to " << v_to_w[i].target_read_id << std::endl;
            	continue;
            }
        	if(v == v_to_w[i].query_read_id){
            	mark[v_to_w[i].target_read_id] = 1;
        	} else {
        		mark[v_to_w[i].query_read_id] = 1;
        	}
        }
        //Match::sort_matches(v_to_w);
        for (int i = 0; i < v_to_w.size(); ++i)
        {
        	if(v_to_w[i].reduce){
        		continue;
        	}
        	if(v_to_w[i].length > flongest && v_to_w[i].orientation > 0){
                flongest = v_to_w[i].length;
            }
            if(v_to_w[i].length > rlongest && v_to_w[i].orientation < 0){
                rlongest = v_to_w[i].length;
            }
        }
        flongest += fuzz;
        rlongest += fuzz;
        for (int i = 0; i < v_to_w.size(); ++i)
        {
        	int longest = 0;
        	if(v_to_w[i].reduce){
        		//std::cerr << v_to_w[i].target_read_id << " " << v_to_w[i].query_read_id	<< " " << v_to_w[i].type << std::endl;
        		continue;
        	}
        	if(v_to_w[i].orientation < 0){
        		longest = rlongest;
        	} else {
        		longest = flongest;
        	}
            // if w in v->w in play
            std::string w;
            if(v_to_w[i].query_read_id == v){
            	w = v_to_w[i].target_read_id;
            } else {
            	w  = v_to_w[i].query_read_id;
            }
            if(mark[w] == 1){
                //std::cerr << "  " << w << std::endl;
                std::vector<Match> w_to_x = edge_lists[w];
        		//MatchUtils::get_all_matches_for_node(all_matches, w, w_to_x);
        		//for (int j = 0; j < w_to_x.size(); j++)
                //{
                	//w_to_x[j].set_length(w);
        			//w_to_x[j].set_orientation(w);
                	//std::cerr << v << " to " << w_to_x[j].query_read_id << " " << w_to_x[j].target_read_id << " " << v_to_w[i].suffix_length + w_to_x[j].suffix_length << " " << longest << std::endl;
                //}
                //Match::sort_matches(w_to_x);
                for (int j = 0; j < w_to_x.size(); ++j)
                {
		            if(w_to_x[j].reduce){ // ||w_to_x[j].orientation != v_to_w[i].orientation){
		            	// v to w is in one direction, but w to x is in the wrong direction, so ignore it
		            	//mark[w_to_x[j].target_read_id] = -1;
		            	//std::cerr << "A0 Oreintation off " << v << " to " << w_to_x[j].target_read_id << std::endl;
		            	//std::cerr << w_to_x[j].target_read_id << " " << w_to_x[j].query_read_id	<< " " << w_to_x[j].type << std::endl;
		            	continue;
		            }
		            if(v_to_w[i].strand == '-'){
		            	// w is now flipped, so any w to x that is + is really - in the eyes of v 
		            	if(w_to_x[j].orientation == v_to_w[i].orientation){
		            		continue;
		            	}
		            } else if(w_to_x[j].orientation != v_to_w[i].orientation){
		            	continue;
		            }
                	if(w == w_to_x[j].query_read_id){
                		//std::cerr << v << " to " << w_to_x[j].query_read_id << " " << w_to_x[j].target_read_id << std::endl;
                		if(v == w_to_x[j].target_read_id){
		            		//Same edge just reversed, Don't need to remove yourself
		                  //  std::cerr << "Reversed " << v << std::endl;
		            		continue;
		            	}
		                //std::cerr << "For " << v << " Check " << w << " to " << all_matches[w][j].target_read_id << std::endl;
		            	if(v_to_w[i].length + w_to_x[j].length <= longest){
		                	if(mark[w_to_x[j].target_read_id] == 1){
		                    	mark[w_to_x[j].target_read_id] = -1;
		                    //    std::cerr << "A1. Marking: " << v << " to " << w_to_x[j].target_read_id << " For reduction" << std::endl;
		                	}
		            	} else {
		            		break;
		            	}
                	} else {
                		//std::cerr << v << " to " << w_to_x[j].target_read_id << " " << w_to_x[j].query_read_id << std::endl;
                		if(v == w_to_x[j].query_read_id){
		            		//Same edge just reversed, Don't need to remove yourself
		                    //std::cerr << "Reversed " << v << std::endl;
		            		continue;
		            	}
		                //std::cerr << "For " << v << " Check " << w << " to " << all_matches[w][j].target_read_id << std::endl;
		            	if(v_to_w[i].length + w_to_x[j].length <= longest){
		                	if(mark[w_to_x[j].query_read_id] == 1){
		                    	mark[w_to_x[j].query_read_id] = -1;
		                       //std::cerr << "A2. Marking: " << v << " to " << w_to_x[j].query_read_id << " For reduction" << std::endl;
		                	}
		            	} else {
		            		break;
		            	}
                	}

                }
            }
        }
        // This second for loop is skipped in miniasm
        /*
        for (int i = 0; i < v_to_w.size(); ++i)
        {
        	//if(v_to_w[i].orientation < 0){
        		//continue;
        	//}
        	if(v_to_w[i].reduce){
        		//std::cerr << v_to_w[i].target_read_id << " " << v_to_w[i].query_read_id	<< v_to_w[i].type << std::endl;
        		continue;
        	}
            std::string w;
            if(v_to_w[i].query_read_id == v){
            	w = v_to_w[i].target_read_id;
            } else {
            	w  = v_to_w[i].query_read_id;
            }
            // compute smallest edge out of w
            int wmin = 10000000;
            std::vector<Match> w_to_x = edge_lists[w];
        	//MatchUtils::get_all_matches_for_node(all_matches, w, w_to_x);
        	//for (int j = 0; j < w_to_x.size(); j++)
            //{
                // First compute the suffix/edge length and store it in the match, relative to a specific read
	        	// Needed so we can sort by increasing distance
	        //	w_to_x[j].set_length(w);
        	//	w_to_x[j].set_orientation(w);
	        //}
	        //Match::sort_matches(w_to_x);
            for (int j = 0; j < w_to_x.size(); j++){
            	if(w_to_x[j].reduce){
            		continue;
            	}
                if(w_to_x[j].suffix_length < wmin && w_to_x[j].orientation == v_to_w[i].orientation) {
                    wmin = w_to_x[j].suffix_length;
                }
            }
            for (int j = 0; j < w_to_x.size() && (w_to_x[j].suffix_length < fuzz || w_to_x[j].suffix_length == wmin); ++j)
            {
            	if(w_to_x[j].reduce){
            		continue;
            	}
            	if(w == w_to_x[j].query_read_id && w_to_x[j].orientation == v_to_w[i].orientation){
                	if(mark[w_to_x[j].target_read_id] == 1){
                    	mark[w_to_x[j].target_read_id] = -1;
                    	//std::cerr << "B1. Marking: " << v << " to " << w_to_x[j].target_read_id << " For reduction" << std::endl;
                	}
                } else if(w_to_x[j].orientation == v_to_w[i].orientation) {
                	if(mark[w_to_x[j].query_read_id] == 1){
                    	mark[w_to_x[j].query_read_id] = -1;
                    	//std::cerr << "B2. Marking: " << v << " to " << w_to_x[j].query_read_id << " For reduction" << std::endl;
                	}
                }
            }
        }
        */
        count += MatchUtils::mark_matches_for_node(all_matches, v, mark);
        for (int i = 0; i < v_to_w.size(); ++i)
        {
        	if(v_to_w[i].reduce){
        		continue;
        	}
        	if(v == v_to_w[i].query_read_id){
		        mark[v_to_w[i].target_read_id] = 0;
		    } else {
		    	mark[v_to_w[i].query_read_id] = 0;
		    }
        }
    }
    std::cerr << "\tReduced " << count << " edges" << std::endl;
}

void MatchUtils::toGfa(std::map<std::string, std::vector<Match> >& all_matches, std::map<std::string, int>& read_lengths, std::string file_name, std::map<std::string, std::vector<std::string> >& read_indegree, std::map<std::string, std::vector<std::string> >& read_outdegree, std::map<std::string, std::string>& read_names, std::map<std::string, std::string>& colours)
{
	std::ofstream gfaOutput;
  	gfaOutput.open(file_name);
  	gfaOutput << "H\tVN:Z:Test\n";
  	for (std::map<std::string, int>::iterator it=read_lengths.begin(); it!=read_lengths.end(); ++it){
        if(read_outdegree[it->first].size() > 0 || read_indegree[it->first].size() > 0){
  		    gfaOutput << "S\t" << read_names[it->first] <<"\t*\tLN:i:" << it->second << "\tCL:z:" << colours[it->first] << "\tC2:z:" << colours[it->first] << "\n";
        }
  	}
  	for (std::map<std::string, std::vector<Match> >::iterator it=all_matches.begin(); it!=all_matches.end(); ++it){
  		for (int i = 0; i < it->second.size(); ++i)
  		{
  			//std::cerr << it->second[i].query_read_id << " To " << it->second[i].target_read_id << " " << it->second[i].reduce << std::endl;
  			if(!it->second[i].reduce){
  				if(it->second[i].strand == '+'){
					if(it->second[i].query_read_start > it->second[i].target_read_start){
						gfaOutput << "L\t" << read_names[it->second[i].query_read_id] << "\t+\t" << read_names[it->second[i].target_read_id] << "\t+\t"<< it->second[i].cigar <<"\n";
					} else {
						gfaOutput << "L\t" << read_names[it->second[i].target_read_id] << "\t+\t" << read_names[it->second[i].query_read_id] << "\t+\t"<< it->second[i].cigar <<"\n";
					}
				} else {
					if(it->second[i].query_read_start > it->second[i].query_read_length - it->second[i].query_read_end && it->second[i].target_read_start > it->second[i].target_read_length - it->second[i].target_read_end){
						gfaOutput << "L\t" << read_names[it->second[i].query_read_id] << "\t+\t" << read_names[it->second[i].target_read_id] << "\t-\t"<< it->second[i].cigar <<"\n";
					} else {
						gfaOutput << "L\t" << read_names[it->second[i].target_read_id] << "\t-\t" << read_names[it->second[i].query_read_id] << "\t+\t"<< it->second[i].cigar <<"\n";
					}
				}
                
  			}
  		}
  	}
  	gfaOutput.close();
}

