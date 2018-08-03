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

#include "MatchUtils.hpp"

int MatchUtils::read_paf_file(std::map<std::string, std::vector<Match> >& edge_lists, std::map<std::string, std::vector<Match> >& all_matches, std::map<std::string, std::vector<Match> >& raw_matches, std::set<std::string>& read_ids, std::map<std::string, int>& read_lengths, std::string file_name, std::set<std::string>& chimeric_reads, std::map<std::string, Read>& read_classification, bool gfa)
{
	using namespace std;
    io::filtering_istream in_filter;
    in_filter.push(io::gzip_decompressor());
    in_filter.push(io::file_source(file_name));

	//ifstream inputFile_filter(file_name);
	string line;
    set<string> to_drop;
    cerr << "Checking for Contained Reads" << endl;
    int count = 0;
	while (getline(in_filter, line, '\n'))
	{
		istringstream lin(line);
    	string c1, c6, meta, cg;
    	char c5;
    	int c2, c3, c4, c7, c8, c9, c10, c11;
    	lin >> c1 >> c2 >> c3 >> c4 >> c5 >> c6 >> c7 >> c8 >> c9 >> c10 >> c11;
        getline(lin, meta);
        Match tmpLine(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,0,0,"");
        read_ids.insert(c1);
        read_ids.insert(c6);
        if(to_drop.count(c1) > 0 || to_drop.count(c6) > 0 || tmpLine.internal_edge()){
        	continue;
        }
        int contained = tmpLine.check_match_contained();
	    if(contained == 1 && to_drop.count(tmpLine.target_read_id) == 0){
	        to_drop.insert(tmpLine.target_read_id);
	    } else if (contained == -1 && to_drop.count(tmpLine.query_read_id) == 0) {
	        to_drop.insert(tmpLine.query_read_id);
	    }
        if(chimeric_reads.count(c1) != 0){
            to_drop.insert(c1);
        }
        if(chimeric_reads.count(c6) != 0){
            to_drop.insert(c6);
        }
    }
	cerr << "Found: " << to_drop.size() << " Contained Reads and "<< read_ids.size() << " total reads" << endl;
	cerr << "Reading in all valid Matches" << endl;
	io::filtering_istream in;
    in.push(io::gzip_decompressor());
    in.push(io::file_source(file_name));
    count = 0;
    vector<int> sizes;
    read_ids.clear();
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
        if(c1 == c6 || to_drop.count(c1) > 0 || to_drop.count(c6) > 0) {
            if(c1 == c6){
                count++;
            }
            continue;
        }
        // Match is not long enough
        if((c4 - c3) < 2000 || (c9 - c8) < 2000 || c10 < 100){
        	//cerr << c4 - c3 << " " << c8 - c7  << endl;
            count++;
        	continue;
        }
        // Check that the alignments are proper, and at the ends of reads, if not edge reduction will fail    	
        // Generate Unflitered Matches for the raw gfa file
        Match tmpLine(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,0,0,cg);
        if(!tmpLine.internal_edge() && tmpLine.check_match_contained() == 0){
        	if(raw_matches.count(c1) == 0){
        		vector<Match> tmp;
        		raw_matches.insert(pair<string, vector<Match> >(c1,tmp));
        	}
        	if(raw_matches.count(c6) == 0){
        		vector<Match> tmp;
        		raw_matches.insert(pair<string, vector<Match> >(c6,tmp));
        	}
	        // Determine which list to store under based on lexographic comparison
	        // Should never have equality here, check done above already 
	        if(c1 < c6) {
	        	raw_matches[c1].push_back(tmpLine);
	        } else {
	        	raw_matches[c6].push_back(tmpLine);
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
                count += 1;
	        }
	        if(read_lengths.count(c6) == 0){
	            read_lengths.insert(pair<string, int>(c6, c7));
                sizes.push_back(c7);
                count += 1;
	        }
    	} else {
            count++;
        }

    }
    cerr << count << " Self Matches and/or Internal Matches removed" << endl;
	all_matches = raw_matches;
    // Drop any reads found to be contained, Drop all entries to them from all_matches
    for (map<string, vector<Match> >::iterator it2=all_matches.begin(); it2!=all_matches.end(); ++it2)
    {
        for (int i = 0; i < it2->second.size(); ++i)
        {
            if(to_drop.count(it2->second[i].target_read_id) != 0){
                it2->second[i].reduce = true;
            }
            if(to_drop.count(it2->second[i].query_read_id) != 0){
                it2->second[i].reduce = true;
            }
        }
    }
    for (map<string, vector<Match> >::iterator it2=edge_lists.begin(); it2!=edge_lists.end(); ++it2)
    {
    	Match::sort_matches(it2->second);
        for (int i = 0; i < it2->second.size(); ++i)
        {
            if(to_drop.count(it2->second[i].target_read_id) != 0){
                it2->second[i].reduce = true;
            }
            if(to_drop.count(it2->second[i].query_read_id) != 0){
                it2->second[i].reduce = true;
            }
        }
    }
    return accumulate( sizes.begin(), sizes.end(), 0)/sizes.size();
}


void MatchUtils::subset_matches(std::map<std::string, std::vector<Match> >& all_matches, std::map<std::string, std::vector<Match> >& edge_lists, std::map<std::string, std::vector<Match> >& species_matches, std::map<std::string, std::vector<Match> >&  species_edge_lists, std::set<std::string> ids_to_use){
    for (map<string, vector<Match> >::iterator it=all_matches.begin(); it!=all_matches.end(); ++it)
    {
        for (int i = 0; i < it2->second.size(); ++i)
        {
            if(ids_to_use.count(it->second[i].query_read_id) > 0 && ids_to_use.count(it->second[i].target_read_id) > 0){
                if(species_matches.count(it->first) == 0){
                    vector<Match> tmp;
                    species_matches.insert(pair<string, vector<Match> >(it->first,tmp));
                }
                species_matches[it->first].push_back(it->second[i]);
            }
        }
    }

    for (map<string, vector<Match> >::iterator it=edge_lists.begin(); it!=edge_lists.end(); ++it)
    {
        for (int i = 0; i < it2->second.size(); ++i)
        {
            if(ids_to_use.count(it->second[i].query_read_id) > 0 && ids_to_use.count(it->second[i].target_read_id) > 0){
                if(species_edge_lists.count(it->first) == 0){
                    vector<Match> tmp;
                    species_edge_lists.insert(pair<string, vector<Match> >(it->first,tmp));
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
    std::cerr << "Removed Edges " << count << std::endl;
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
    std::cerr << "Reduced " << count << " edges" << std::endl;
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

