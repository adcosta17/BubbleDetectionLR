#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include <fstream>
#include <algorithm>
#include <list>
#include <utility>
#include <cmath>
#include <cstddef>  
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>    
namespace io = boost::iostreams;

#include <dirent.h>

#include "MatchUtils.hpp"

int char_to_hex_int[103] = {
    0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,1,
    2,3,4,5,6,7,8,9,0,0,
    0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,10,11,12,
    13,14,15};

std::string hex_int_string_to_char_string[16] = {
    "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "a", "b", "c", "d", "e", "f"};

std::string MatchUtils::get_hex_string(std::string& str, bool use_names_as_is){
    // Convert ReadId to ACSII represented HexString
    if(use_names_as_is){
        return str;
    }
    std::string to_ret = "";
    int count = 0;
    int tmp = 0;
    for(char& c : str) {
        if(c == '-'){
            continue;
        }
        if(count == 0){
            count +=1;
            tmp = char_to_hex_int[(int)c];
            continue;
        }
        //std::cout << (tmp<<4) << std::endl;
        //std::cout << char_to_hex_int[c] << std::endl;
        //std::cout << (tmp<<4 | char_to_hex_int[c]) << std::endl;
        to_ret += (char) (tmp<<4 | char_to_hex_int[(int)c]);
        count = 0;    
    }
    if(count == 1){
        // Have an odd number of characters
        to_ret += (char) tmp<<4;
        std::string to_ret2 = "0";
        to_ret2 += to_ret;
        return to_ret2;
    }
    return to_ret;
}

std::string MatchUtils::get_read_string(std::string& str, bool use_names_as_is){
    if(use_names_as_is){
        return str;
    }
    // Convert ASCII represented HexString back to String
    std::string to_ret = "";

    for(char& c : str) {
        unsigned char tmp_c = static_cast<unsigned char>(c);
        //std::cout << (unsigned int) tmp_c << std::endl;
        unsigned int v1 = ((unsigned int) tmp_c)>>4;
        //std::cout << v1 << std::endl;

        unsigned int v2 = (((unsigned int) tmp_c)<<28)>>28;
        //std::cout << v2 << std::endl;

        to_ret += hex_int_string_to_char_string[v1];
        to_ret += hex_int_string_to_char_string[v2];
    }

    return to_ret;
}

void get_contained_and_chimeric_reads(std::unordered_set<std::string>& to_drop, std::unordered_set<std::string>& chimeric_reads, std::unordered_set<std::string>& read_ids, std::string file_name, bool reads, bool use_names_as_is){
	using namespace std;
    string suffix = ".paf.gz";
    if(file_name.find(suffix) != string::npos){
        io::filtering_istream in_filter;
        in_filter.push(io::gzip_decompressor());
        in_filter.push(io::file_source(file_name));

    	string line;
    	while (getline(in_filter, line, '\n'))
    	{
    		istringstream lin(line);
        	string c1, c6, meta, cg;
        	char c5;
        	int c2, c3, c4, c7, c8, c9, c10, c11;
        	lin >> c1 >> c2 >> c3 >> c4 >> c5 >> c6 >> c7 >> c8 >> c9 >> c10 >> c11;
            c1 = MatchUtils::get_hex_string(c1, use_names_as_is);
            c6 = MatchUtils::get_hex_string(c6, use_names_as_is);
            //getline(lin, meta);
            if(reads){
            	read_ids.insert(c1);
            	read_ids.insert(c6);
        	}
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
    else {
        ifstream in_filter(file_name);
        string line;
        while (getline(in_filter, line, '\n'))
        {
            istringstream lin(line);
            string c1, c6, meta, cg;
            char c5;
            int c2, c3, c4, c7, c8, c9, c10, c11;
            lin >> c1 >> c2 >> c3 >> c4 >> c5 >> c6 >> c7 >> c8 >> c9 >> c10 >> c11;
            c1 = MatchUtils::get_hex_string(c1, use_names_as_is);
            c6 = MatchUtils::get_hex_string(c6, use_names_as_is);
            //getline(lin, meta);
            if(reads){
                read_ids.insert(c1);
                read_ids.insert(c6);
            }
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
}

std::vector<int> get_all_matches_for_file(std::unordered_map<std::string, std::unordered_map<std::string,Match> >& all_matches, std::unordered_map<std::string, std::vector<std::string>>& matches_indexed, std::unordered_map<std::string, std::unordered_map<std::string,Match> >& raw_matches, std::unordered_set<std::string>& read_ids, std::unordered_map<std::string, int>& read_lengths, std::string file_name, std::unordered_map<std::string, Read>& read_classification, std::unordered_set<std::string>& to_drop, std::string name, bool raw, bool use_names_as_is){
	using namespace std;

    string suffix = ".paf.gz";
    vector<int> sizes;
    string line;
    if(file_name.find(suffix) != string::npos){

    	io::filtering_istream in;
        in.push(io::gzip_decompressor());
        in.push(io::file_source(file_name));
    	while (getline(in, line, '\n'))
    	{
            // Each line of input is split into columns
            // Each line coresponds to a called overlap
    		istringstream lin(line);
        	string c1, c6, meta, cg;
        	char c5;
        	int c2, c3, c4, c7, c8, c9, c10, c11;
        	lin >> c1 >> c2 >> c3 >> c4 >> c5 >> c6 >> c7 >> c8 >> c9 >> c10 >> c11;
            c1 = MatchUtils::get_hex_string(c1, use_names_as_is);
            c6 = MatchUtils::get_hex_string(c6, use_names_as_is);
            //getline(lin, meta);
            //size_t idx = meta.find("cg:");
            cg = "";
            //if(idx != string::npos){
            //    cg = meta.substr(idx+5);
            //}
            //cerr << cg << endl;
            // Check for self alignments && contained reads
            Match tmpLine(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,0,0,cg);
            if(raw){
    	        if(raw_matches.count(c1) == 0){
    	            unordered_map<string,Match> tmp;
    	            raw_matches.insert(pair<string, unordered_map<string,Match> >(c1,tmp));
    	        }
    	        if(raw_matches.count(c6) == 0){
    	            unordered_map<string,Match> tmp;
    	            raw_matches.insert(pair<string, unordered_map<string,Match> >(c6,tmp));
    	        }
            	if(c1 < c6) {
                	raw_matches[c1].insert(pair<string, Match>(c6, tmpLine));
            	} else {
                	raw_matches[c6].insert(pair<string, Match>(c1, tmpLine));
            	}
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
            		unordered_map<string,Match> tmp;
            		all_matches.insert(pair<string, unordered_map<string,Match> >(c1,tmp));
            	}
            	if(all_matches.count(c6) == 0){
            		unordered_map<string,Match> tmp;
            		all_matches.insert(pair<string, unordered_map<string,Match> >(c6,tmp));
            	}
    	        // Determine which list to store under based on lexographic comparison
    	        // Should never have equality here, check done above already
                tmpLine.overlap_species = name;
    	        if(c1 < c6) {
    	        	all_matches[c1].insert(pair<string, Match>(c6, tmpLine));
    	        	if(matches_indexed.count(c6) == 0){
    	        		vector<string> tmp;
    	        		matches_indexed.insert(pair<string, vector<string> >(c6, tmp));
    	        	}
    	        	matches_indexed[c6].push_back(c1);
    	        } else {
    	        	all_matches[c6].insert(pair<string, Match>(c1, tmpLine));
    	        	if(matches_indexed.count(c1) == 0){
    	        		vector<string> tmp;
    	        		matches_indexed.insert(pair<string, vector<string> >(c1, tmp));
    	        	}
    	        	matches_indexed[c1].push_back(c6);
    	        }
    	        read_ids.insert(c1);
    	        read_ids.insert(c6);
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
    }
    else {
        ifstream in(file_name);
        while (getline(in, line, '\n'))
        {
            // Each line of input is split into columns
            // Each line coresponds to a called overlap
            istringstream lin(line);
            string c1, c6, meta, cg;
            char c5;
            int c2, c3, c4, c7, c8, c9, c10, c11;
            lin >> c1 >> c2 >> c3 >> c4 >> c5 >> c6 >> c7 >> c8 >> c9 >> c10 >> c11;
            c1 = MatchUtils::get_hex_string(c1, use_names_as_is);
            c6 = MatchUtils::get_hex_string(c6, use_names_as_is);
            //getline(lin, meta);
            //size_t idx = meta.find("cg:");
            cg = "";
            //if(idx != string::npos){
            //    cg = meta.substr(idx+5);
            //}
            //cerr << cg << endl;
            // Check for self alignments && contained reads
            Match tmpLine(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,0,0,cg);
            if(raw){
                if(raw_matches.count(c1) == 0){
                    unordered_map<string,Match> tmp;
                    raw_matches.insert(pair<string, unordered_map<string,Match> >(c1,tmp));
                }
                if(raw_matches.count(c6) == 0){
                    unordered_map<string,Match> tmp;
                    raw_matches.insert(pair<string, unordered_map<string,Match> >(c6,tmp));
                }
                if(c1 < c6) {
                    raw_matches[c1].insert(pair<string, Match>(c6, tmpLine));
                } else {
                    raw_matches[c6].insert(pair<string, Match>(c1, tmpLine));
                }
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
                    unordered_map<string,Match> tmp;
                    all_matches.insert(pair<string, unordered_map<string,Match> >(c1,tmp));
                }
                if(all_matches.count(c6) == 0){
                    unordered_map<string,Match> tmp;
                    all_matches.insert(pair<string, unordered_map<string,Match> >(c6,tmp));
                }
                // Determine which list to store under based on lexographic comparison
                // Should never have equality here, check done above already
                tmpLine.overlap_species = name;
                if(c1 < c6) {
                    all_matches[c1].insert(pair<string, Match>(c6, tmpLine));
                    if(matches_indexed.count(c6) == 0){
                        vector<string> tmp;
                        matches_indexed.insert(pair<string, vector<string> >(c6, tmp));
                    }
                    matches_indexed[c6].push_back(c1);
                } else {
                    all_matches[c6].insert(pair<string, Match>(c1, tmpLine));
                    if(matches_indexed.count(c1) == 0){
                        vector<string> tmp;
                        matches_indexed.insert(pair<string, vector<string> >(c1, tmp));
                    }
                    matches_indexed[c1].push_back(c6);
                }
                /*
                if(edge_lists.count(c1) == 0){
                    vector<Match> tmp;
                    edge_lists.insert(pair<string, vector<Match> >(c1,tmp));
                }
                if(edge_lists.count(c6) == 0){
                    vector<Match> tmp;
                    edge_lists.insert(pair<string, vector<Match> >(c6,tmp));
                }
                */
                //tmpLine.cigar = "";
                //edge_lists[c1].push_back(tmpLine);
                //Match tmpLine2(c6,c7,c8,c9,c5,c1,c2,c3,c4,c10,c11,0,0,"");
                //tmpLine2.check_match_contained();
                //edge_lists[c6].push_back(tmpLine2);

                read_ids.insert(c1);
                read_ids.insert(c6);
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
    }
    return sizes;
}

int MatchUtils::read_and_assemble_paf_dir_binned(std::unordered_map<std::string, std::unordered_map<std::string,Match> >& all_matches, std::vector<std::string>& n50_values, std::unordered_set<std::string>& read_ids, std::unordered_map<std::string, int>& read_lengths, std::string file_name, std::unordered_set<std::string>& chimeric_reads, std::unordered_map<std::string, Read>& read_classification, int fuzz, int iterations, int threshold, bool use_names_as_is){
    using namespace std;
    unordered_set<string> paf_files;
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

    vector<int> average_read_lens;
    vector<int> total_reads;
    unordered_set<string> to_drop;
    // Take the list of paf_files and then for each of them read in the file
    cerr << "Caculating Contained Reads" << endl;
    for (unordered_set<string>::iterator it = paf_files.begin(); it != paf_files.end(); ++it) {

        // get species name
        string tmp = *it;
        unordered_set<string> tmp_read_ids;
        get_contained_and_chimeric_reads(to_drop, chimeric_reads, tmp_read_ids, tmp, false, use_names_as_is);
        
    }

    for (unordered_set<string>::iterator it = paf_files.begin(); it != paf_files.end(); ++it) {
        string tmp = *it;
        size_t found = tmp.find_last_of("/");
        string name(tmp.substr(found+1).substr(0,tmp.substr(found+1).size()-11));
        cerr << "Assembling " << name << endl;
        unordered_set<string> tmp_read_ids;
        unordered_map<string, unordered_map<string,Match>> tmp_all_matches;
        unordered_map<string, unordered_map<string,Match>> tmp_raw_matches;
        unordered_map<string, vector<string>> matches_indexed;
        unordered_map<string, vector<string> > read_indegree;
        unordered_map<string, vector<string> > read_outdegree;
        cerr << "\tReading in overlaps" << endl;
        vector<int> sizes = get_all_matches_for_file(tmp_all_matches, matches_indexed, tmp_raw_matches, tmp_read_ids, read_lengths, tmp, read_classification, to_drop, name, false, use_names_as_is);
        if(sizes.size() == 0){
            continue;
        }
        int mean_read_length = accumulate( sizes.begin(), sizes.end(), 0)/sizes.size();
        total_reads.push_back(sizes.size());
        average_read_lens.push_back(mean_read_length);

        reduce_edges(tmp_all_matches, matches_indexed, tmp_read_ids, fuzz);
        read_indegree.clear();
        read_outdegree.clear();
        compute_in_out_degree(tmp_all_matches, tmp_read_ids, read_indegree, read_outdegree);
        cerr << "\tDropping Reduced Edges" << endl;
        clean_matches(tmp_all_matches);
        cerr << "\tPruning Dead Ends" << endl;
        int rm_edge = 0;
        for (int j = 0; j < iterations; ++j)
        {
            unordered_set<string> de_ids;
            unordered_map<string, vector<string> > de_paths;
            read_indegree.clear();
            read_outdegree.clear();
            compute_in_out_degree(tmp_all_matches, tmp_read_ids, read_indegree, read_outdegree);
            // Prune Dead End Reads
            compute_dead_ends(tmp_read_ids,read_indegree, read_outdegree, de_ids, de_paths);
            rm_edge += prune_dead_paths(tmp_all_matches, tmp_read_ids, read_indegree, read_outdegree, de_paths, mean_read_length, threshold);
            clean_matches(tmp_all_matches);
        }
        cerr << "\tRemoved " << rm_edge << " Edges" << endl;
        read_indegree.clear();
        read_outdegree.clear();
        compute_in_out_degree(tmp_all_matches, tmp_read_ids, read_indegree, read_outdegree);

        string n50_val = compute_n50(tmp_all_matches, read_indegree, read_outdegree, tmp_read_ids, threshold);
        if(n50_val != ""){
            n50_values.push_back(name + "\n" + n50_val);
        }

        for (unordered_set<string>::iterator it2=tmp_read_ids.begin(); it2!=tmp_read_ids.end(); ++it2){
            if(read_outdegree[*it2].size() > 0 || read_indegree[*it2].size() > 0){
                // Read is in connected componenet, add to set to drop
                read_ids.insert(*it2);
            }
        }

        // Then we also need to take the overlaps for this species and add it to the master set of overalps to look at that are valid
        // This valid set is what we will call bubbles on
        for (unordered_map<string, unordered_map<string,Match> >::iterator it2 = tmp_all_matches.begin(); it2 != tmp_all_matches.end(); ++it2)
        {
            if(all_matches.count(it2->first) == 0){
                unordered_map<string, Match> tmp;
                all_matches.insert(pair<string, unordered_map<string,Match> >(it2->first,tmp));
            }
            for(unordered_map<string, Match>::iterator it3=it2->second.begin(); it3 != it2->second.end(); ++it3){
            	if(all_matches[it2->first].count(it3->first) >0){
            		continue;
            	}
            	all_matches[it2->first].insert(pair<string,Match>(it3->first, it3->second));
            }
        }
    }

    if(average_read_lens.size() > 0){
        int n = accumulate( total_reads.begin(), total_reads.end(), 0);
        int to_ret = 0;
        for(size_t i = 0; i < average_read_lens.size(); i++){
            to_ret += average_read_lens[i] * static_cast<float>(total_reads[i])/n;
        }
        return to_ret;
    }

    return 0;
}



int MatchUtils::read_and_assemble_paf_dir(std::unordered_map<std::string, std::unordered_map<std::string,Match> >& all_matches, std::unordered_map<std::string, std::vector<std::string>>& matches_indexed, std::vector<std::string>& n50_values, std::unordered_set<std::string>& read_ids, std::unordered_map<std::string, int>& read_lengths, std::string file_name, std::unordered_set<std::string>& chimeric_reads, std::unordered_map<std::string, Read>& read_classification, int fuzz, int iterations, int threshold, bool use_names_as_is){
	using namespace std;
	unordered_set<string> paf_files;
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
    vector<int> average_read_lens;
    vector<int> total_reads;
    unordered_set<string> to_drop;
    cerr << "Caculating Contained Reads" << endl;
	for (unordered_set<string>::iterator it = paf_files.begin(); it != paf_files.end(); ++it) {
        unordered_set<string> tmp_read_ids;
        string tmp = *it;
		get_contained_and_chimeric_reads(to_drop, chimeric_reads, tmp_read_ids, tmp, false, use_names_as_is);
        cerr << "Got Contained for " << tmp << endl;
    }
    cerr << "Found " << to_drop.size() << " Contained Reads" << endl;
    unordered_map<string, vector<string> > read_indegree;
    unordered_map<string, vector<string> > read_outdegree;
    int mean_read_length = 0;
    for (unordered_set<string>::iterator it = paf_files.begin(); it != paf_files.end(); ++it) {
        // get species name
        string tmp = *it;
        size_t found = tmp.find_last_of("/");
        string name(tmp.substr(found+1).substr(0,tmp.substr(found+1).size()-11));
        cerr << "Assembling " << name << endl;
        unordered_set<string> tmp_read_ids;
        unordered_map<string, unordered_map<string,Match>> tmp_raw_matches;
        cerr << "\tReading in overlaps" << endl;
		vector<int> sizes = get_all_matches_for_file(all_matches, matches_indexed, tmp_raw_matches, read_ids, read_lengths, tmp, read_classification, to_drop, name, false, use_names_as_is);
        if(sizes.size() == 0){
            continue;
        }
        mean_read_length = accumulate( sizes.begin(), sizes.end(), 0)/sizes.size();
        total_reads.push_back(sizes.size());
        average_read_lens.push_back(mean_read_length);
    }

    if(average_read_lens.size() > 0){
        int n = accumulate( total_reads.begin(), total_reads.end(), 0);
        int to_ret = 0;
        for(size_t i = 0; i < average_read_lens.size(); i++){
            to_ret += average_read_lens[i] * static_cast<float>(total_reads[i])/n;
        }
        mean_read_length = to_ret;
    }

    // Reduce the edges of the combined all_matches set
    reduce_edges(all_matches, matches_indexed, read_ids, fuzz);
    read_indegree.clear();
    read_outdegree.clear();
    compute_in_out_degree(all_matches, read_ids, read_indegree, read_outdegree);
    cerr << "\tDropping Reduced Edges" << endl;
    MatchUtils::clean_matches(all_matches);
    cerr << "\tPruning Dead Ends" << endl;
    int rm_edge = 0;
    for (int j = 0; j < iterations; ++j)
    {
        unordered_set<string> de_ids;
        unordered_map<string, vector<string> > de_paths;
        read_indegree.clear();
        read_outdegree.clear();
        MatchUtils::compute_in_out_degree(all_matches, read_ids, read_indegree, read_outdegree);
        // Prune Dead End Reads
        MatchUtils::compute_dead_ends(read_ids, read_indegree, read_outdegree, de_ids, de_paths);
        rm_edge += MatchUtils::prune_dead_paths(all_matches, read_ids, read_indegree, read_outdegree, de_paths, mean_read_length, threshold);
        MatchUtils::clean_matches(all_matches);
    }
    cerr << "\tRemoved " << rm_edge << " Edges" << endl;
    read_indegree.clear();
    read_outdegree.clear();
    MatchUtils::compute_in_out_degree(all_matches, read_ids, read_indegree, read_outdegree);


    return mean_read_length;
}

int MatchUtils::read_paf_file(std::unordered_map<std::string, std::unordered_map<std::string,Match> >& all_matches, std::unordered_map<std::string, std::vector<std::string>>& matches_indexed, std::unordered_map<std::string, std::unordered_map<std::string,Match> >& raw_matches, std::unordered_set<std::string>& read_ids, std::unordered_map<std::string, int>& read_lengths, std::string file_name, std::unordered_set<std::string>& chimeric_reads, std::unordered_map<std::string, Read>& read_classification, bool gfa, bool use_names_as_is)
{
	using namespace std;
	string line;
    unordered_set<string> to_drop;
    cerr << "Checking for Contained Reads" << endl;
    get_contained_and_chimeric_reads(to_drop, chimeric_reads, read_ids, file_name, true, use_names_as_is);
	cerr << "Found: " << to_drop.size() << " Contained Reads and "<< read_ids.size() << " total reads" << endl;
	cerr << "Reading in all valid Matches" << endl;

	read_ids.clear();
	vector<int> sizes = get_all_matches_for_file(all_matches, matches_indexed, raw_matches, read_ids, read_lengths, file_name, read_classification, to_drop, "", false, use_names_as_is);

    return accumulate( sizes.begin(), sizes.end(), 0)/sizes.size();
}

void MatchUtils::collapseBubble(std::vector<std::vector<std::string> >& arms, std::unordered_map<std::string, std::unordered_map<std::string,Match> >& all_matches){
    // Keep the first arm and remove the rest
    if(arms.size() < 2){
    	//std::cout << "Not enough arms" << std::endl;
        return;
    }
    // Want to remove the smaller arms
    // Compute the size of each arm first, and keep the longest one
    int arm_to_keep = 0;
    std::vector<int> avg_len;
    for (std::size_t i = 0; i < arms.size(); ++i)
    {
        int tmp = 0;
        for (std::size_t j = 0; j < arms[i].size()-1; ++j)
        {
            // Need to find the edge between these 2 reads
            Match tmp_match;
            if(arms[i][j] < arms[i][j+1]){
                // Edge will be on j
                tmp_match = all_matches[arms[i][j]][arms[i][j+1]];
            } else {
                // Edge will be on j+1
                tmp_match = all_matches[arms[i][j+1]][arms[i][j]];
            }
            if(j == 0){
                if(tmp_match.target_read_id == arms[i][j]){
                    tmp += tmp_match.target_read_length;
                } else {
                    tmp += tmp_match.query_read_length;
                }
            }
            tmp += tmp_match.length;
        }
        avg_len.push_back(tmp);
    }
    for (std::size_t i = 0; i < avg_len.size(); ++i)
    {
        if(avg_len[i] > avg_len[arm_to_keep]){
            arm_to_keep = i;
        }
    }
    for (std::size_t i = 0; i < arms.size(); ++i)
    {
        if((int) i == arm_to_keep){
            continue;
        }
        //Remove this arm
        for (std::size_t j = 0; j < arms[i].size()-1; ++j)
        {
            // Edges are ordered so we should be able to find the respective edge in all matches and remove it
            std::string first = arms[i][j];
            std::string second = arms[i][j+1];
            if(first < second){
                all_matches[first][second].reduce = true;
            } else {
                all_matches[second][first].reduce = true;
            }
        }
    }
}

float MatchUtils::getArmLengthRatio(std::vector<std::vector<std::string> >& arms, std::unordered_map<std::string, std::unordered_map<std::string,Match> >& all_matches){
    std::vector<float> avg_len;
    for (std::size_t i = 0; i < arms.size(); ++i)
    {
        float tmp = 0.0;
        for (std::size_t j = 0; j < arms[i].size()-1; ++j)
        {
            // Need to find the edge between these 2 reads
            Match tmp_match;
            if(arms[i][j] < arms[i][j+1]){
                // Edge will be on j
                tmp_match = all_matches[arms[i][j]][arms[i][j+1]];
            } else {
                // Edge will be on j+1
                tmp_match = all_matches[arms[i][j+1]][arms[i][j]];
            }
            if(j == 0){
	            if(tmp_match.target_read_id == arms[i][j]){
	                tmp += tmp_match.target_read_length;
	            } else {
	                tmp += tmp_match.query_read_length;
	            }
	        }
	        tmp += tmp_match.length;
        }
        avg_len.push_back(tmp);
    }
    return avg_len[0]/avg_len[1];
}

float MatchUtils::validBubbleCov(std::vector<std::vector<std::string> >& arms, std::unordered_map<std::string, float>& read_coverage, std::unordered_map<std::string, int>& read_lengths){
    // Can see if there is drastic differences in coverage between the two arms, Assuming that short bubbles can be caused by a sequencing error
    // Smaller sequencing errors should have a lower coverage than true variation
    // We now have each arm of the bubble. Now can compute the coverage for each arm, as the average of the coverage for each read in the arm
    std::vector<float> avg_cov;
    for (std::size_t i = 0; i < arms.size(); ++i)
    {
        float tmp = 0.0;
        int arm_size = 0;
        for (std::size_t j = 0; j < arms[i].size(); ++j)
        {
            tmp += (read_coverage[arms[i][j]] * read_lengths[arms[i][j]]);
            arm_size += read_lengths[arms[i][j]];
        }
        avg_cov.push_back(tmp/arm_size);
    }
    // know there should only be two arms
    float ratio = avg_cov[0]/avg_cov[1];
    return ratio;
}

void MatchUtils::remove_internal_bubbles(std::map<std::pair<std::string,std::string>, std::unordered_set<std::string> >& bubble_sets, std::unordered_map<std::string, std::unordered_map<std::string,Match> >& all_matches, std::unordered_map<std::string,std::vector<std::string> >& read_indegree, std::unordered_map<std::string,std::vector<std::string> >& read_outdegree){
	using namespace std;
	set<pair<string, string> > seen_bubbles;
	for (map<pair<string,string>, unordered_set<string> >::iterator it=bubble_sets.begin(); it!=bubble_sets.end(); ++it)
    {    
       // First check if Bubble is fully contained by a larger bubble
       // If so mark it as seen and remove the edges of it appropraitely based on the in/out degree and start and end postions
        for (map<pair<string,string>, unordered_set<string> >::iterator it2=bubble_sets.begin(); it2!=bubble_sets.end(); ++it2)
        { 
            if(it->first == it2->first || std::make_pair((it->first).second, (it->first).first) == it2->first){
                //same bubble
                continue;
            }
            if(seen_bubbles.count(it->first) != 0 && seen_bubbles.count(std::make_pair((it->first).second, (it->first).first)) != 0){
            	continue;
            }
            if(seen_bubbles.count(it2->first) != 0 && seen_bubbles.count(std::make_pair((it2->first).second, (it2->first).first)) != 0){
            	continue;
            }
            // Now check to see which set is larger
            if(it->second.size() < it2->second.size()){
                unordered_set<string> intersect;
                set_intersection(it->second.begin(),it->second.end(),it2->second.begin(),it2->second.end(),inserter(intersect,intersect.begin()));
                if(intersect.size() == it->second.size()){
                    // All nodes are contained in this larger bubble
                    // Should now evaluate each node and its edges to see if they should be kept or removed 
                    for(unordered_set<string>::iterator it3=it->second.begin(); it3 != it->second.end(); it3++){
                        // First check to see if read is at the start or end of the bigger bubble, if so then don't do anything to it
                        if(*it3 == (it2->first).first || *it3 == (it2->first).second){
                            continue;
                        }
                        // If the read is cleanly part of the bubble then leave it
                        // It Could cleanly be part of the larger bubble as well
                        if(read_outdegree[*it3].size() == 1 && read_indegree[*it3].size() == 1){
                            continue;
                        }
                        // Now read must either be a dead end, which we won't touch in this case
                        // OR the read must have either 2 reads in or 2 out. 
                        // Only Modify the edges if the other read in the pair is also in the bubble
                        if(read_indegree[*it3].size() >= 2 && read_outdegree[*it3].size() != 0){
                            for(std::size_t i = 0; i < read_indegree[*it3].size(); i++){
                                string read_target = read_indegree[*it3][i];
                                // Skip read if it is outside the big bubble
                                if(it2->second.count(read_target) == 0 || read_target == (it2->first).first || read_target == (it2->first).second){
                                    continue;
                                }
                                // Dont look at read_target's full indegree, only bubble specific
                                int rt_in = 0;
                                int rt_out = 0;
                                for(std::size_t j = 0; j < read_indegree[read_target].size(); j++){
                                    rt_in += it2->second.count(read_indegree[read_target][j]);
                                }
                                for(std::size_t j = 0; j < read_outdegree[read_target].size(); j++){
                                    rt_out += it2->second.count(read_outdegree[read_target][j]);
                                }
                                // If this read has multiple edges in and one of those edges is our current node
                                if(rt_in > 1 && std::find(read_indegree[read_target].begin(), read_indegree[read_target].end(), *it3) != read_indegree[read_target].end()){
                                    // Remove this edge
                                    MatchUtils::remove_edge(all_matches, read_indegree, read_outdegree, *it3, read_target);
                                }
                                else if(rt_out > 1 && std::find(read_outdegree[read_target].begin(), read_outdegree[read_target].end(), *it3) != read_outdegree[read_target].end()){
                                    // Remove this edge
                                    MatchUtils::remove_edge(all_matches, read_indegree, read_outdegree, *it3, read_target);
                                }
                            }
                        }
                        if(read_indegree[*it3].size() != 0 && read_outdegree[*it3].size() >= 2){
                            for(std::size_t i = 0; i < read_outdegree[*it3].size(); i++){
                                string read_target = read_outdegree[*it3][i];
                                // Skip read if it is outside the big bubble
                                if(it2->second.count(read_target) == 0 || read_target == (it2->first).first || read_target == (it2->first).second){
                                    continue;
                                }
                                // Dont look at read_target's full indegree, only bubble specific
                                int rt_in = 0;
                                int rt_out = 0;
                                for(std::size_t j = 0; j < read_indegree[read_target].size(); j++){
                                    rt_in += it2->second.count(read_indegree[read_target][j]);
                                }
                                for(std::size_t j = 0; j < read_outdegree[read_target].size(); j++){
                                    rt_out += it2->second.count(read_outdegree[read_target][j]);
                                }
                                // If this read has multiple edges in and one of those edges is our current node
                                if(rt_in > 1 && std::find(read_indegree[read_target].begin(), read_indegree[read_target].end(), *it3) != read_indegree[read_target].end()){
                                    // Remove this edge
                                    MatchUtils::remove_edge(all_matches, read_indegree, read_outdegree, *it3, read_target);
                                }
                                else if(rt_out > 1 && std::find(read_outdegree[read_target].begin(), read_outdegree[read_target].end(), *it3) != read_outdegree[read_target].end()){
                                    // Remove this edge
                                    MatchUtils::remove_edge(all_matches, read_indegree, read_outdegree, *it3, read_target);
                                }
                            }
                        }
                    }
                    seen_bubbles.insert(it->first);
                    seen_bubbles.insert(it2->first);
                } 
            }
        }
    }
}

bool MatchUtils::validBubbleTax(std::vector<std::vector<std::string> >& arms, std::unordered_map<std::string, std::string>& read_lowest_taxonomy){
    // Can use the taxonomy information to see if the arms have distinct species or subspecies
    // Don't want to have shared lowest level of taxinomic classification between arms
    // I.E. if each arm has a subspecies classification and share reads that have the same species classification, ok
    // But if they share a subspecies classification, then need to see what percentage of each arm is what subspecies
    // Base it being a bubble on if it has a majority of one subspecies in each arm (> 75%)
    using namespace std;

    bool valid = false;
    vector<unordered_set<string> > arm_classifcation;
    for (std::size_t i = 0; i < arms.size(); ++i)
    {
        unordered_set<string> tmp;
        for (std::size_t j = 0; j < arms[i].size(); ++j)
        {
            tmp.insert(read_lowest_taxonomy[arms[i][j]]);
        }
        arm_classifcation.push_back(tmp);
    }
    // Now compare the species in each arm to see if there are any that are shared. We should have at least one unique classification in each, ie. distinct sets
    for(std::size_t i = 0; i < arm_classifcation.size(); i++){
        for (std::size_t j = 0; j < arm_classifcation.size(); ++j)
        {
            if(i == j){
                continue;
            }
            unordered_set<string> res1;
            unordered_set<string> res2;
            set_difference( arm_classifcation[i].begin(), arm_classifcation[i].end(), arm_classifcation[j].begin(), arm_classifcation[j].end(), inserter(res1, res1.begin()));
            set_difference( arm_classifcation[j].begin(), arm_classifcation[j].end(), arm_classifcation[i].begin(), arm_classifcation[i].end(), inserter(res2, res2.begin()));
            if(res1.size() > 0 && res2.size() > 0){
                valid = true;
            }
        }
    }
    return valid;
}

std::vector<float> MatchUtils::validBubbleTaxCov(std::vector<std::vector<std::string> >& arms, std::unordered_map<std::string, float>& read_coverage, std::unordered_map<std::string, float>& per_species_coverage, std::unordered_map<std::string, std::string>& read_levels, std::unordered_map<std::string, int>& read_lengths, bool binned, std::unordered_map<std::string, std::unordered_map<std::string,Match> >& all_matches){
    // Can use both the coverage info we have and the taxonomy to get per species/subspecies coverage
    // And then see if it matches what coverage the arms have
    bool print = false;
    // We now have each arm of the bubble. Now can compute the coverage for each arm, as the average of the coverage for each read in the arm
    if(arms[0][0] == "4776" || arms[0][0] == "157"){
        //print = true;
    }
    std::vector<float> avg_cov;
    for (std::size_t i = 0; i < arms.size(); ++i)
    {
        float tmp = 0.0;
        int arm_size = 0;
        for (std::size_t j = 0; j < arms[i].size(); ++j)
        {
            if(print){
                std::cerr << "Coverage of " << arms[i][j] << " " << read_lengths[arms[i][j]] << " " << read_coverage[arms[i][j]] * read_lengths[arms[i][j]] << std::endl;
            }
            tmp += (read_coverage[arms[i][j]] * read_lengths[arms[i][j]]);
            arm_size += read_lengths[arms[i][j]];
            //std::cout << arms[i][j] << " " << read_coverage[arms[i][j]] << std::endl;
        }
        if(print){
            std::cerr << "Arm Coverage " << tmp/arm_size << std::endl;
        }
        avg_cov.push_back(tmp/arm_size);
    }
    // Now need to compute the average coverage for each species/supspecies found in the arm
    // First get the species/subspecies in each arm using read_full_taxonomy
    // Want to also get the length of each arm as a contig at this spot at well
    std::vector<float> arm_tax_cov;
    for (std::size_t i = 0; i < arms.size(); ++i)
    {
        float tmp = 0.0;
        int arm_size = 0;
        if(binned){
            // Need to use all_matches here
            for (std::size_t j = 0; j < arms[i].size()-1; ++j)
            {
                // First find the overlap between read j and j + 1
                std::string species = "";
                if(arms[i][j] < arms[i][j+1]){
                	species = all_matches[arms[i][j]][arms[i][j+1]].overlap_species;
                } else {
                    species = all_matches[arms[i][j+1]][arms[i][j]].overlap_species;
                }
                tmp += per_species_coverage[species] * (read_lengths[arms[i][j]]/2 + read_lengths[arms[i][j+1]]/2);
                arm_size += (read_lengths[arms[i][j]]/2 + read_lengths[arms[i][j+1]]/2);
            }     
        } else {
            for (std::size_t j = 0; j < arms[i].size(); ++j)
            {
                if(j == 0 || j == arms[i].size()-1){
                    tmp += per_species_coverage[read_levels[arms[i][j]]] * read_lengths[arms[i][j]]/2;
                    arm_size += read_lengths[arms[i][j]]/2;
                    continue;
                }
                // Otherwise get the edge lengths of each overlap, and then weight each read's coverage. Each read should have a classification. 
                // If it doesn't then use the read coverage
                if(print){
                    std::cerr << "Species Coverage of " << arms[i][j] << " " << read_lengths[arms[i][j]] << " " << per_species_coverage[read_levels[arms[i][j]]] * read_lengths[arms[i][j]] << std::endl;
                }
                tmp += per_species_coverage[read_levels[arms[i][j]]] * read_lengths[arms[i][j]];
                arm_size += read_lengths[arms[i][j]];
            }
        }
        if(print){
            std::cerr << "Arm Species Coverage " << tmp/arm_size << std::endl;
        }
        arm_tax_cov.push_back(tmp/arm_size);
    }
    // Now that we have the coverage for each arm, and the average coverage based on the taxonomic classifications of the reads in the arm, we can compare and see if they match what we should be seeing
    // Assumption is that each arm should have similar coverage as the species it comes from. If it does we call it a bubble
    std::vector<float> ratio_cov;
    for (std::size_t i = 0; i < avg_cov.size(); ++i)
    {
        // Compare the average coverages between the arm itself and the species
        ratio_cov.push_back(avg_cov[i]/arm_tax_cov[i]);
        //std::cout << avg_cov[i] << " " << arm_tax_cov[i] << std::endl;
        if(print){
            std::cerr << i << " " << avg_cov[i]/arm_tax_cov[i] << std::endl;
        }
    }
    return ratio_cov;
}

std::string MatchUtils::compute_n50(std::unordered_map<std::string, std::unordered_map<std::string,Match> >& all_matches, std::unordered_map<std::string,std::vector<std::string> >& read_indegree, std::unordered_map<std::string,std::vector<std::string> >& read_outdegree, std::unordered_set<std::string>& read_ids, int threshold){
    // Need to evaluate the assembly stats (N50, Overall length, and total number of contigs.)
    // Don't want to fully visualize as we need to see where each read comes from
    // Last step of Layout phase in OLC
    // Look for nodes that are branches. Iterate along each branch and add overlaps to list of edges against respective contig
    // Keep going until we get to a node that is a branch, in either direction.

    using namespace std;

    unordered_map<int, vector<Match> > contigs_sets;
    int contig_num = 0;
    for (unordered_set<string>::iterator it=read_ids.begin(); it!=read_ids.end(); ++it)
    {
        // Look for nodes that have at least 2 valid neighbours, or are dead ends
        // For each one we need to iterate over each neighbour and find path to next branch node
        if(read_indegree[*it].size() >= 2 || read_outdegree[*it].size() >= 2 || (read_indegree[*it].size() == 0 && read_outdegree[*it].size() >= 1) || (read_outdegree[*it].size() == 0 && read_indegree[*it].size() >= 1)) {
            contig_num = compute_contigs(*it, all_matches, read_indegree, read_outdegree, contigs_sets, contig_num);
        }
    }
    if(contig_num == 0){
        // No contigs found
        return "";
    }
    vector<int> contig_lengths;
    int total_length = 0;
    // Now compute the length of each contig
    for (unordered_map<int, vector<Match> >::iterator it = contigs_sets.begin(); it != contigs_sets.end(); ++it)
    {
        int len = 0;
        if(it->second.size() == 1){
            // Pick longer of two reads plus added length
            len = std::max(it->second[0].target_read_length, it->second[0].query_read_length) + it->second[0].length;
        } else {
            for (std::size_t i = 0; i < it->second.size(); ++i)
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
        }
        // Want to have a filter on contig lengths. Don't include contigs that are smaller than the threshold
        if(len >= threshold){
            contig_lengths.push_back(len);
            total_length += len;
        }
    }

    // Compute N50 and report number of contigs and total length
    sort(contig_lengths.begin(), contig_lengths.end(), greater<int>());
    int n50 = 0;
    int sum_len = 0;
    for (std::size_t i = 0; i < contig_lengths.size(); ++i)
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

bool add_overlap(std::vector<std::string>& rin, std::vector<std::string>& rout, std::string target1, std::string target2){
    if(find(rin.begin(), rin.end(), target1) != rin.end() &&
        find(rin.begin(), rin.end(), target2) != rin.end()){
            return false;
    }
    if(find(rout.begin(), rout.end(), target1) != rout.end() &&
        find(rout.begin(), rout.end(), target2) != rout.end()){
            return false;
    }
    return true;
}

bool overlap_orientation(std::vector<std::string>& rin, std::vector<std::string>& rout, std::string target1, std::string target2){
    if(find(rin.begin(), rin.end(), target1) != rin.end() &&
        find(rout.begin(), rout.end(), target2) != rout.end()){
            return true;
    }
    return false;
}


void MatchUtils::collapse_contigs(std::unordered_map<std::string, std::unordered_map<std::string,Match> >& all_matches, std::unordered_map<std::string,std::vector<std::string> >& read_indegree, std::unordered_map<std::string,std::vector<std::string> >& read_outdegree, std::unordered_set<std::string>& read_ids, std::unordered_map<std::string, std::string>& colours, std::unordered_map<std::string, float>& read_coverage, std::string outputFileName, int threshold){
    using namespace std;

    unordered_map<int, vector<Match> > contigs_sets;
    int contig_num = 0;
    for (unordered_set<string>::iterator it=read_ids.begin(); it!=read_ids.end(); ++it)
    {
        // Look for nodes that have at least 2 valid neighbours, or are dead ends
        // For each one we need to iterate over each neighbour and find path to next branch node
        if(read_indegree[*it].size() >= 2 || read_outdegree[*it].size() >= 2 || (read_indegree[*it].size() == 0 && read_outdegree[*it].size() >= 1) || (read_outdegree[*it].size() == 0 && read_indegree[*it].size() >= 1)) {
            //cerr << "Checking " << *it << endl;
            contig_num = MatchUtils::compute_contigs(*it, all_matches, read_indegree, read_outdegree, contigs_sets, contig_num);
        }
    }
    if(contig_num == 0){
        // No contigs found
        return;
    }
    cerr << "Found " << contig_num << " contigs" << endl;
    unordered_map<int, int> contig_lengths;
    // Now compute the length of each contig
    for (unordered_map<int, vector<Match> >::iterator it = contigs_sets.begin(); it != contigs_sets.end(); ++it)
    {
        int len = 0;
        if(it->second.size() == 1){
            len = std::max(it->second[0].target_read_length, it->second[0].query_read_length) + it->second[0].length;
        } else {
            for (std::size_t i = 0; i < it->second.size(); ++i)
            {
                if(i == 0){
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
        }
        if(len >= threshold){
            contig_lengths.insert(make_pair(it->first,len));
        }
    }
    cerr << "Computed contig lengths" << endl;
    // Now need to collapse each contig and generate a set of overlaps between contigs
    // First compute each contig start and end read
    // Then compute the overlaps between contigs
    // Finally output a list of contigs, their ids, lengths and overlaps between them to GFA file, with appropriate colours
    unordered_map<int, vector<string> > contig_ends;
    for (unordered_map<int, int>::iterator it = contig_lengths.begin(); it != contig_lengths.end(); ++it)
    {
        vector<string> tmp;
        if(contigs_sets[it->first].size() == 1){
            tmp.push_back(contigs_sets[it->first][0].query_read_id);
            tmp.push_back(contigs_sets[it->first][0].target_read_id);
            tmp.push_back(contigs_sets[it->first][0].target_read_id);
            tmp.push_back(contigs_sets[it->first][0].query_read_id);
            contig_ends.insert(make_pair(it->first, tmp));
            //cerr << it->first << " " << tmp[0] << " " << tmp[1] << " " << tmp[2] << " " << tmp[3] << endl;
            continue;
        }
        if(contigs_sets[it->first][0].query_read_id == contigs_sets[it->first][1].query_read_id || contigs_sets[it->first][0].query_read_id == contigs_sets[it->first][1].target_read_id){
            tmp.push_back(contigs_sets[it->first][0].target_read_id);
            tmp.push_back(contigs_sets[it->first][0].query_read_id);
        } else {
            tmp.push_back(contigs_sets[it->first][0].query_read_id);
            tmp.push_back(contigs_sets[it->first][0].target_read_id);
        }
        if(contigs_sets[it->first][contigs_sets[it->first].size()-1].query_read_id == contigs_sets[it->first][contigs_sets[it->first].size()-2].query_read_id || 
            contigs_sets[it->first][contigs_sets[it->first].size()-1].query_read_id == contigs_sets[it->first][contigs_sets[it->first].size()-2].target_read_id){
            tmp.push_back(contigs_sets[it->first][contigs_sets[it->first].size()-1].target_read_id);
            tmp.push_back(contigs_sets[it->first][contigs_sets[it->first].size()-1].query_read_id);
        } else {
            tmp.push_back(contigs_sets[it->first][contigs_sets[it->first].size()-1].query_read_id);
            tmp.push_back(contigs_sets[it->first][contigs_sets[it->first].size()-1].target_read_id);
        }
        contig_ends.insert(make_pair(it->first, tmp));
        //cerr << it->first << " " << tmp[0] << " " << tmp[1] << " " << tmp[2] << " " << tmp[3] << endl;
    }
    cerr << "Computed contig ends" << endl;
    // Compute the colours of each contig based on the colours of the reads in it. Exclude grey. 
    // Compute a set of colours. If only one colour in set at end then use it
    // Otherwise use grey. If no colours in set use grey as well
    unordered_map<int, string> contig_colours;
    for (unordered_map<int, int>::iterator it = contig_lengths.begin(); it != contig_lengths.end(); ++it)
    {
        if(colours.size() > 0){
            unordered_set<string> cols;
            for(std::size_t i=0; i<contigs_sets[it->first].size(); i++){
                if(colours[contigs_sets[it->first][i].query_read_id] != "#bebebe"){
                    cols.insert(colours[contigs_sets[it->first][i].query_read_id]);
                }
                if(colours[contigs_sets[it->first][i].target_read_id] != "#bebebe"){
                    cols.insert(colours[contigs_sets[it->first][i].target_read_id]);
                }
            }
            if(cols.size() == 1){
                contig_colours.insert(make_pair(it->first, *cols.begin()));
            } else {
                contig_colours.insert(make_pair(it->first, "#bebebe"));
            }
        } else {
            contig_colours.insert(make_pair(it->first, "#bebebe"));
        }
    }
    cerr << "Computed contig colours" << endl;
    // Now pull up actual overlaps between the start and ends of contigs.
    // unordered_map an overlap between two contigs pair<string, string> to a Match
    // Then when outputing orient edge based on Match orientation as done before for reads
    set<pair<int, int>> contig_matches;
    for (unordered_map<int,vector<string> >::iterator it = contig_ends.begin(); it != contig_ends.end(); ++it)
    {
        for (unordered_map<int, vector<string>>::iterator it2 = contig_ends.begin(); it2 != contig_ends.end(); ++it2)
        {
            if(it->first <= it2->first){
                continue;
            }
            bool it_single_overlap = false;
            bool it2_single_overlap = false;
            if(it->second[0] == it->second[3] || it->second[0] == it->second[2] || it->second[1] == it->second[2] || it->second[1] == it->second[3]){
                it_single_overlap = true;
            }
            if(it2->second[0] == it2->second[3] || it2->second[0] == it2->second[2] || it2->second[1] == it2->second[2] || it2->second[1] == it2->second[3]){
                it2_single_overlap = true;
            }
            // Check all 16 combinations of overlaps
            if(it->second[0] == it2->second[0] && add_overlap(read_indegree[it->second[0]], read_outdegree[it->second[0]], it->second[1], it2->second[1])){
                // Need to check for each pair how to orient the edges
                // If it is in and it2 is out, then it -> it2
                if(overlap_orientation(read_indegree[it->second[0]], read_outdegree[it->second[0]], it->second[1], it2->second[1])){
                    contig_matches.insert(make_pair(it->first, it2->first));
                } else {
                    contig_matches.insert(make_pair(it2->first, it->first));
                }
                //cout << "Adding Edge between " << it->first << " " << it2->first << " " << it->second[0] << " " << it->second[1] << " " << it2->second[1] << endl;
            }
            if(it->second[0] == it2->second[1] && add_overlap(read_indegree[it->second[0]], read_outdegree[it->second[0]], it->second[1], it2->second[0])){
                if(overlap_orientation(read_indegree[it->second[0]], read_outdegree[it->second[0]], it->second[1], it2->second[0])){
                    contig_matches.insert(make_pair(it->first, it2->first));
                } else {
                    contig_matches.insert(make_pair(it2->first, it->first));
                }
                //cout << "Adding Edge between " << it->first << " " << it2->first << " " << it->second[0] << " " << it->second[1] << " " << it2->second[0] << endl;
            }
            if(it->second[1] == it2->second[0] && add_overlap(read_indegree[it->second[1]], read_outdegree[it->second[1]], it->second[0], it2->second[1])){
                if(overlap_orientation(read_indegree[it->second[1]], read_outdegree[it->second[1]], it->second[0], it2->second[1])){
                    contig_matches.insert(make_pair(it->first, it2->first));
                } else {
                    contig_matches.insert(make_pair(it2->first, it->first));
                }
                //cout << "Adding Edge between " << it->first << " " << it2->first << " " << it->second[1] << " " << it->second[0] << " " << it2->second[1] << endl;
            }
            if(it->second[1] == it2->second[1] && add_overlap(read_indegree[it->second[1]], read_outdegree[it->second[1]], it->second[0], it2->second[0])){
                if(overlap_orientation(read_indegree[it->second[1]], read_outdegree[it->second[1]], it->second[0], it2->second[0])){
                    contig_matches.insert(make_pair(it->first, it2->first));
                } else {
                    contig_matches.insert(make_pair(it2->first, it->first));
                }
                //cout << "Adding Edge between " << it->first << " " << it2->first << " " << it->second[1] << " " << it->second[0] << " " << it2->second[0] << endl;
            }
            if(!it2_single_overlap){
                if(it->second[0] == it2->second[2] && add_overlap(read_indegree[it->second[0]], read_outdegree[it->second[0]], it->second[1], it2->second[3])){
                    
                    if(overlap_orientation(read_indegree[it->second[0]], read_outdegree[it->second[0]], it->second[1], it2->second[3])){
                        contig_matches.insert(make_pair(it->first, it2->first));
                    } else {
                        contig_matches.insert(make_pair(it2->first, it->first));
                    }

                    //cout << "Adding Edge between " << it->first << " " << it2->first << " " << it->second[0] << " " << it->second[1] << " " << it2->second[3] << endl;
                }
                if(it->second[0] == it2->second[3] && add_overlap(read_indegree[it->second[0]], read_outdegree[it->second[0]], it->second[1], it2->second[2])){
                    
                    if(overlap_orientation(read_indegree[it->second[0]], read_outdegree[it->second[0]], it->second[1], it2->second[2])){
                        contig_matches.insert(make_pair(it->first, it2->first));
                    } else {
                        contig_matches.insert(make_pair(it2->first, it->first));
                    }

                    //cout << "Adding Edge between " << it->first << " " << it2->first << " " << it->second[0] << " " << it->second[1] << " " << it2->second[2] << endl;
                }
                if(it->second[1] == it2->second[2] && add_overlap(read_indegree[it->second[1]], read_outdegree[it->second[1]], it->second[0], it2->second[3])){
                    
                    if(overlap_orientation(read_indegree[it->second[1]], read_outdegree[it->second[1]], it->second[0], it2->second[3])){
                        contig_matches.insert(make_pair(it->first, it2->first));
                    } else {
                        contig_matches.insert(make_pair(it2->first, it->first));
                    }

                    //cout << "Adding Edge between " << it->first << " " << it2->first << " " << it->second[1] << " " << it->second[0] << " " << it2->second[3] << endl;
                }
                if(it->second[1] == it2->second[3] && add_overlap(read_indegree[it->second[1]], read_outdegree[it->second[1]], it->second[0], it2->second[2])){
                    
                    if(overlap_orientation(read_indegree[it->second[1]], read_outdegree[it->second[1]], it->second[0], it2->second[2])){
                        contig_matches.insert(make_pair(it->first, it2->first));
                    } else {
                        contig_matches.insert(make_pair(it2->first, it->first));
                    }

                    //cout << "Adding Edge between " << it->first << " " << it2->first << " " << it->second[1] << " " << it->second[0] << " " << it2->second[2] << endl;
                }
            }
            if(!it_single_overlap){
                if(it->second[2] == it2->second[0] && add_overlap(read_indegree[it->second[2]], read_outdegree[it->second[2]], it->second[3], it2->second[1])){
                    
                    if(overlap_orientation(read_indegree[it->second[2]], read_outdegree[it->second[2]], it->second[3], it2->second[1])){
                        contig_matches.insert(make_pair(it->first, it2->first));
                    } else {
                        contig_matches.insert(make_pair(it2->first, it->first));
                    }

                    //cout << "Adding Edge between " << it->first << " " << it2->first << " " << it->second[2] << " " << it->second[3] << " " << it2->second[1] << endl;
                }
                if(it->second[2] == it2->second[1] && add_overlap(read_indegree[it->second[2]], read_outdegree[it->second[2]], it->second[3], it2->second[0])){
                    
                    if(overlap_orientation(read_indegree[it->second[2]], read_outdegree[it->second[2]], it->second[3], it2->second[0])){
                        contig_matches.insert(make_pair(it->first, it2->first));
                    } else {
                        contig_matches.insert(make_pair(it2->first, it->first));
                    }

                    //cout << "Adding Edge between " << it->first << " " << it2->first << " " << it->second[2] << " " << it->second[3] << " " << it2->second[0] << endl;
                }
                if(it->second[3] == it2->second[0] && add_overlap(read_indegree[it->second[3]], read_outdegree[it->second[3]], it->second[2], it2->second[1])){
                    
                    if(overlap_orientation(read_indegree[it->second[3]], read_outdegree[it->second[3]], it->second[2], it2->second[1])){
                        contig_matches.insert(make_pair(it->first, it2->first));
                    } else {
                        contig_matches.insert(make_pair(it2->first, it->first));
                    }

                    //cout << "Adding Edge between " << it->first << " " << it2->first << " " << it->second[3] << " " << it->second[2] << " " << it2->second[1] << endl;
                }
                if(it->second[3] == it2->second[1] && add_overlap(read_indegree[it->second[3]], read_outdegree[it->second[3]], it->second[2], it2->second[0])){
                    
                    if(overlap_orientation(read_indegree[it->second[3]], read_outdegree[it->second[3]], it->second[2], it2->second[0])){
                        contig_matches.insert(make_pair(it->first, it2->first));
                    } else {
                        contig_matches.insert(make_pair(it2->first, it->first));
                    }

                    //cout << "Adding Edge between " << it->first << " " << it2->first << " " << it->second[3] << " " << it->second[2] << " " << it2->second[0] << endl;
                }
            }
            if(!it_single_overlap && !it2_single_overlap){
                if(it->second[2] == it2->second[2] && add_overlap(read_indegree[it->second[2]], read_outdegree[it->second[2]], it->second[3], it2->second[3])){
                    
                    if(overlap_orientation(read_indegree[it->second[2]], read_outdegree[it->second[2]], it->second[3], it2->second[3])){
                        contig_matches.insert(make_pair(it->first, it2->first));
                    } else {
                        contig_matches.insert(make_pair(it2->first, it->first));
                    }

                    //cout << "Adding Edge between " << it->first << " " << it2->first << " " << it->second[2] << " " << it->second[3] << " " << it2->second[3] << endl;
                }
                if(it->second[2] == it2->second[3] && add_overlap(read_indegree[it->second[2]], read_outdegree[it->second[2]], it->second[3], it2->second[2])){
                    
                    if(overlap_orientation(read_indegree[it->second[2]], read_outdegree[it->second[2]], it->second[3], it2->second[2])){
                        contig_matches.insert(make_pair(it->first, it2->first));
                    } else {
                        contig_matches.insert(make_pair(it2->first, it->first));
                    }

                    //cout << "Adding Edge between " << it->first << " " << it2->first << " " << it->second[2] << " " << it->second[3] << " " << it2->second[2] << endl;
                }
                if(it->second[3] == it2->second[2] && add_overlap(read_indegree[it->second[3]], read_outdegree[it->second[3]], it->second[2], it2->second[3])){
                    
                    if(overlap_orientation(read_indegree[it->second[3]], read_outdegree[it->second[3]], it->second[2], it2->second[3])){
                        contig_matches.insert(make_pair(it->first, it2->first));
                    } else {
                        contig_matches.insert(make_pair(it2->first, it->first));
                    }

                    //cout << "Adding Edge between " << it->first << " " << it2->first << " " << it->second[3] << " " << it->second[2] << " " << it2->second[3] << endl;
                }
                if(it->second[3] == it2->second[3] && add_overlap(read_indegree[it->second[3]], read_outdegree[it->second[3]], it->second[2], it2->second[2])){
                    
                    if(overlap_orientation(read_indegree[it->second[3]], read_outdegree[it->second[3]], it->second[2], it2->second[2])){
                        contig_matches.insert(make_pair(it->first, it2->first));
                    } else {
                        contig_matches.insert(make_pair(it2->first, it->first));
                    }

                    //cout << "Adding Edge between " << it->first << " " << it2->first << " " << it->second[3] << " " << it->second[2] << " " << it2->second[2] << endl;
                }
            }
        }
    }
    cerr << "Computed contig overlaps" << endl;
    contig_ends.clear();

    unordered_map<int, float> contig_cov;
    for (unordered_map<int, int>::iterator it = contig_lengths.begin(); it != contig_lengths.end(); ++it)
    {
        float tmp = 0;
        for (std::size_t i = 0; i < contigs_sets[it->first].size(); ++i)
        {
            tmp += contigs_sets[it->first][i].length * read_coverage[contigs_sets[it->first][i].query_read_id] * 0.5;
            tmp += contigs_sets[it->first][i].length * read_coverage[contigs_sets[it->first][i].target_read_id] * 0.5;
        }
        contig_cov.insert(make_pair(it->first, tmp/it->second));
    }
    // know there should only be two arms
    cerr << "Output as GFA" << endl;
    // Now just have to take contig lengths colours and overlaps and output to GFA
    ofstream gfaOutput;
    gfaOutput.open(outputFileName);
    gfaOutput << "H\tVN:Z:Collapsed\n";
    for (unordered_map<int, int>::iterator it=contig_lengths.begin(); it!=contig_lengths.end(); ++it){
        gfaOutput << "S\t" << it->first <<"\t*\tLN:i:" << it->second << "\tKC:i:"<< static_cast<int>(it->second*contig_cov[it->first]) << "\tCL:z:" << contig_colours[it->first] << "\tC2:z:" << contig_colours[it->first] << "\n";
    }
    for (set<pair<int, int>>::iterator it = contig_matches.begin(); it!=contig_matches.end(); ++it)
    {
        //cout << "Output Edge For " << it->first << " " << it->second << endl;
        gfaOutput << "L\t" << it->first << "\t+\t" << it->second << "\t+\t\n";
    }
    gfaOutput.close();
}

std::string MatchUtils::compute_ng50(std::unordered_map<std::string, std::unordered_map<std::string,Match> >& all_matches, std::unordered_map<std::string,std::vector<std::string> >& read_indegree, std::unordered_map<std::string,std::vector<std::string> >& read_outdegree, std::unordered_set<std::string>& read_ids, int genome_size, int threshold){
    // Need to evaluate the assembly stats (N50, Overall length, and total number of contigs.)
    // Don't want to fully visualize as we need to see where each read comes from
    // Last step of Layout phase in OLC
    // Look for nodes that are branches. Iterate along each branch and add overlaps to list of edges against respective contig
    // Keep going until we get to a node that is a branch, in either direction.

    using namespace std;

    unordered_map<int, vector<Match> > contigs_sets;
    int contig_num = 0;
    for (unordered_set<string>::iterator it=read_ids.begin(); it!=read_ids.end(); ++it)
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
    for (unordered_map<int, vector<Match> >::iterator it = contigs_sets.begin(); it != contigs_sets.end(); ++it)
    {
        int len = 0;
        if(it->second.size() == 1){
            len = std::max(it->second[0].target_read_length, it->second[0].query_read_length) + it->second[0].length;
        }
        else {
            for (std::size_t i = 0; i < it->second.size(); ++i)
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
        }
        if(len >= threshold){
            contig_lengths.push_back(len);
            total_length += len;
        }
    }

    // Compute N50 and report number of contigs and total length
    sort(contig_lengths.begin(), contig_lengths.end(), greater<int>());
    int ng50 = 0;
    int sum_len = 0;
    for (std::size_t i = 0; i < contig_lengths.size(); ++i)
    {
        sum_len += contig_lengths[i];
        //cout << sum_len << " " << contig_lengths[i] << " " << total_length << " " << sum_len/float(total_length) << endl;
        if(sum_len/float(genome_size) > 0.5){
            ng50 = contig_lengths[i];
            break;
        }
    }

    return string("") + to_string(ng50) + "\t" + to_string(total_length) + "\t" + to_string(contig_lengths.size()) + "\t" + to_string(contig_lengths[0]);
}


void MatchUtils::subset_matches(std::unordered_map<std::string, std::vector<Match> >& all_matches, std::unordered_map<std::string, std::vector<Match> >& edge_lists, std::unordered_map<std::string, std::vector<Match> >& species_matches, std::unordered_map<std::string, std::vector<Match> >&  species_edge_lists, std::unordered_set<std::string> ids_to_use){
    for (std::unordered_map<std::string, std::vector<Match> >::iterator it=all_matches.begin(); it!=all_matches.end(); ++it)
    {
        for (std::size_t i = 0; i < it->second.size(); ++i)
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

    for (std::unordered_map<std::string, std::vector<Match> >::iterator it=edge_lists.begin(); it!=edge_lists.end(); ++it)
    {
        for (std::size_t i = 0; i < it->second.size(); ++i)
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

void MatchUtils::remove_edge(std::unordered_map<std::string, std::unordered_map<std::string,Match> >& all_matches, std::unordered_map<std::string,std::vector<std::string> >& read_indegree, std::unordered_map<std::string,std::vector<std::string> >& read_outdegree, std::string start, std::string end){
    // First remove from all_matches
    //std::cout << "collapseBubble between "<< start << " " << end << std::endl;
    if(start < end){
    	all_matches[start][end].reduce = true;
    } else {
    	all_matches[end][start].reduce = true;
    }

    // Remove edge from read_indegree and out_degree
    if(std::find(read_indegree[start].begin(), read_indegree[start].end(), end) != read_indegree[start].end()){
        std::vector<std::string> tmp;
        for(std::size_t i = 0; i < read_indegree[start].size(); i++){
            if(read_indegree[start][i] != end){
                tmp.push_back(read_indegree[start][i]);
            }
        }
        read_indegree[start] = tmp;
    } else if(std::find(read_outdegree[start].begin(), read_outdegree[start].end(), end) != read_outdegree[start].end()){
        std::vector<std::string> tmp;
        for(std::size_t i = 0; i < read_outdegree[start].size(); i++){
            if(read_outdegree[start][i] != end){
                tmp.push_back(read_outdegree[start][i]);
            }
        }
        read_outdegree[start] = tmp;
    }
    if(std::find(read_indegree[end].begin(), read_indegree[end].end(), start) != read_indegree[end].end()){
        std::vector<std::string> tmp;
        for(std::size_t i = 0; i < read_indegree[end].size(); i++){
            if(read_indegree[end][i] != start){
                tmp.push_back(read_indegree[end][i]);
            }
        }
        read_indegree[end] = tmp;
    } else if(std::find(read_outdegree[end].begin(), read_outdegree[end].end(), start) != read_outdegree[end].end()){
        std::vector<std::string> tmp;
        for(std::size_t i = 0; i < read_outdegree[end].size(); i++){
            if(read_outdegree[end][i] != start){
                tmp.push_back(read_outdegree[end][i]);
            }
        }
        read_outdegree[end] = tmp;
    }
}


void MatchUtils::find_bubble(std::string start, std::unordered_map<std::string,std::vector<std::string> >& read_indegree, std::unordered_map<std::string,std::vector<std::string> >& read_outdegree, std::map<std::pair<std::string,std::string>, std::unordered_set<std::string> >& bubble_sets, std::vector<std::string>& start_ids){

    std::unordered_set<std::string> visited;
    std::list<std::pair<std::string,std::string> > q;
    std::map<std::pair<std::string,std::string>, std::unordered_set<std::string>> bubbles;
    std::pair<std::string,std::string> end = std::make_pair("","");
    visited.insert(start);
    bool print = false;
    //if(start == "1916" || start == "1807"){
    //    print = true;
    //}
    for (std::size_t i = 0; i < start_ids.size(); ++i)
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
            for (std::size_t i = 0; i < read_outdegree[s.first].size(); ++i)
            {
                if(visited.count(read_outdegree[s.first][i]) == 0){
                    visited.insert(read_outdegree[s.first][i]);
                    q.push_back(std::make_pair(read_outdegree[s.first][i], s.first));
                } else {
                    // Found a node we have already seen, possible bubble
                    end = std::make_pair(read_outdegree[s.first][i],s.first);
                    if(bubbles.size() < 2){
                        bubbles.insert(std::make_pair(end, visited));
                    } else {
                        q.clear();
                        break;
                    }
                }
            }
        } else {
            for (std::size_t i = 0; i < read_indegree[s.first].size(); ++i)
            {
                if(visited.count(read_indegree[s.first][i]) == 0){
                    visited.insert(read_indegree[s.first][i]);
                    q.push_back(std::make_pair(read_indegree[s.first][i], s.first));
                } else {
                    // Found a node we have already seen, possible bubble
                    end = std::make_pair(read_indegree[s.first][i],s.first);
                    if(bubbles.size() < 2){
                        bubbles.insert(std::make_pair(end, visited));
                    } else {
                        q.clear();
                        break;
                    }
                }
            }
        }
    }

    // Once we have a node that we have already seen, backtrack and compute the set of visited nodes via a BFS from node we've already seen to the start node
    // Know that 2 paths exist, as proven in the first half so we should be able to find them.
    // Can take the two visited sets and do a intersection, giving us the set of nodes that is on the path between the two exactly
    for(std::map<std::pair<std::string,std::string>, std::unordered_set<std::string>>::iterator it = bubbles.begin(); it != bubbles.end(); it++){
        visited = it->second;
        end = it->first;
        if(end.first != ""){
            if(print){
                std::cout << "Possible Bubble between: " << start << " and " << end.first << std::endl;
            }
            q.clear();
            std::unordered_set<std::string> visited_back;
            visited_back.insert(end.first);
            // Need to check on the direction of our end node, see what the predacesor was and if it is in indegreee or outdegree
            if(std::find(read_indegree[end.first].begin(), read_indegree[end.first].end(), end.second) != read_indegree[end.first].end()){
                for (std::size_t i = 0; i < read_indegree[end.first].size(); ++i)
                {   
                    visited_back.insert(read_indegree[end.first][i]);
                    q.push_back(std::make_pair(read_indegree[end.first][i], end.first));
                }
            } else {
                for (std::size_t i = 0; i < read_outdegree[end.first].size(); ++i)
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
                    for (std::size_t i = 0; i < read_outdegree[e.first].size(); ++i)
                    {
                        if(visited_back.count(read_outdegree[e.first][i]) == 0){
                            visited_back.insert(read_outdegree[e.first][i]);
                            q.push_back(std::make_pair(read_outdegree[e.first][i], e.first));
                        } else {
                            // Found a node we have already seen, possible bubble
                            if(visited_back.count(start) > 0){
                                bubble_end = read_outdegree[e.first][i];
                                if(print){
                                    std::cout << "Adding Bubble between: " << bubble_end << " and " << end.first << std::endl;
                                }
                                std::unordered_set<std::string> combined;
                                std::set_intersection(visited.begin(), visited.end(), visited_back.begin(), visited_back.end(), std::inserter(combined, combined.begin()));
                                bubble_sets.insert(std::pair<std::pair<std::string,std::string>, std::unordered_set<std::string>>(make_pair(bubble_end, end.first), combined));
                                q.clear();
                                break;
                            }
                        }
                    }
                } else {
                    for (std::size_t i = 0; i < read_indegree[e.first].size(); ++i)
                    {
                        if(visited_back.count(read_indegree[e.first][i]) == 0){
                            visited_back.insert(read_indegree[e.first][i]);
                            q.push_back(std::make_pair(read_indegree[e.first][i], e.first));
                        } else {
                            // Found a node we have already seen, possible bubble
                            if(visited_back.count(start) > 0){
                                bubble_end = read_indegree[e.first][i];
                                if(print){
                                    std::cout << "Adding Bubble between: " << bubble_end << " and " << end.first << std::endl;
                                }
                                std::unordered_set<std::string> combined;
                                std::set_intersection(visited.begin(), visited.end(), visited_back.begin(), visited_back.end(), std::inserter(combined, combined.begin()));
                                bubble_sets.insert(std::pair<std::pair<std::string,std::string>, std::unordered_set<std::string>>(make_pair(bubble_end, end.first), combined));
                                q.clear();
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
}

int recurse_bubble_arm(std::string id, std::string end, std::unordered_set<std::string> reads, std::unordered_map<std::string,std::vector<std::string> >& read_indegree, std::unordered_map<std::string,std::vector<std::string> >& read_outdegree, std::vector<std::vector<std::string>>& tmp_arms, int path_num, std::string prev){
    // Check the path that we are on to make sure it doesn't exceed the read set size
    bool print = false;
    //if(end == "2484"){
    //    print = true;
    //}
    if(tmp_arms[path_num].size() > reads.size()){
        return path_num;
    }
    if(std::find(tmp_arms[path_num].begin(), tmp_arms[path_num].end(), id) == tmp_arms[path_num].end()){
		tmp_arms[path_num].push_back(id);
        if(print){
            std::cerr << "Adding " << id << " to path" << std::endl;
        }
	} else {
        if(print){
            std::cerr << "Already Seen "<< id << std::endl;
        }
        return path_num;
    }
	// Next check each neighbour of the id
	// There should only be one neighbour that is in the set of reads
	// If this read is our end, then return this arm, otherwise recurse to it
	std::vector<std::string> ids_to_take;
    // First check which direction we came in on, ie look at the last node in tmp_arms[path_num]
    // Want to leave in opposite direction, If we came in via reads indegree, leave on outdgree
    if(std::find(read_indegree[id].begin(), read_indegree[id].end(), prev) == read_indegree[id].end()){
    	for(std::size_t i = 0; i < read_indegree[id].size(); i++){
            if(reads.count(read_indegree[id][i]) > 0){
                if(std::find(tmp_arms[path_num].begin(), tmp_arms[path_num].end(), read_indegree[id][i]) == tmp_arms[path_num].end()){
    				// Id is in our read set and we haven't taken it
                    if(print){
                        std::cerr << "Adding I " << read_indegree[id][i] << " to ids_to_take" << std::endl;
                    }
    				ids_to_take.push_back(read_indegree[id][i]);
    			}           
            }
        }
    } else {
	    for(std::size_t i = 0; i < read_outdegree[id].size(); i++){
	        if(reads.count(read_outdegree[id][i]) > 0){
	            if(std::find(tmp_arms[path_num].begin(), tmp_arms[path_num].end(), read_outdegree[id][i]) == tmp_arms[path_num].end()){
					// Id is in our read set and we haven't taken it
                    if(print){
                        std::cerr << "Adding O " << read_outdegree[id][i] << " to ids_to_take" << std::endl;
                    }
					ids_to_take.push_back(read_outdegree[id][i]);
				}
	        }
	    }
	}
    if(ids_to_take.size() == 1){
        // Either have another read in path or we have reached the end
        if(ids_to_take[0] == end){
            if(print){
                std::cerr << "Found end, returning" << std::endl;
            }
            tmp_arms[path_num].push_back(end);
            return path_num;
        }
        // Continue on same path with next read
        if(print){
            std::cerr << "Continue to " << ids_to_take[0] << std::endl;
        }
        return recurse_bubble_arm(ids_to_take[0], end, reads, read_indegree, read_outdegree, tmp_arms, path_num, id);
	}
	else if(ids_to_take.size() > 1){
        // Now need to go through each possible branch we have
        int tmp_path_num = path_num;
        std::vector<std::string> path_so_far = tmp_arms[path_num];
        if(ids_to_take[0] == end){
            if(print){
                std::cerr << "Found end, returning" << std::endl;
            }
            tmp_arms[path_num].push_back(end);
        } else {
            if(print){
                std::cerr << "Multiple Path, going to " << ids_to_take[0] << std::endl;
            }
            tmp_path_num = recurse_bubble_arm(ids_to_take[0], end, reads, read_indegree, read_outdegree, tmp_arms, tmp_path_num, id);
        }
        for(std::size_t i = 1; i < ids_to_take.size(); i++){
            if((int) tmp_arms.size() == tmp_path_num + 1){
                // Need to create a new vector that is a copy of current path and
                std::vector<std::string> tmp = path_so_far;
                tmp_arms.push_back(tmp);
            }
            if(ids_to_take[i] == end){
                if(print){
                    std::cerr << "Found end, returning" << std::endl;
                }
                tmp_arms[tmp_path_num + 1].push_back(end);
                tmp_path_num++;
            } else {
                if(print){
                    std::cerr << "Multiple Path, going to " << ids_to_take[0] << std::endl;
                }
                tmp_path_num = recurse_bubble_arm(ids_to_take[i], end, reads, read_indegree, read_outdegree, tmp_arms, tmp_path_num + 1, id);
            }
        }
        return tmp_path_num;
	}
    return path_num;

}

void MatchUtils::get_bubble_arms(std::string start, std::string end, std::unordered_set<std::string> reads, std::unordered_map<std::string,std::vector<std::string> >& read_indegree, std::unordered_map<std::string,std::vector<std::string> >& read_outdegree, std::vector<std::vector<std::string> >& arms){
    // Know that if we get here there is a path between start and end that is seperate from one another.
    // Check indegree and then outdgree. If we have an edge from start to a read in our read set, then recurse over it, until we get to end
    bool print = false;
    //if(start == "3219"){
        //print = true;
    //    std::cerr << "Start "<< start << " End "<< end << std::endl;
    //}
    if(read_indegree[start].size() > 1){
        std::vector<std::vector<std::string>> tmp_arms;
        for(std::size_t i = 0; i < read_indegree[start].size(); i++){
            if(reads.count(read_indegree[start][i]) > 0){
                // Edge is valid, take it
                std::vector<std::string> tmp;
                tmp.push_back(start);
                tmp_arms.push_back(tmp);
                recurse_bubble_arm(read_indegree[start][i], end, reads, read_indegree, read_outdegree, tmp_arms, tmp_arms.size() - 1, start);
            }
        }
        std::vector<bool> to_keep;
        int keep_count = 0;
        std::unordered_map<int, int> ms_count;
        if(print){
            std::cerr << "Arms:" << std::endl;
            for(std::size_t i = 0; i < tmp_arms.size(); i++){
                for(std::size_t j = 0; j < tmp_arms[i].size(); j++){
                    std::cerr << tmp_arms[i][j] << " ";
                }
                std::cerr << std::endl;
            }
            std::cerr << std::endl;
        }
        for(std::size_t i = 0; i < tmp_arms.size(); i++){
            // Check if we get to the end in each arm of temp_arms
            ms_count.insert(std::make_pair(i,0));
            if(std::find(tmp_arms[i].begin(), tmp_arms[i].end(), end) == tmp_arms[i].end()){
                to_keep.push_back(false);
            } else {
                to_keep.push_back(true);
                keep_count++;
            }
        }

        if(tmp_arms.size() < 2){
            return;
        }
        if(tmp_arms.size() == 2){
            if(to_keep[0] == true && to_keep[1] == true){
                arms = tmp_arms;
            }
            return;
        }
        if(keep_count < 2){
            return;
        }
        if(keep_count == 2){
            for(std::size_t i = 0; i < tmp_arms.size(); i++){
                if(to_keep[i]){
                    arms.push_back(tmp_arms[i]);
                }
            }
            return;
        }
        // now check for mutually exclusive arms, remove start and end points from each arm and see
        for(std::size_t i = 0; i< tmp_arms.size(); i++){
            for(std::size_t j = 1; j< tmp_arms.size(); j++){
                if(j <= i){
                    // either seen this combination or both are same
                    continue;
                }
                if(!to_keep[i] || !to_keep[j]){
                    continue;
                }
                std::vector<std::string> v1 = tmp_arms[i];
                std::vector<std::string> v2 = tmp_arms[j];
                std::vector<std::string> v3;
                std::sort(v1.begin(), v1.end());
                std::sort(v2.begin(), v2.end());
                std::set_intersection(v1.begin(),v1.end(),v2.begin(),v2.end(),std::back_inserter(v3));
                if(v3.size() == 2){
                    // should have just the start and end reads that are the same
                    ms_count[i]++;
                    ms_count[j]++;
                }
            }
        }
        // Now find all arms that have the maximum count
        int max_val = 0;
        for(std::size_t i = 0; i < tmp_arms.size(); i++){
            if(ms_count[i] > max_val){
                max_val = ms_count[i];
            }
        }
        for(std::size_t i = 0; i < tmp_arms.size(); i++){
            if(ms_count[i] == max_val){
                arms.push_back(tmp_arms[i]);
            }
        }
        return;
    }
    if(read_outdegree[start].size() > 1 && arms.size() == 0){
        // If there are arms then we must have hit the end from paths above, no need to check outdegree for this bubble
        std::vector<std::vector<std::string>> tmp_arms;
        for(std::size_t i = 0; i < read_outdegree[start].size(); i++){
            if(reads.count(read_outdegree[start][i]) > 0){
                // Edge is valid, take it
                std::vector<std::string> tmp;
                tmp.push_back(start);
                tmp_arms.push_back(tmp);
                recurse_bubble_arm(read_outdegree[start][i], end, reads, read_indegree, read_outdegree, tmp_arms, tmp_arms.size() - 1, start);              
            }
        }
        std::vector<bool> to_keep;
        int keep_count = 0;
        std::unordered_map<int, int> ms_count;
        if(print){
            std::cerr << "Arms:" << std::endl;
            for(std::size_t i = 0; i < tmp_arms.size(); i++){
                for(std::size_t j = 0; j < tmp_arms[i].size(); j++){
                    std::cerr << tmp_arms[i][j] << " ";
                }
                std::cerr << std::endl;
            }
            std::cerr << std::endl;
        }
        for(std::size_t i = 0; i < tmp_arms.size(); i++){
            // Check if we get to the end in each arm of temp_arms
            ms_count.insert(std::make_pair(i,0));
            if(std::find(tmp_arms[i].begin(), tmp_arms[i].end(), end) == tmp_arms[i].end()){
                to_keep.push_back(false);
            } else {
                to_keep.push_back(true);
                keep_count++;
            }
        }

        if(tmp_arms.size() < 2){
            return;
        }
        if(tmp_arms.size() == 2){
            if(to_keep[0] == true && to_keep[1] == true){
                arms = tmp_arms;
            }
            return;
        }
        if(keep_count < 2){
            return;
        }
        if(keep_count == 2){
            for(std::size_t i = 0; i < tmp_arms.size(); i++){
                if(to_keep[i]){
                    arms.push_back(tmp_arms[i]);
                }
            }
            return;
        }
        for(std::size_t i = 0; i< tmp_arms.size(); i++){
            for(std::size_t j = 1; j< tmp_arms.size(); j++){
                if(j <= i){
                    // either seen this combination or both are same
                    continue;
                }
                if(!to_keep[i] || !to_keep[j]){
                    continue;
                }
                std::vector<std::string> v3;

                std::sort(tmp_arms[i].begin(), tmp_arms[i].end());
                std::sort(tmp_arms[j].begin(), tmp_arms[j].end());

                std::set_intersection(tmp_arms[i].begin(),tmp_arms[i].end(),tmp_arms[j].begin(),tmp_arms[j].end(),std::back_inserter(v3));
                if(v3.size() == 2){
                    // should have just the start and end reads that are the same
                    ms_count[i]++;
                    ms_count[j]++;
                }
            }
        }
        // Now find all arms that have the maximum count
        int max_val = 0;
        for(std::size_t i = 0; i < tmp_arms.size(); i++){
            if(ms_count[i] > max_val){
                max_val = ms_count[i];
            }
        }
        for(std::size_t i = 0; i < tmp_arms.size(); i++){
            if(ms_count[i] == max_val){
                arms.push_back(tmp_arms[i]);
            }
        }
    }
}

bool MatchUtils::check_bubble(std::string start, std::string end, std::unordered_set<std::string> reads, std::unordered_map<std::string,std::vector<std::string> >& read_indegree, std::unordered_map<std::string,std::vector<std::string> >& read_outdegree)
{
    // other than start and end reads any reads in set should only have in/out edges to other reads in set
    for (std::unordered_set<std::string>::iterator it=reads.begin(); it!=reads.end(); ++it)
    {
        if(start == *it || end == *it){
            // These are the exception nodes, they should be connected to other things in our graph
            continue;
        }
        for (std::size_t i = 0; i < read_indegree[*it].size(); ++i)
        {
            if(reads.count(read_indegree[*it][i]) == 0){
                // read at *it has an edge to something that isn't in our read set
                // Thus not a bubble
                return false;
            }
        }
        for (std::size_t i = 0; i < read_outdegree[*it].size(); ++i)
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


void MatchUtils::compute_sets(std::unordered_map<std::string, std::vector<Match> >& all_matches, std::string start, std::string current, std::string end, std::map<std::pair<std::string,std::string>, std::unordered_set<std::string> >& bubble_sets, std::unordered_set<std::string> seen){
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
    for (std::size_t i = 0; i < all_matches[current].size(); ++i)
    {
        if(!all_matches[current][i].reduce){
            MatchUtils::compute_sets(all_matches, start, all_matches[current][i].target_read_id, end, bubble_sets, seen);
        }
    }
}

void MatchUtils::compute_in_out_degree(std::unordered_map<std::string, std::unordered_map<std::string, Match> >& all_matches, std::unordered_set<std::string>& read_ids, std::unordered_map<std::string, std::vector<std::string> >& read_indegree, std::unordered_map<std::string, std::vector<std::string> >& read_outdegree)
{
    for (std::unordered_set<std::string>::iterator it=read_ids.begin(); it!=read_ids.end(); ++it)
    {
    	std::vector<std::string> tmp;
        read_indegree.insert(std::pair<std::string, std::vector<std::string> >(*it,tmp));
        read_outdegree.insert(std::pair<std::string, std::vector<std::string> >(*it,tmp));
    }
    for (std::unordered_map<std::string, std::unordered_map<std::string, Match> >::iterator it=all_matches.begin(); it!=all_matches.end(); ++it)
    {
        for (std::unordered_map<std::string, Match>::iterator it2=it->second.begin(); it2!=it->second.end(); ++it2)
        {
            if(!it2->second.reduce){
                if(it2->second.strand == '+'){
					if(it2->second.query_read_start > it2->second.target_read_start){
						read_indegree[it2->second.target_read_id].push_back(it2->second.query_read_id);
                		read_outdegree[it2->second.query_read_id].push_back(it2->second.target_read_id);
					} else {
						read_outdegree[it2->second.target_read_id].push_back(it2->second.query_read_id);
                		read_indegree[it2->second.query_read_id].push_back(it2->second.target_read_id);
					}
				} else {
					if(it2->second.query_read_start > it2->second.query_read_length - it2->second.query_read_end && it2->second.target_read_start > it2->second.target_read_length - it2->second.target_read_end){
						read_outdegree[it2->second.target_read_id].push_back(it2->second.query_read_id);
                		read_outdegree[it2->second.query_read_id].push_back(it2->second.target_read_id);
					} else {
						read_indegree[it2->second.target_read_id].push_back(it2->second.query_read_id);
                		read_indegree[it2->second.query_read_id].push_back(it2->second.target_read_id);
					}
				}
            }
        }
    }
}

void get_path_to_branch(std::vector<std::string>& path, std::string id, std::string prev, std::unordered_map<std::string, std::vector<std::string> >& read_indegree, std::unordered_map<std::string, std::vector<std::string> >& read_outdegree){
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
	        get_path_to_branch(path, read_outdegree[id][0] ,id, read_indegree, read_outdegree);
	        return;
	    }
	    std::vector<std::string> tmp;
	    for (std::size_t i = 0; i < read_outdegree[id].size(); ++i)
	    {
            if(tmp.size() == 0){
                get_path_to_branch(tmp, read_outdegree[id][i], id, read_indegree, read_outdegree);
            } else {
                std::vector<std::string> tmp2;
                get_path_to_branch(tmp2, read_outdegree[id][i], id, read_indegree, read_outdegree);
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
	        get_path_to_branch(path, read_indegree[id][0] ,id, read_indegree, read_outdegree);
	        return;
	    }
	    std::vector<std::string> tmp;
	    for (std::size_t i = 0; i < read_indegree[id].size(); ++i)
	    {
            if(tmp.size() == 0){
                get_path_to_branch(tmp, read_indegree[id][i], id, read_indegree, read_outdegree);
            } else {
                std::vector<std::string> tmp2;
                get_path_to_branch(tmp2, read_indegree[id][i], id, read_indegree, read_outdegree);
                if(tmp2.size() < tmp.size()){
                    tmp = tmp2;
                }
            }
	    }
	    path.insert(path.end(),tmp.begin(),tmp.end());
	}
}

void MatchUtils::compute_dead_ends(std::unordered_set<std::string>& read_ids, std::unordered_map<std::string, std::vector<std::string> >& read_indegree, std::unordered_map<std::string, std::vector<std::string> >& read_outdegree, std::unordered_set<std::string>& de_ids, std::unordered_map<std::string, std::vector<std::string> >& de_paths)
{
    for (std::unordered_set<std::string>::iterator it=read_ids.begin(); it!=read_ids.end(); ++it)
    {
        if(read_indegree[*it].size() == 0 && read_outdegree[*it].size() == 0)
        {
            // Read isn't connected to anything, so It cant be a dead end
            //std::cout<< "Size 0 " << *it << std::endl;
            continue;
        } else if(read_indegree[*it].size() == 0 || read_outdegree[*it].size() == 0) 
        {
            //Read must be a dead end, either on in our out direction
            de_ids.insert(*it);
            //std::cout<< "DE " << *it << std::endl;
        }
    }

    for (std::unordered_set<std::string>::iterator it=de_ids.begin(); it!=de_ids.end(); ++it)
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
                get_path_to_branch(tmp, read_outdegree[*it][0], *it, read_indegree, read_outdegree);
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
                get_path_to_branch(tmp, read_indegree[*it][0], *it, read_indegree, read_outdegree);
            }
            // If query has more than one valid target it is the branch node itself, so no need to find path
            de_paths.insert(std::pair<std::string, std::vector<std::string> >(*it, tmp));
        }
    }
}


int MatchUtils::prune_dead_paths(std::unordered_map<std::string, std::unordered_map<std::string,Match> >& all_matches, std::unordered_set<std::string>& read_ids, std::unordered_map<std::string, std::vector<std::string> >& read_indegree, std::unordered_map<std::string, std::vector<std::string> >& read_outdegree, std::unordered_map<std::string, std::vector<std::string> >& de_paths, int mean_read_length, int threshold)
{
    std::unordered_map<std::string, std::vector<int> > branch_count;
    int count = 0;
    for (std::unordered_map<std::string, std::vector<std::string> >::iterator it=de_paths.begin(); it!=de_paths.end(); ++it)
    {
    	if(branch_count.count(it->second[it->second.size() - 1]) == 0){
            std::vector<int> tmp;
    		branch_count.insert(make_pair(it->second[it->second.size() - 1], tmp));
    	}
    	branch_count[it->second[it->second.size() - 1]].push_back(it->second.size());
    }
    // Need to check two things, First if a branch has more than one dead end coming back to it, remove the shorter paths, only keep the longest one
    // If there is only one path back to the branch, see how long the path is, if it is short(ie caused by seqencing error ) If sum of suffixes in path is less than 2x mean read length, remove it
    for (std::unordered_map<std::string, std::vector<std::string> >::iterator it=de_paths.begin(); it!=de_paths.end(); ++it)
    {
        //std::cout << it->first << std::endl;
        if(it->second.size() == 1){
            // Branch is a dead end, look at all the neighbors is has, if any of them are dead ends too, remove the edge
            for (std::size_t i = 0; i < read_indegree[it->first].size(); ++i)
            {
                // for each branch check to see if its neighbors are also dead ends 
                if(read_indegree[read_indegree[it->first][i]].size() == 0 || read_outdegree[read_indegree[it->first][i]].size() == 0){
                    if(it->first < read_indegree[it->first][i]){
                        all_matches[it->first][read_indegree[it->first][i]].reduce = true;
                        count++;
                    } else {
                        all_matches[read_indegree[it->first][i]][it->first].reduce = true;
                        count++;
                    }
                }
            }
            for (std::size_t i = 0; i < read_outdegree[it->first].size(); ++i)
            {
                if(read_outdegree[read_outdegree[it->first][i]].size() == 0 || read_indegree[read_outdegree[it->first][i]].size() == 0){
                    if(it->first < read_outdegree[it->first][i]){
                        all_matches[it->first][read_outdegree[it->first][i]].reduce = true;
                        count++;
                    } else {
                        all_matches[read_outdegree[it->first][i]][it->first].reduce = true;
                        count++;
                    }
                }
            }
            continue;
        } else if(branch_count[it->second[it->second.size() - 1]].size() > 1){
            // Have more than one dead end at this branch, check if it is the largest
            if((int) it->second.size() < *std::max_element( branch_count[it->second[it->second.size() - 1]].begin(), branch_count[it->second[it->second.size() - 1]].end() )){
                // Not the largest path to a dead end out of this branch, so remove it
                for (std::size_t i = 0; i < it->second.size()-1; ++i)
                {
                    if(it->second[i] < it->second[i+1]){
                    	all_matches[it->second[i]][it->second[i+1]].reduce = true;
                    	count++;
                    } else {
                    	all_matches[it->second[i+1]][it->second[i]].reduce = true;
                    	count++;
                    }
                }
                continue;
            }
        }
        // check the length of the branch
        int size = 0;
        bool larger = false;
        for (std::size_t i = 0; i < it->second.size()-1; ++i)
        {   
            if(size > threshold*mean_read_length){
                larger = true;
                break;
            }
            if(it->second[i] < it->second[i+1]){
            	size += all_matches[it->second[i]][it->second[i+1]].length;
            } else {
                size += all_matches[it->second[i+1]][it->second[i]].length;
            }
        }
        if(!larger){
            for (std::size_t i = 0; i < it->second.size()-1; ++i)
            {
                if(it->second[i] < it->second[i+1]){
                    all_matches[it->second[i]][it->second[i+1]].reduce = true;
                    count++;
                } else {
                    all_matches[it->second[i+1]][it->second[i]].reduce = true;
                    count++;
                }
            }
        }
    }
    return count;
}


void MatchUtils::clean_matches(std::unordered_map<std::string, std::unordered_map<std::string, Match> >& all_matches){

    for (std::unordered_map<std::string, std::unordered_map<std::string, Match> >::iterator it=all_matches.begin(); it!=all_matches.end(); ++it){
        std::unordered_set<std::string> to_remove;
        for (std::unordered_map<std::string, Match>::iterator it2=it->second.begin(); it2!=it->second.end(); ++it2){
            if(it2->second.reduce){
                to_remove.insert(it2->first);
            }
        }
        for (std::unordered_set<std::string>::iterator it2=to_remove.begin(); it2!=to_remove.end(); ++it2){
            it->second.erase(*it2);
        }
    }
}


int mark_matches_for_node(std::unordered_map<std::string, std::unordered_map<std::string, Match> >& all_matches, std::string id, std::unordered_map<std::string, int>& mark){
	int count = 0;
	for (std::unordered_map<std::string, int>::iterator it=mark.begin(); it!=mark.end(); ++it){
		if(it->second != -1){
			// only care about values that have to be reduced
			continue;
		}
		// compare v and w. If v < w then deal with it in the second section
		// if v > w then we iterate through w's edges and find v and reduce it
		if(id < it->first){
            all_matches[id][it->first].reduce = true;
            count ++;
		} else {
            all_matches[it->first][id].reduce = true;
            count ++;
        }
	}
	return count;
}


bool check_contigs_for_match(std::string& id, std::string& id2, std::unordered_map<int, std::vector<Match> >& contig_unordered_map){
    for (std::unordered_map<int,std::vector<Match> >::iterator it = contig_unordered_map.begin(); it != contig_unordered_map.end(); ++it)
    {
        for(std::size_t i = 0; i < it->second.size(); i++){
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

void get_matches_for_contig(std::string& id, std::string& id2, std::unordered_map<std::string, std::unordered_map<std::string,Match> >& all_matches, std::unordered_map<std::string,std::vector<std::string> >& read_indegree, std::unordered_map<std::string,std::vector<std::string> >& read_outdegree, std::vector<Match>& tmp){
    // Always going from id to id2
    // First check that we haven't already visted this edge before for the same contig (Circular contigs specifically)
    for(std::size_t i = 0; i < tmp.size(); i++){
        if ((id == tmp[i].target_read_id && id2 == tmp[i].query_read_id) || (id2 == tmp[i].target_read_id && id == tmp[i].query_read_id)){
            return;
        }
    }
    // Second locate the edge from id to id2, and add it to the vector
    // Sorted by lexographically smaller value
    if(id < id2){
    	tmp.push_back(all_matches[id][id2]);
    } else {
        tmp.push_back(all_matches[id2][id]);
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

int MatchUtils::compute_contigs(std::string id, std::unordered_map<std::string, std::unordered_map<std::string,Match> >& all_matches, std::unordered_map<std::string,std::vector<std::string> >& read_indegree, std::unordered_map<std::string,std::vector<std::string> >& read_outdegree, std::unordered_map<int, std::vector<Match> >& contig_unordered_map, int contig_number){
    // Look at indegree and outdegree of id
    // If it is a branch or a dead end, iterate until we get to a branch
    if((read_indegree[id].size() >= 2) || (read_outdegree[id].size() == 0 && read_indegree[id].size() >=1)){
        // Iterate over each branch
        //std::cerr << "Indegree of 2 or more or dead end on outdegree for " << id << std::endl;
        int count_already_matched = 0;
        for(std::size_t i = 0; i < read_indegree[id].size(); i++){
            // Check if branch has already been taken from other direction
            if(!check_contigs_for_match(id, read_indegree[id][i], contig_unordered_map)){
                // Get vector of overlaps and save as a contig
                std::vector<Match> tmp;
                //std::cerr << "\t"<< id << " " << read_indegree[id][i] << std::endl;
                get_matches_for_contig(id, read_indegree[id][i], all_matches, read_indegree, read_outdegree, tmp);
                //std::cerr << tmp.size() << std::endl;
                contig_unordered_map.insert(make_pair(contig_number, tmp));
                contig_number++;
            } else {
                count_already_matched++;
            }
        }
        // Have the case where the outdegree can be > 0.
        if(read_outdegree[id].size() >= 1){
            for(std::size_t i = 0; i < read_outdegree[id].size(); i++){
                // Check if branch has already been taken from other direction
                if(!check_contigs_for_match(id, read_outdegree[id][i], contig_unordered_map)){
                    // Get vector of overlaps and save as a contig
                    std::vector<Match> tmp;
                    //std::cerr << "\t"<< id << " " << read_outdegree[id][i] << std::endl;
                    get_matches_for_contig(id, read_outdegree[id][i], all_matches, read_indegree, read_outdegree, tmp);
                    //std::cerr << tmp.size() << std::endl;
                    contig_unordered_map.insert(make_pair(contig_number, tmp));
                    contig_number++;
                }
            }
        }
        // Here we have the case where a read is a branch and dead end. ie the end of a bubble. But it needs to be its own contig
        if((read_outdegree[id].size() == 0 && (int) read_indegree[id].size() >=2) && count_already_matched == (int) read_indegree[id].size()){
            std::string id2 = read_indegree[id][0];
            std::vector<Match> tmp;
            Match tmp_match;
            if(id < id2){
                tmp_match = all_matches[id][id2];
            } else {
				tmp_match = all_matches[id2][id];
			}
            // Want to make sure we use the correct length
            if(id == tmp_match.target_read_id){
                tmp_match.length_to_use = tmp_match.target_read_length;
            } else {
                tmp_match.length_to_use = tmp_match.query_read_length;
            }
            tmp.push_back(tmp_match);
            contig_unordered_map.insert(make_pair(contig_number, tmp));
            contig_number++;
        }
    }
    if((read_outdegree[id].size() >= 2) || (read_indegree[id].size() == 0 && read_outdegree[id].size() >=1)){
        //std::cerr << "Outdegree of 2 or more or dead end on indegree for " << id << std::endl;
        int count_already_matched = 0;
        for(std::size_t i = 0; i < read_outdegree[id].size(); i++){
            // Check if branch has already been taken from other direction
            if(!check_contigs_for_match(id, read_outdegree[id][i], contig_unordered_map)){
                // Get vector of overlaps and save as a contig
                std::vector<Match> tmp;
                //std::cerr << "\t"<< id << " " << read_outdegree[id][i] << std::endl;
                get_matches_for_contig(id, read_outdegree[id][i], all_matches, read_indegree, read_outdegree, tmp);
                contig_unordered_map.insert(make_pair(contig_number, tmp));
                //std::cerr << tmp.size() << std::endl;
                contig_number++;
            } else {
                count_already_matched++;
            }
        }
        // Have the case where the outdegree can be > 0.
        if(read_indegree[id].size() >= 1){
            for(std::size_t i = 0; i < read_indegree[id].size(); i++){
                // Check if branch has already been taken from other direction
                if(!check_contigs_for_match(id, read_indegree[id][i], contig_unordered_map)){
                    // Get vector of overlaps and save as a contig
                    std::vector<Match> tmp;
                    //std::cerr << "\t"<< id << " " << read_indegree[id][i] << std::endl;
                    get_matches_for_contig(id, read_indegree[id][i], all_matches, read_indegree, read_outdegree, tmp);
                    contig_unordered_map.insert(make_pair(contig_number, tmp));
                    //std::cerr << tmp.size() << std::endl;
                    contig_number++;
                }
            }
        }
        if((read_indegree[id].size() == 0 && (int) read_outdegree[id].size() >=2) && count_already_matched == (int) read_outdegree[id].size()){
            std::string id2 = read_outdegree[id][0];
            std::vector<Match> tmp;
            Match tmp_match;
            if(id < id2){
                tmp_match = all_matches[id][id2];
            } else {
                tmp_match = all_matches[id2][id];
            }
            // Want to make sure we use the correct length
            if(id == tmp_match.target_read_id){
                tmp_match.length_to_use = tmp_match.target_read_length;
            } else {
                tmp_match.length_to_use = tmp_match.query_read_length;
            }
            tmp.push_back(tmp_match);
            contig_unordered_map.insert(make_pair(contig_number, tmp));
            contig_number++;
        }
    }
    return contig_number;
}

std::vector<Match> get_all_matches_for_node(std::unordered_map<std::string, std::unordered_map<std::string, Match> >& all_matches, std::unordered_map<std::string, std::vector<std::string>>& matches_indexed, std::string v){
    std::vector<Match> tmp;
    // First get all v
    for (std::unordered_map<std::string, Match>::iterator it= all_matches[v].begin(); it!=all_matches[v].end(); ++it){
    	tmp.push_back(it->second);
    }
    for(std::size_t j = 0; j < matches_indexed[v].size(); j++){
        tmp.push_back(all_matches[matches_indexed[v][j]][v]);
    }
    return tmp;
}

void MatchUtils::reduce_edges(std::unordered_map<std::string, std::unordered_map<std::string, Match> >& all_matches, std::unordered_map<std::string, std::vector<std::string>>& matches_indexed, std::unordered_set<std::string>& read_ids, int fuzz)
{
	//int fuzz = 3500;
	int count = 0;
	std::unordered_map<std::string, int> mark;
    // Myers Transitive Reduction Alg
    //Mark each node as vacant (0) & reduce = false
    for (std::unordered_set<std::string>::iterator it=read_ids.begin(); it!=read_ids.end(); ++it){
        mark.insert(std::pair<std::string,int>(*it, 0));
    }

    for (std::unordered_set<std::string>::iterator it=read_ids.begin(); it!=read_ids.end(); ++it){
        std::string v = *it;
        // Get a vector that has all the edges that come in/out of v (bidirected edges)
        //std::vector<Match> v_to_w = get_all_matches_for_node(all_matches, matches_indexed ,v);
        int flongest = 0;
        int rlongest = 0;
        //for (std::size_t i = 0; i < v_to_w.size(); ++i)
        //{	
        	// First compute the suffix/edge length and store it in the match, relative to a specific read
        	// Needed so we can sort by increasing distance
        	//v_to_w[i].set_length(v);
        	//v_to_w[i].set_orientation(v);
            // Set all w for each edge v->w as in play
            //if(v_to_w[i].reduce){
            //	std::cerr << "Found Already reduced edge from: " << v_to_w[i].query_read_id << " to " << v_to_w[i].target_read_id << std::endl;
            //	continue;
            //}
        //	if(v == v_to_w[i].query_read_id){
        //    	mark[v_to_w[i].target_read_id] = 1;
       // 	} else {
       // 		mark[v_to_w[i].query_read_id] = 1;
        //	}
        //}
        for (std::unordered_map<std::string, Match>::iterator it2= all_matches[v].begin(); it2!=all_matches[v].end(); ++it2){
        	if(v == it2->second.query_read_id){
            	mark[it2->second.target_read_id] = 1;
        	} else {
        		mark[it2->second.query_read_id] = 1;
        	}
        }
        for(std::size_t j = 0; j < matches_indexed[v].size(); j++){
        	if(v == all_matches[matches_indexed[v][j]][v].query_read_id){
            	mark[all_matches[matches_indexed[v][j]][v].target_read_id] = 1;
        	} else {
        		mark[all_matches[matches_indexed[v][j]][v].query_read_id] = 1;
        	}
    	}
        //Match::sort_matches(v_to_w);
        //for (std::size_t i = 0; i < v_to_w.size(); ++i)
        //{
        	//if(v_to_w[i].reduce){
        	//	continue;
        	//}
        //	if(v_to_w[i].length > flongest && v_to_w[i].orientation > 0){
        //        flongest = v_to_w[i].length;
        //    }
         //   if(v_to_w[i].length > rlongest && v_to_w[i].orientation < 0){
         //       rlongest = v_to_w[i].length;
         //   }
        //}
        for (std::unordered_map<std::string, Match>::iterator it2= all_matches[v].begin(); it2!=all_matches[v].end(); ++it2){
        	if(v == it2->second.query_read_id){
	        	if(it2->second.length > flongest && it2->second.orientation > 0){
	                flongest = it2->second.length;
	            }
	            if(it2->second.length > rlongest && it2->second.orientation < 0){
	                rlongest = it2->second.length;
	            }
	        } else{
	        	if(it2->second.length > flongest && it2->second.orientation > 0){
	                flongest = it2->second.prefix_length;
	            }
	            if(it2->second.length > rlongest && it2->second.orientation < 0){
	                rlongest = it2->second.prefix_length;
	            }
	        }
        }
        for(std::size_t j = 0; j < matches_indexed[v].size(); j++){
        	if(v == all_matches[matches_indexed[v][j]][v].query_read_id){
	        	if(all_matches[matches_indexed[v][j]][v].length > flongest && all_matches[matches_indexed[v][j]][v].orientation > 0){
	                flongest = all_matches[matches_indexed[v][j]][v].length;
	            }
	            if(all_matches[matches_indexed[v][j]][v].length > rlongest && all_matches[matches_indexed[v][j]][v].orientation < 0){
	                rlongest = all_matches[matches_indexed[v][j]][v].length;
	            }
	        } else {
	        	if(all_matches[matches_indexed[v][j]][v].length > flongest && all_matches[matches_indexed[v][j]][v].orientation > 0){
	                flongest = all_matches[matches_indexed[v][j]][v].prefix_length;
	            }
	            if(all_matches[matches_indexed[v][j]][v].length > rlongest && all_matches[matches_indexed[v][j]][v].orientation < 0){
	                rlongest = all_matches[matches_indexed[v][j]][v].prefix_length;
	            }
	        }
    	}
        flongest += fuzz;
        rlongest += fuzz;

        for (std::unordered_map<std::string, Match>::iterator it2= all_matches[v].begin(); it2!=all_matches[v].end(); ++it2)
        {
        	int longest = 0;
        	if(it2->second.orientation < 0){
        		longest = rlongest;
        	} else {
        		longest = flongest;
        	}
            // if w in v->w in play
            std::string w;
            if(it2->second.query_read_id == v){
            	w = it2->second.target_read_id;
            } else {
            	w  = it2->second.query_read_id;
            }
            if(mark[w] == 1){
                for (std::unordered_map<std::string, Match>::iterator it3= all_matches[w].begin(); it3!=all_matches[w].end(); ++it3)
                {
		            if(it2->second.strand == '-'){
		            	// w is now flipped, so any w to x that is + is really - in the eyes of v 
		            	if(it3->second.orientation == it2->second.orientation){
		            		continue;
		            	}
		            } else if(it3->second.orientation != it2->second.orientation){
		            	continue;
		            }
                	if(w == it3->second.query_read_id){
                		if(v == it3->second.target_read_id){
		            		//Same edge just reversed, Don't need to remove yourself
		            		continue;
		            	}
		            	if(it2->second.length + it3->second.length <= longest){
		                	if(mark[it3->second.target_read_id] == 1){
		                    	mark[it3->second.target_read_id] = -1;
		                	}
		            	}
                	} else {
                		if(v == it3->second.query_read_id){
		            		//Same edge just reversed, Don't need to remove yourself
		            		continue;
		            	}
		            	if(it2->second.length + it3->second.prefix_length <= longest){
		                	if(mark[it3->second.query_read_id] == 1){
		                    	mark[it3->second.query_read_id] = -1;
		                	}
		            	}
                	}
                }

                for(std::size_t j = 0; j < matches_indexed[w].size(); j++)
                {
		            if(it2->second.strand == '-'){
		            	// w is now flipped, so any w to x that is + is really - in the eyes of v 
		            	if(all_matches[matches_indexed[w][j]][w].orientation == it2->second.orientation){
		            		continue;
		            	}
		            } else if(all_matches[matches_indexed[w][j]][w].orientation != it2->second.orientation){
		            	continue;
		            }
                	if(w == all_matches[matches_indexed[w][j]][w].query_read_id){
                		if(v == all_matches[matches_indexed[w][j]][w].target_read_id){
		            		//Same edge just reversed, Don't need to remove yourself
		            		continue;
		            	}
		            	if(it2->second.length + all_matches[matches_indexed[w][j]][w].length <= longest){
		                	if(mark[all_matches[matches_indexed[w][j]][w].target_read_id] == 1){
		                    	mark[all_matches[matches_indexed[w][j]][w].target_read_id] = -1;
		                	}
		            	} 
                	} else {
                		if(v == all_matches[matches_indexed[w][j]][w].query_read_id){
		            		//Same edge just reversed, Don't need to remove yourself
		            		continue;
		            	}
		            	if(it2->second.length + all_matches[matches_indexed[w][j]][w].prefix_length <= longest){
		                	if(mark[all_matches[matches_indexed[w][j]][w].query_read_id] == 1){
		                    	mark[all_matches[matches_indexed[w][j]][w].query_read_id] = -1;
		                	}
		            	}
                	}
                }
            }
        }

        for(std::size_t i = 0; i < matches_indexed[v].size(); i++)
        {
        	int longest = 0;
        	if(all_matches[matches_indexed[v][i]][v].orientation < 0){
        		longest = rlongest;
        	} else {
        		longest = flongest;
        	}
            // if w in v->w in play
            std::string w;
            if(all_matches[matches_indexed[v][i]][v].query_read_id == v){
            	w = all_matches[matches_indexed[v][i]][v].target_read_id;
            } else {
            	w  = all_matches[matches_indexed[v][i]][v].query_read_id;
            }
            if(mark[w] == 1){
                for (std::unordered_map<std::string, Match>::iterator it3= all_matches[w].begin(); it3!=all_matches[w].end(); ++it3)
                {
		            if(all_matches[matches_indexed[v][i]][v].strand == '-'){
		            	// w is now flipped, so any w to x that is + is really - in the eyes of v 
		            	if(it3->second.orientation == all_matches[matches_indexed[v][i]][v].orientation){
		            		continue;
		            	}
		            } else if(it3->second.orientation != all_matches[matches_indexed[v][i]][v].orientation){
		            	continue;
		            }
                	if(w == it3->second.query_read_id){
                		//std::cerr << v << " to " << w_to_x[j].query_read_id << " " << w_to_x[j].target_read_id << std::endl;
                		if(v == it3->second.target_read_id){
		            		//Same edge just reversed, Don't need to remove yourself
		                  //  std::cerr << "Reversed " << v << std::endl;
		            		continue;
		            	}
		                //std::cerr << "For " << v << " Check " << w << " to " << all_matches[w][j].target_read_id << std::endl;
		            	if(all_matches[matches_indexed[v][i]][v].length + it3->second.length <= longest){
		                	if(mark[it3->second.target_read_id] == 1){
		                    	mark[it3->second.target_read_id] = -1;
		                    //    std::cerr << "A1. Marking: " << v << " to " << w_to_x[j].target_read_id << " For reduction" << std::endl;
		                	}
		            	}
                	} else {
                		//std::cerr << v << " to " << w_to_x[j].target_read_id << " " << w_to_x[j].query_read_id << std::endl;
                		if(v == it3->second.query_read_id){
		            		//Same edge just reversed, Don't need to remove yourself
		                    //std::cerr << "Reversed " << v << std::endl;
		            		continue;
		            	}
		                //std::cerr << "For " << v << " Check " << w << " to " << all_matches[w][j].target_read_id << std::endl;
		            	if(all_matches[matches_indexed[v][i]][v].length + it3->second.prefix_length <= longest){
		                	if(mark[it3->second.query_read_id] == 1){
		                    	mark[it3->second.query_read_id] = -1;
		                       //std::cerr << "A2. Marking: " << v << " to " << w_to_x[j].query_read_id << " For reduction" << std::endl;
		                	}
		            	}
                	}
                }

                for(std::size_t j = 0; j < matches_indexed[w].size(); j++)
                {
		            if(all_matches[matches_indexed[v][i]][v].strand == '-'){
		            	// w is now flipped, so any w to x that is + is really - in the eyes of v 
		            	if(all_matches[matches_indexed[w][j]][w].orientation == all_matches[matches_indexed[v][i]][v].orientation){
		            		continue;
		            	}
		            } else if(all_matches[matches_indexed[w][j]][w].orientation != all_matches[matches_indexed[v][i]][v].orientation){
		            	continue;
		            }
                	if(w == all_matches[matches_indexed[w][j]][w].query_read_id){
                		//std::cerr << v << " to " << w_to_x[j].query_read_id << " " << w_to_x[j].target_read_id << std::endl;
                		if(v == all_matches[matches_indexed[w][j]][w].target_read_id){
		            		//Same edge just reversed, Don't need to remove yourself
		                  //  std::cerr << "Reversed " << v << std::endl;
		            		continue;
		            	}
		                //std::cerr << "For " << v << " Check " << w << " to " << all_matches[w][j].target_read_id << std::endl;
		            	if(all_matches[matches_indexed[v][i]][v].length + all_matches[matches_indexed[w][j]][w].length <= longest){
		                	if(mark[all_matches[matches_indexed[w][j]][w].target_read_id] == 1){
		                    	mark[all_matches[matches_indexed[w][j]][w].target_read_id] = -1;
		                    //    std::cerr << "A1. Marking: " << v << " to " << w_to_x[j].target_read_id << " For reduction" << std::endl;
		                	}
		            	} 
                	} else {
                		//std::cerr << v << " to " << w_to_x[j].target_read_id << " " << w_to_x[j].query_read_id << std::endl;
                		if(v == all_matches[matches_indexed[w][j]][w].query_read_id){
		            		//Same edge just reversed, Don't need to remove yourself
		                    //std::cerr << "Reversed " << v << std::endl;
		            		continue;
		            	}
		                //std::cerr << "For " << v << " Check " << w << " to " << all_matches[w][j].target_read_id << std::endl;
		            	if(all_matches[matches_indexed[v][i]][v].length + all_matches[matches_indexed[w][j]][w].prefix_length <= longest){
		                	if(mark[all_matches[matches_indexed[w][j]][w].query_read_id] == 1){
		                    	mark[all_matches[matches_indexed[w][j]][w].query_read_id] = -1;
		                       //std::cerr << "A2. Marking: " << v << " to " << w_to_x[j].query_read_id << " For reduction" << std::endl;
		                	}
		            	}
                	}
                }
            }
        }

        count += mark_matches_for_node(all_matches, v, mark);
        for (std::unordered_map<std::string, Match>::iterator it2= all_matches[v].begin(); it2!=all_matches[v].end(); ++it2){
        	if(it2->second.reduce){
            	continue;
        	}
        	if(v == it2->second.query_read_id){
		        mark[it2->second.target_read_id] = 0;
		    } else {
		    	mark[it2->second.query_read_id] = 0;
		    }
        }
        for(std::size_t j = 0; j < matches_indexed[v].size(); j++){
        	if(all_matches[matches_indexed[v][j]][v].reduce){
            	continue;
        	}
        	if(v == all_matches[matches_indexed[v][j]][v].query_read_id){
		        mark[all_matches[matches_indexed[v][j]][v].target_read_id] = 0;
		    } else {
		    	mark[all_matches[matches_indexed[v][j]][v].query_read_id] = 0;
		    }
    	}
    }
    std::cerr << "\tReduced " << count << " edges" << std::endl;
}

void MatchUtils::toGfa(std::unordered_map<std::string, std::unordered_map<std::string,Match> >& all_matches, std::unordered_map<std::string, int>& read_lengths, std::string file_name, std::unordered_map<std::string, std::vector<std::string> >& read_indegree, std::unordered_map<std::string, std::vector<std::string> >& read_outdegree, std::unordered_map<std::string, std::string>& read_names, std::unordered_map<std::string, std::string>& colours, std::unordered_map<std::string, float>& read_coverage)
{
	std::ofstream gfaOutput;
  	gfaOutput.open(file_name);
  	gfaOutput << "H\tVN:Z:Test\n";
  	for (std::unordered_map<std::string, int>::iterator it=read_lengths.begin(); it!=read_lengths.end(); ++it){
        if(read_outdegree[it->first].size() > 0 || read_indegree[it->first].size() > 0){
  		    gfaOutput << "S\t" << read_names[it->first] <<"\t*\tLN:i:" << it->second << "\tKC:i:"<< static_cast<int>(it->second*read_coverage[it->first]) << "\tCL:z:" << colours[it->first] << "\tC2:z:" << colours[it->first] << "\n";
        }
  	}
  	for (std::unordered_map<std::string, std::unordered_map<std::string,Match> >::iterator it=all_matches.begin(); it!=all_matches.end(); ++it){
  		for (std::unordered_map<std::string,Match>::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
  		{
  			//std::cerr << it->second[i].query_read_id << " To " << it->second[i].target_read_id << " " << it->second[i].reduce << std::endl;
  			if(!it2->second.reduce){
  				if(it2->second.strand == '+'){
					if(it2->second.query_read_start > it2->second.target_read_start){
						gfaOutput << "L\t" << read_names[it2->second.query_read_id] << "\t+\t" << read_names[it2->second.target_read_id] << "\t+\t"<< it2->second.cigar <<"\n";
					} else {
						gfaOutput << "L\t" << read_names[it2->second.target_read_id] << "\t+\t" << read_names[it2->second.query_read_id] << "\t+\t"<< it2->second.cigar <<"\n";
					}
				} else {
					if(it2->second.query_read_start > it2->second.query_read_length - it2->second.query_read_end && it2->second.target_read_start > it2->second.target_read_length - it2->second.target_read_end){
						gfaOutput << "L\t" << read_names[it2->second.query_read_id] << "\t+\t" << read_names[it2->second.target_read_id] << "\t-\t"<< it2->second.cigar <<"\n";
					} else {
						gfaOutput << "L\t" << read_names[it2->second.target_read_id] << "\t-\t" << read_names[it2->second.query_read_id] << "\t+\t"<< it2->second.cigar <<"\n";
					}
				}
                
  			}
  		}
  	}
  	gfaOutput.close();
}

