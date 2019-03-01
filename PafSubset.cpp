#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include "MatchUtils.hpp"

using namespace std;

void helper(set<string>& ids, int levels, string source_read, unordered_map<string, unordered_map<string,Match> >& raw_matches)
{
    if(levels >= 0){
        ids.insert(source_read);
        if(levels >= 1){
            // Want to get some number of levels from the source
            for (unordered_map<string, unordered_map<string,Match> >::iterator it=raw_matches.begin(); it!=raw_matches.end(); ++it)
            {
                for (unordered_map<string,Match>::iterator it2=it->second.begin(); it2!=it->second.end(); ++it2)
                {
                    if(it2->second.query_read_id == source_read){
                        if(ids.count(it2->second.target_read_id) == 0){
                            helper(ids, levels -1, it2->second.target_read_id, raw_matches);
                        }
                    }
                    if(it2->second.target_read_id == source_read){
                        if(ids.count(it2->second.query_read_id) == 0){
                            helper(ids, levels -1, it2->second.query_read_id, raw_matches);
                        }
                    }
                }
            }
        }
        // if levels is just 0 only get edges to and from the source read
    }    
}

int main(int argc, char** argv)
{
    if(argc != 4){
        cout << "USAGE: PafSubset file.paf read_id levels" << endl;
        return 0;
    }
    // Read in and Parse input file
	unordered_map<string, unordered_map<string,Match> > all_matches;
    unordered_map<string, unordered_map<string,Match> > raw_matches;
    set<string> read_ids;
    set<string> chimeric_reads;
    unordered_map<string, Read> read_classification;
    unordered_map<string, vector<string>> matches_indexed;
    unordered_map<string, int> read_lengths;
    MatchUtils::read_paf_file(all_matches, matches_indexed, raw_matches, read_ids, read_lengths, argv[1], chimeric_reads, read_classification, true, true);
    cout << raw_matches.size() << endl;
    string source_read = argv[2];
    int levels = atoi(argv[3]);
    set<string> ids;
    helper(ids, levels, source_read, raw_matches);
    for (unordered_map<string, unordered_map<string,Match> >::iterator it=raw_matches.begin(); it!=raw_matches.end(); ++it)
    {
        for (unordered_map<string,Match>::iterator it2=it->second.begin(); it2!=it->second.end(); ++it2)
        {
            //Iterate over every edge, if edge contains an id in id set then output it in paf form
            if(find(ids.begin(), ids.end(), it2->second.query_read_id) != ids.end() || find(ids.begin(), ids.end(), it2->second.target_read_id) != ids.end()) {
                //Output data in paf format with cigar strings if possible
                cout << it2->second.query_read_id << "\t" << it2->second.query_read_length << "\t" << it2->second.query_read_start << "\t" << it2->second.query_read_end << "\t";
                cout << it2->second.strand << "\t" << it2->second.target_read_id << "\t" << it2->second.target_read_length << "\t" << it2->second.target_read_start << "\t";
                cout << it2->second.target_read_end << "\t" << it2->second.residue_matches << "\t" << it2->second.alignment_block << "\t0\tcg:Z:" << it2->second.cigar << endl;
            }
        }
    }
    return 0;
}
