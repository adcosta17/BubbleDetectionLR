#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include "MatchUtils.hpp"

using namespace std;

void helper(set<string>& ids, int levels, string source_read, map<string, vector<Match> >& raw_matches)
{
    if(levels >= 0){
        ids.insert(source_read);
        if(levels >= 1){
            // Want to get some number of levels from the source
            for (map<string, vector<Match> >::iterator it=raw_matches.begin(); it!=raw_matches.end(); ++it)
            {
                for (int i = 0; i < it->second.size(); ++i)
                {
                    if(it->second[i].query_read_id == source_read){
                        if(ids.count(it->second[i].target_read_id) == 0){
                            helper(ids, levels -1, it->second[i].target_read_id, raw_matches);
                        }
                    }
                    if(it->second[i].target_read_id == source_read){
                        if(ids.count(it->second[i].query_read_id) == 0){
                            helper(ids, levels -1, it->second[i].query_read_id, raw_matches);
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
	map<string, vector<Match> > all_matches;
    map<string, vector<Match> > raw_matches;
    set<string> read_ids;
    set<string> chimeric_reads;
    map<string, Read> read_classification;
    map<string, vector<string>> matches_indexed;
    map<string, int> read_lengths;
    MatchUtils::read_paf_file(all_matches, matches_indexed, raw_matches, read_ids, read_lengths, argv[1], chimeric_reads, read_classification, true, true);
    cout << raw_matches.size() << endl;
    string source_read = argv[2];
    int levels = atoi(argv[3]);
    set<string> ids;
    helper(ids, levels, source_read, raw_matches);
    for (map<string, vector<Match> >::iterator it=raw_matches.begin(); it!=raw_matches.end(); ++it)
    {
        for (int i = 0; i < it->second.size(); ++i)
        {
            //Iterate over every edge, if edge contains an id in id set then output it in paf form
            if(find(ids.begin(), ids.end(), it->second[i].query_read_id) != ids.end() || find(ids.begin(), ids.end(), it->second[i].target_read_id) != ids.end()) {
                //Output data in paf format with cigar strings if possible
                cout << it->second[i].query_read_id << "\t" << it->second[i].query_read_length << "\t" << it->second[i].query_read_start << "\t" << it->second[i].query_read_end << "\t";
                cout << it->second[i].strand << "\t" << it->second[i].target_read_id << "\t" << it->second[i].target_read_length << "\t" << it->second[i].target_read_start << "\t";
                cout << it->second[i].target_read_end << "\t" << it->second[i].residue_matches << "\t" << it->second[i].alignment_block << "\t0\tcg:Z:" << it->second[i].cigar << endl;
            }
        }
    }
    return 0;
}
