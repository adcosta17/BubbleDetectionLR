#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include <set>
#include <zlib.h>
#include <dirent.h>
#include <iomanip>
#include <string>
#include <stdexcept>
#include <map>
#include <vector>
#include <iterator>
#include "Read_Alignment.hpp"

using namespace std;


// take from http://stackoverflow.com/a/236803/248823
void split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
}

int main(int argc, char** argv) {

	if(argc != 5){
        cerr << "USAGE: ProcessBubbleList inputBubbleList Prefix1 Prefix2 outputfile" << endl;
        return 0;
    }

    map<string, Read_Alignment> all_alignments;

    // Read Pafs first. Create a Read_Alignment for each read
    // Then as you read the paf create an alignment for each line

    string input = "" + string(argv[2])+ string(".paf");
    ifstream inputFilePaf1(input);
    string line;
    while (getline(inputFilePaf1, line))
    {
        istringstream lin(line);
        string id, c6, meta, cg;
        char c5;
        int c2, c3, c4, c7, c8, c9, c10, c11, qual;
        lin >> id >> c2 >> c3 >> c4 >> c5 >> c6 >> c7 >> c8 >> c9 >> c10 >> c11 >> qual;

        Alignment tmp_aln(id, c11, c3, c4, qual);

        if(all_alignments.count(id) == 0){
            Read_Alignment tmp(id, c2);
            vector<Alignment> tmp_vec;
            tmp.species1_aln = tmp_vec;
            all_alignments.insert(make_pair(id, tmp));
        }
        all_alignments.at(id).species1_aln.push_back(tmp_aln);
    }

    input = "" + string(argv[3])+ string(".paf");
    ifstream inputFilePaf2(input);
    while (getline(inputFilePaf2, line))
    {
        istringstream lin(line);
        string id, c6, meta, cg;
        char c5;
        int c2, c3, c4, c7, c8, c9, c10, c11, qual;
        lin >> id >> c2 >> c3 >> c4 >> c5 >> c6 >> c7 >> c8 >> c9 >> c10 >> c11 >> qual;

        Alignment tmp_aln(id, c11, c3, c4, qual);

        if(all_alignments.count(id) == 0){
            Read_Alignment tmp(id, c2);
            vector<Alignment> tmp_vec;
            tmp.species2_aln = tmp_vec;
            all_alignments.insert(make_pair(id, tmp));
        }
        all_alignments.at(id).species2_aln.push_back(tmp_aln);
    }

    input = "" + string(argv[2]) + string(".sam");
    ifstream inputFileSam1(input);
    while (getline(inputFileSam1, line))
    {
        istringstream lin(line);
        string id;
        int mapscore;
        lin >> id >> mapscore;
        if(mapscore >= 2048){
            all_alignments.at(id).split_1 = true;
        }
    }

    input = "" + string(argv[3])+ string(".sam");
    ifstream inputFileSam2(input);
    while (getline(inputFileSam2, line))
    {
        istringstream lin(line);
        string id;
        int mapscore;
        lin >> id >> mapscore;
        if(mapscore >= 2048){
            all_alignments.at(id).split_2 = true;
        }
    }

    ofstream output;
    output.open(argv[4]);

    ifstream infile(argv[1]);
    while (getline(infile, line))
    {
        vector<string> row;
        split(line, '\t', row);
        float score = stof(row[0]);
        int num_reads = stoi(row[1]);
        int diff_tax = stoi(row[5]);

        cout << score << endl;

        string read;
        stringstream data1(row[10]);
        vector<string> arm1;
        while(getline(data1,read,' '))
        {
            arm1.push_back(read);
        }

        stringstream data2(row[11]);
        vector<string> arm2;
        while(getline(data2,read,' '))
        {
            arm2.push_back(read);
        }

        // Want to count the number of reads in each arm that mapped uniquely to each species and the ones that mapped to both or neither
        // Use 30 as cutoff for alignment
        int arm1_uniuqe_s1 = 0;
        int arm1_uniuqe_s2 = 0;
        int arm2_uniuqe_s1 = 0;
        int arm2_uniuqe_s2 = 0;
        int arm1_both = 0;
        int arm1_neither = 0;
        int arm2_neither = 0;
        int arm2_both = 0;
        int start_in_s1 = 0;
        int start_in_s2 = 0;
        int end_in_s1 = 0;
        int end_in_s2 = 0;

        for(std::size_t i = 0; i < arm1.size(); i++){
            bool s1 = false;
            bool s2 = false;
            // check arm[i] against species 1
            if(all_alignments.at(arm1[i]).align_species_1()){
                s1 = true;
            }
            if(all_alignments.at(arm1[i]).align_species_2()){
                s2 = true;
            }
            if(s1 && s2){
                arm1_both += 1;
            } else if(s1){
                arm1_uniuqe_s1 += 1;
            } else if(s2){
                arm1_uniuqe_s2 += 1;
            } else {
                arm1_neither += 1;
            }
            if(i == 0){
                if(s1){
                    start_in_s1 = 1;
                }
                if(s2){
                    start_in_s2 = 1;
                }
            } 
            if(i == arm1.size() - 1){
                if(s1){
                    end_in_s1 = 1;
                }
                if(s2){
                    end_in_s2 = 1;
                }
            }
        }

        for(std::size_t i = 0; i < arm2.size(); i++){
            bool s1 = false;
            bool s2 = false;
            if(all_alignments.at(arm2[i]).align_species_1()){
                s1 = true;
            }
            if(all_alignments.at(arm2[i]).align_species_2()){
                s2 = true;
            }
            if(s1 && s2){
                arm2_both += 1;
            } else if(s1){
                arm2_uniuqe_s1 += 1;
            } else if(s2){
                arm2_uniuqe_s2 += 1;
            } else {
                arm2_neither += 1;
            }
        }


        // bubble start and end will be same for both arm
        output << score << "\t" << num_reads << "\t" << diff_tax << "\t" << arm1_uniuqe_s1 << "\t" << arm1_uniuqe_s2 << "\t" << arm1_both << "\t" << arm1_neither << "\t" << arm2_uniuqe_s1 << "\t" << arm2_uniuqe_s2 << "\t" << arm2_both << "\t" << arm2_neither << "\t" << start_in_s1 << "\t" << start_in_s2 << "\t" << end_in_s1 << "\t" << end_in_s2 <<"\n";
     }

     output.close();

    return 0;
}
