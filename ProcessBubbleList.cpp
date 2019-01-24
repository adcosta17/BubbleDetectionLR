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
        cerr << "USAGE: ProcessBubbleList inputBubbleList Sam1 Sam2 outputfile" << endl;
        return 0;
    }

    // Read Sams and store the best map score for each read
    map<string, int> sam1;
    map<string, int> sam2;

    ifstream inputFileSam1(argv[2]);
    string line;
    while (getline(inputFileSam1, line))
    {
        istringstream lin(line);
        string id, tmp1, tmp2, tmp3;
        int mapscore;
        lin >> id >> tmp1 >> tmp2 >> tmp3 >> mapscore;
        if(sam1.count(id) == 0){
            sam1.insert(make_pair(id, mapscore));
        } else if(sam1[id] < mapscore){
            sam1[id] = mapscore;
        }
        cout << sam1[id] << endl;
    }

    ifstream inputFileSam2(argv[3]);
    while (getline(inputFileSam2, line))
    {
        istringstream lin(line);
        string id, tmp1, tmp2, tmp3;
        int mapscore;
        lin >> id >> tmp1 >> tmp2 >> tmp3 >> mapscore;
        if(sam2.count(id) == 0){
            sam2.insert(make_pair(id, mapscore));
        } else if(sam2[id] < mapscore){
            sam2[id] = mapscore;
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

        // Compute the average mapscore for arm1 and arm2 for each sam
        float a1ms1 =0;
        float a1ms2 =0;
        float a2ms1 =0;
        float a2ms2 =0;
        for(int i = 0; i < arm1.size() -1; i++){
            if(i == 0) {
                continue;
            }
            if(sam1.count(arm1[i]) > 0){
                a1ms1 += static_cast<float>(sam1[arm1[i]]);
            }
            if(sam2.count(arm1[i]) > 0){
                a1ms2 += static_cast<float>(sam2[arm1[i]]);
            }
        }
        if(arm1.size()-2 > 0){
            a1ms1 = a1ms1/(arm1.size()-2);
            a1ms2 = a1ms2/(arm1.size()-2);
        }

        for(int i = 0; i < arm2.size() -1; i++){
            if(i == 0) {
                continue;
            }
            if(sam1.count(arm2[i]) > 0){
                a2ms1 += static_cast<float>(sam1[arm2[i]]);
            }
            if(sam2.count(arm2[i]) > 0){
                a2ms2 += static_cast<float>(sam2[arm2[i]]);
            }
        }
        if(arm2.size()-2){
            a2ms1 = a2ms1/(arm2.size()-2);
            a2ms2 = a2ms2/(arm2.size()-2);
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

        for(int i = 0; i < arm1.size(); i++){
            bool s1 = false;
            bool s2 = false;
            if(sam1.count(arm1[i]) > 0 && sam1[arm1[i]] >= 30){
                s1 = true;
            }
            if(sam2.count(arm1[i]) > 0 && sam2[arm1[i]] >= 30){
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

        for(int i = 0; i < arm2.size(); i++){
            bool s1 = false;
            bool s2 = false;
            if(sam1.count(arm2[i]) > 0 && sam1[arm2[i]] >= 30){
                s1 = true;
            }
            if(sam2.count(arm2[i]) > 0 && sam2[arm2[i]] >= 30){
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
        output << score << "\t" << num_reads << "\t" << a1ms1 << "\t" << a1ms2 << "\t" << a2ms1 << "\t" << a2ms2 << "\t" << arm1_uniuqe_s1 << "\t" << arm1_uniuqe_s2 << "\t" << arm1_both << "\t" << arm1_neither << "\t" << arm2_uniuqe_s1 << "\t" << arm2_uniuqe_s2 << "\t" << arm2_both << "\t" << arm2_neither << "\t" << start_in_s1 << "\t" << start_in_s2 << "\t" << end_in_s1 << "\t" << end_in_s2 <<"\n";
     }

     output.close();

    return 0;
}
