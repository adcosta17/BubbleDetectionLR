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

#include "Match.hpp"

using namespace std;

int main(int argc, char** argv)
{
    if(argc != 2){
        cout << "USAGE: ComputeContained file.paf.gz" << endl;
        return 0;
    }

    io::filtering_istream in_filter;
    in_filter.push(io::gzip_decompressor());
    in_filter.push(io::file_source(argv[1]));

	//ifstream inputFile_filter(file_name);
	string line;
	while (getline(in_filter, line, '\n'))
	{
		istringstream lin(line);
    	string c1, c6, meta, cg;
    	char c5;
    	int c2, c3, c4, c7, c8, c9, c10, c11;
    	lin >> c1 >> c2 >> c3 >> c4 >> c5 >> c6 >> c7 >> c8 >> c9 >> c10 >> c11;
        Match tmpLine(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,0,0,"");
        if(tmpLine.internal_edge()){
        	continue;
        }
        int contained = tmpLine.check_match_contained();
	    if(contained == 1){
	        cout << tmpLine.target_read_id << endl;
	    } else if (contained == -1 ) {
	        cout << tmpLine.query_read_id << endl;
	    }
    }

    return 0;
}
