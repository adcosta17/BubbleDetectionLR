#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include <set>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <zlib.h>
#include <dirent.h>
namespace io = boost::iostreams;

using namespace std;


int main(int argc, char** argv)
{
    if(argc != 3){
        cout << "USAGE: CombinePaf input_paf_dir output_file.paf" << endl;
        return 0;
    }
  
    set<string> overlaps;
    string file_name = argv[1];

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
    cerr << "Combining files" << endl;
    int count = 0;
    string output_file = argv[2];
    gzFile tmp_file = gzopen(output_file.c_str(),"ab");
    for (set<string>::iterator it = paf_files.begin(); it != paf_files.end(); ++it) {

        // get species name
        string tmp = *it;
        size_t found = tmp.find_last_of("/");
        string name(tmp.substr(found+1).substr(0,tmp.substr(found+1).size()-11));
        cerr << "Combining " << name << endl;

        io::filtering_istream in;
        in.push(io::gzip_decompressor());
        in.push(io::file_source(tmp));
        vector<int> sizes;
        string line;
        while (getline(in, line, '\n'))
        {
            // Each line of input is split into columns
            // Each line coresponds to a called overlap
            istringstream lin(line);
            string c1, c6;
            char c5;
            int c2, c3, c4;
            lin >> c1 >> c2 >> c3 >> c4 >> c5 >> c6;
            string tmp;
            if(c1 < c6){
                tmp = c1 + "_:_" + c6;
            } else {
                tmp = c6 + "_:_" + c1;
            }
            //cerr << tmp << endl;
            if(overlaps.count(tmp) == 0){
                overlaps.insert(tmp);
                string val = line +"\n";
                //cerr << line << endl;
                gzprintf(tmp_file,val.c_str());
            }
        }
    }
    gzclose_w(tmp_file);
    cerr << "Written Output to " << output_file << endl;
    return 0;
}
