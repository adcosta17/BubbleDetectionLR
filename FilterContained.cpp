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
#include <iomanip>
#include <string>
#include <stdexcept>

namespace io = boost::iostreams;

using namespace std;

string compress_string(const string& str, int compressionlevel = Z_BEST_COMPRESSION) {
    z_stream zs;                        // z_stream is zlib's control structure
    memset(&zs, 0, sizeof(zs));

    if (deflateInit(&zs, compressionlevel) != Z_OK)
        throw(runtime_error("deflateInit failed while compressing."));

    zs.next_in = (Bytef*)str.data();
    zs.avail_in = str.size();           // set the z_stream's input

    int ret;
    char outbuffer[32768];
    string outstring;

    // retrieve the compressed bytes blockwise
    do {
        zs.next_out = reinterpret_cast<Bytef*>(outbuffer);
        zs.avail_out = sizeof(outbuffer);

        ret = deflate(&zs, Z_FINISH);

        if (outstring.size() < zs.total_out) {
            // append the block to the output string
            outstring.append(outbuffer,
                             zs.total_out - outstring.size());
        }
    } while (ret == Z_OK);

    deflateEnd(&zs);

    if (ret != Z_STREAM_END) {          // an error occurred that was not EOF
        ostringstream oss;
        oss << "Exception during zlib compression: (" << ret << ") " << zs.msg;
        throw(runtime_error(oss.str()));
    }

    return outstring;
}

int main(int argc, char** argv)
{
    if(argc != 4){
        cerr << "USAGE: FilterContained input_dir_of_contained_read_lists input_file.paf.gz output_file.gz" << endl;
        return 0;
    }
  
    set<string> overlaps;
    string file_name = argv[1];

    set<string> txt_files;
    DIR *dir;
    struct dirent *ent;
    string suffix = ".txt";
    if ((dir = opendir (file_name.c_str())) != NULL) {
      /* print all the files and directories within directory */
      while ((ent = readdir (dir)) != NULL) {
        string name(ent->d_name);
        //cerr << name << endl;
        if (name.find(suffix) != string::npos){
            txt_files.insert(file_name + "/" + name);
            //cerr << file_name + "/" + name << endl;
        }
      }
      closedir (dir);
    } else {
      /* could not open directory */
      cerr << "ERROR: Could not open " << file_name << endl; 
      return 0;
    }
    cerr << "Found " << txt_files.size() << " Contained Read Files" << endl;

    // Take the list of paf_files and then for each of them read in the file
    cerr << "Reading in Contained Reads" << endl;
    for (set<string>::iterator it = txt_files.begin(); it != txt_files.end(); ++it) {
    	string tmp = *it;
        size_t found = tmp.find_last_of("/");
        string name(tmp.substr(found+1));
        cerr << "\tReading in " << name << endl;

        ifstream inputFile(tmp);
        string line;
        while (getline(inputFile, line))
        {   
            istringstream lin(line);
            string sid;
            lin >> sid;
            overlaps.insert(compress_string(sid));
        }
    }

    string output_file = argv[3];
    gzFile tmp_file = gzopen(output_file.c_str(),"ab");

    io::filtering_istream in;
    in.push(io::gzip_decompressor());
    in.push(io::file_source(argv[2]));
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
        string id1 = compress_string(c1);
        string id2 = compress_string(c6);
            
        if(overlaps.count(id1) == 0 && overlaps.count(id2) == 0){
            string val = line +"\n";
            //cout << line << endl;
            gzprintf(tmp_file,val.c_str());
        }
    }
    gzclose_w(tmp_file);
    cerr << "Written Output to " << output_file << endl;
    return 0;
}

