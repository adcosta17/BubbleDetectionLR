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

/** Compress a STL string using zlib with given compression level and return
  * the binary data. */
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

/** Decompress an STL string using zlib and return the original data. */
string decompress_string(const string& str)
{
    z_stream zs;                        // z_stream is zlib's control structure
    memset(&zs, 0, sizeof(zs));

    if (inflateInit(&zs) != Z_OK)
        throw(runtime_error("inflateInit failed while decompressing."));

    zs.next_in = (Bytef*)str.data();
    zs.avail_in = str.size();

    int ret;
    char outbuffer[32768];
    string outstring;

    // get the decompressed bytes blockwise using repeated calls to inflate
    do {
        zs.next_out = reinterpret_cast<Bytef*>(outbuffer);
        zs.avail_out = sizeof(outbuffer);

        ret = inflate(&zs, 0);

        if (outstring.size() < zs.total_out) {
            outstring.append(outbuffer,
                             zs.total_out - outstring.size());
        }

    } while (ret == Z_OK);

    inflateEnd(&zs);

    if (ret != Z_STREAM_END) {          // an error occurred that was not EOF
        ostringstream oss;
        oss << "Exception during zlib decompression: (" << ret << ") "
            << zs.msg;
        throw(runtime_error(oss.str()));
    }

    return outstring;
}



int main(int argc, char** argv)
{
    if(argc != 3){
        cerr << "USAGE: CombinePaf input_paf_dir output_file.paf" << endl;
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
            string ids;
            if(c1 < c6){
                ids = c1 + "_:_" + c6;
            } else {
                ids = c6 + "_:_" + c1;
            }
            ids = compress_string(ids);
            //cerr << ids << endl;
            if(overlaps.count(ids) == 0){
                overlaps.insert(ids);
                string val = line +"\n";
                //cout << line << endl;
                gzprintf(tmp_file,val.c_str());
            }
        }
    }
    gzclose_w(tmp_file);
    cerr << "Written Output to " << output_file << endl;
    return 0;
}

