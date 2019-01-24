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
#include <map>

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
        cerr << "USAGE: ReadCoverageCombined input_dir_of_read_cov_lists output_file.txt" << endl;
        return 0;
    }

    map<string, float> running_cov;

    string file_name = argv[1];

    set<string> txt_files;
    DIR *dir;
    struct dirent *ent;
    string suffix = ".read-cov.txt";
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
    cerr << "Found " << txt_files.size() << " Read Coverage Files" << endl;

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
            string sid, tmp;
            float cov;
            lin >> sid >> tmp >> cov;
            if(running_cov.count(compress_string(sid)) == 0){
            	running_cov.insert(make_pair(compress_string(sid), 0.0));
            }
            running_cov[compress_string(sid)] += cov;
        }
    }

    ofstream myfile(argv[2]);

    for (map<string, float>::iterator it = running_cov.begin(); it != running_cov.end(); ++it)
    {
    	if (myfile.is_open())
		{
		myfile <<  it->first << "\t0\t" << it->second << "\n";
		}
    }
    myfile.close();

    return 0;
}