#include "Match.hpp"
#include "Read.hpp"

class MatchUtils 
{
	public:
	static int read_paf_file(std::map<std::string, std::vector<Match> >& edge_lists, std::map<std::string, std::vector<Match> >& all_matches, std::map<std::string, std::vector<Match> >& raw_matces, std::set<std::string>& read_ids, std::map<std::string, int>& read_lengths, std::string file_name, std::set<std::string>& chimeric_reads, std::map<std::string, Read>& read_classification, bool gfa);
	static bool check_bubble(std::string start, std::string end, std::set<std::string> reads, std::map<std::string,std::vector<std::string> >& read_indegree, std::map<std::string,std::vector<std::string> >& read_outdegree);
	static void reduce_edges(std::map<std::string, std::vector<Match> >& all_matches, std::set<std::string>& read_ids, std::map<std::string, std::vector<Match> >& edge_lists, int fuzz);
	static void toGfa(std::map<std::string, std::vector<Match> >& all_matches, std::map<std::string, int>& read_lengths, std::string file_name, std::map<std::string, std::vector<std::string> >& read_indegree, std::map<std::string, std::vector<std::string> >& read_outdegree, std::map<std::string, std::string>& read_names, std::map<std::string, std::string>& colours);
	static void compute_dead_ends(std::map<std::string, std::vector<Match> >& all_matches, std::set<std::string>& read_ids, std::map<std::string, std::vector<std::string> >& read_indegree, std::map<std::string, std::vector<std::string> >& read_outdegree, std::set<std::string>& de_ids, std::map<std::string, std::vector<std::string> >& de_paths);
	static void compute_in_out_degree(std::map<std::string, std::vector<Match> >& all_matches, std::set<std::string>& read_ids, std::map<std::string, std::vector<std::string> >& read_indegree, std::map<std::string, std::vector<std::string> >& read_outdegree);
	static void prune_dead_paths(std::map<std::string, std::vector<Match> >& all_matches, std::set<std::string>& read_ids, std::map<std::string, std::vector<std::string> >& read_indegree, std::map<std::string, std::vector<std::string> >& read_outdegree, std::map<std::string, std::vector<std::string> >& de_paths, int mean_read_length, int threshold);
	static void compute_sets(std::map<std::string, std::vector<Match> >& all_matches, std::string start, std::string current, std::string end, std::map<std::pair<std::string,std::string>, std::set<std::string> >& bubble_sets, std::set<std::string> seen);
	static void find_bubble(std::string start, std::map<std::string,std::vector<std::string> >& read_indegree, std::map<std::string,std::vector<std::string> >& read_outdegree, std::map<std::pair<std::string,std::string>, std::set<std::string> >& bubble_sets, std::vector<std::string>& start_ids);
	static void filter_dead_end_branches(std::map<std::string, std::vector<Match> >& all_matches, std::set<std::string>& read_ids, std::map<std::string, int>& read_indegree, std::map<std::string, int>& read_outdegree, std::map<std::string, std::vector<std::string> >& de_paths);
	static void prune_branch_dead_ends(std::map<std::string, std::vector<Match> >& all_matches, std::map<std::string, int>& read_indegree, std::map<std::string, int>& read_outdegree,std::set<std::string>& read_ids);
	static void clean_matches(std::map<std::string, std::vector<Match> >& all_matches);
	static int mark_matches_for_node(std::map<std::string, std::vector<Match> >& all_matches, std::string id, std::map<std::string, int>& mark);
	static int compute_contigs(std::string id, std::map<std::string, std::vector<Match> >& all_matches, std::map<std::string,std::vector<std::string> >& read_indegree, std::map<std::string,std::vector<std::string> >& read_outdegree, std::map<int, std::vector<Match> >& contig_map, int contig_number);
	static void subset_matches(std::map<std::string, std::vector<Match> >& all_matches, std::map<std::string, std::vector<Match> >& edge_lists, std::map<std::string, std::vector<Match> >& species_matches, std::map<std::string, std::vector<Match> >&  species_edge_lists, std::set<std::string> ids_to_use);
};