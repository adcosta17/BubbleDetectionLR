#include <string>
#include <vector>
#include <map>
#include <set>

class Match 
{
public:
	std::string query_read_id;
	int query_read_length;
	int query_read_start;
	int query_read_end;
	char strand;
	std::string target_read_id;
	int target_read_length;
	int target_read_start;
	int target_read_end;
	int residue_matches;
	int alignment_block;
	bool reversed;
	int prefix_length;
	int suffix_length;
	std::vector<std::string> reads_contained;
	char fixed_strand;
	int length;
	int orientation;
	std::string cigar;
	bool reduce;
	int length_to_use;
	std::string overlap_species; 

	Match(std::string qrid,
		int qrl,
		int qrs,
		int qre,
		char str,
		std::string trid,
		int trl,
		int trs,
		int tre,
		int rm,
		int ab,
		int pl,
		int sl,
		std::string cg,
		int threshold = 100,
		bool rev = false);

	int getQueryAlignmentLen();
	int getTargetAlignmentLen();
	static void sort_matches(std::vector<Match>& matches);
	//std::pair<int, int> getEdgeLength(int previous_alignment);
	int check_match_contained();
	bool internal_edge();
};
