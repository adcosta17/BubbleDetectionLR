#include "Match.hpp"
#include <algorithm>
#include <iostream>
#include <fstream>

Match::Match(std::string qrid, int qrl, int qrs, int qre,
	char str, std::string trid, int trl, int trs, int tre, int rm,
	int ab, int pl, int sl, std::string cg, int threshold, bool rev)
	: query_read_id(qrid), query_read_length(qrl),
	query_read_start(qrs), query_read_end(qre),
	strand(str), target_read_id(trid), target_read_length(trl),
	target_read_start(trs), target_read_end(tre),
	residue_matches(rm), alignment_block(ab), 
	reversed(rev), prefix_length(pl), suffix_length(sl),
	cigar(cg), reduce(false), length_to_use(0), overlap_species("") {

		if(strand == '+'){
			if(query_read_start > target_read_start){
				length = (target_read_length - target_read_end) - (query_read_length - query_read_end);
				suffix_length = length;
				orientation = 1;
			} else {
				length = target_read_start - query_read_start;
				suffix_length = length;
				orientation = -1;
			}
		} else {
			if(query_read_start > query_read_length - query_read_end && target_read_start > target_read_length - target_read_end){
				length = target_read_start - (query_read_length - query_read_end);
				suffix_length = length;
				orientation = 1;
			} else {
				length = (target_read_length - target_read_end) - query_read_start;
				suffix_length = length;
				orientation = -1;
			}
		}

		if(strand == '+'){
			if(query_read_start < target_read_start){
				prefix_length = (query_read_length - query_read_end) - (target_read_length - target_read_end);
				//orientation = 1;
			} else {
				prefix_length = query_read_start - target_read_start;
				//orientation = -1;
			}
		} else {
			if(query_read_start > query_read_length - query_read_end && target_read_start > target_read_length - target_read_end){
				prefix_length = query_read_start - (target_read_length - target_read_end);
				//orientation = 1;
			} else {
				prefix_length = (query_read_length - query_read_end) - target_read_start;
				//orientation = -1;
			}
		}
	}

bool Match::internal_edge(){
	int tl5, tl3, ext5, ext3;
	int max_hang = 1000;
	if (strand == '-'){
		tl5 = target_read_length - target_read_end;
		tl3 = target_read_start; // tl5: 5'-end overhang (on the query strand); tl3: similar
	} else {
		tl5 = target_read_start;
		tl3 = target_read_length - target_read_end;
	}
	ext5 = query_read_start < tl5? query_read_start : tl5;
	ext3 = query_read_length - query_read_end < tl3? query_read_length - query_read_end : tl3;
	if (ext5 > max_hang || ext3 > max_hang || query_read_end - query_read_start < (query_read_end - query_read_start + ext5 + ext3) * 0.8) {
		return true;
	} else {
		return false;
	}
}

int Match::getQueryAlignmentLen()
{
	return query_read_end - query_read_start;
}

int Match::getTargetAlignmentLen()
{
	return target_read_end - target_read_start;
}

void Match::sort_matches(std::vector<Match>& matches)
{
	std::sort(matches.begin(), matches.end(), [](const Match& lhs, const Match& rhs){
		return lhs.length < rhs.length;
	});
	
	//int num = std::min(5, static_cast<int>(matches.size()));
	//matches.erase(matches.begin()+num, matches.end());
}

int Match::check_match_contained()
{
	int tl5, tl3, ext5, ext3;
	int max_hang = 1000;
	if (strand == '-'){
		tl5 = target_read_length - target_read_end;
		tl3 = target_read_start; // tl5: 5'-end overhang (on the query strand); tl3: similar
	} else {
		tl5 = target_read_start;
		tl3 = target_read_length - target_read_end;
	}
	ext5 = query_read_start < tl5 ? query_read_start : tl5;
	ext3 = query_read_length - query_read_end < tl3 ? query_read_length - query_read_end : tl3;
	if (ext5 > max_hang || ext3 > max_hang || query_read_end - query_read_start < (query_read_end - query_read_start + ext5 + ext3) * 0.8){
		return 0;
	}
	if (query_read_start <= tl5 && query_read_length - query_read_end <= tl3) return -1; // query contained
	else if (query_read_start >= tl5 && query_read_length - query_read_end >= tl3) return 1; // target contained
	if (query_read_end - query_read_start + ext5 + ext3 < 2000 || target_read_end - target_read_start + ext5 + ext3 < 2000) return 2; // short overlap
	//length = l;

	// Filter out reads that have a very high overlap percentage. Minimap2 maybe missing some portions at start or end. Anything with 90% overlap as a fraction fo the read is going to be considered contained
	// To be implemented as needed
	return 0;
}

