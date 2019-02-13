#include "Read_Alignment.hpp"

Read_Alignment::Read_Alignment(std::string i,
		int l)
	: id(i), length(l){
		split_1 = false;
		split_2 = false;
	}


bool Read_Alignment::align_species_1(){
	if(split_1){
		// Read is split between multiple alignments
		// IE Chimeric or doesn't map fully
		return false;
	}
	int base_size = 0;
	for(int i = 0; i < species1_aln.size(); i++){
		if(species1_aln[i].qual > 30){
			base_size += (species1_aln[i].end - species1_aln[i].start);
		}
	}
	if(base_size/static_cast<float>(length) > 0.90){
		return true;
	}
	else{
		return false;
	}
}

bool Read_Alignment::align_species_2(){
	if(split_1){
		// Read is split between multiple alignments
		// IE Chimeric or doesn't map fully
		return false;
	}
	int base_size = 0;
	for(int i = 0; i < species2_aln.size(); i++){
		if(species2_aln[i].qual > 30){
			base_size += (species2_aln[i].end - species2_aln[i].start);
		}
	}
	if(base_size/static_cast<float>(length) > 0.90){
		return true;
	}
	else{
		return false;
	}
}