#include <string>
#include <vector>
#include "Alignment.hpp"

class Read_Alignment 
{
public:
	std::string id;
	int length;
	bool split_1;
	bool split_2;
	std::vector<Alignment> species1_aln;
	std::vector<Alignment> species2_aln;

	Read_Alignment(std::string i,
		int l);

	bool align_species_1();
	bool align_species_2();

};