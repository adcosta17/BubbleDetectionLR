#include <string>
#include <vector>

class Read 
{
public:
	std::string id;
	int length;
	std::string domain;
	std::string kingdom;
	std::string phylum;
	std::string order;
	std::string family;
	std::string clas;
	std::string genus;
	std::string species;
	std::string classifcation;

	Read(std::string i,
		int l);

	//Read();

	void setTaxonomy(std::vector<std::string>& tmp);
	std::string getClassification(char level);
	std::string getClassificationForFileName(char level);
	bool parentLevel(std::string level);
};