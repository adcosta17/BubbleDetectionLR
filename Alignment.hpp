#include <string>

class Alignment 
{
public:
	std::string id;
	int length;
	int start;
	int end;
	int qual;

	Alignment(std::string i,
		int l, int s, int e, int q);

};