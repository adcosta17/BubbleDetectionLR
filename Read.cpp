#include "Read.hpp"


/*
Read::Read(){
		id = "";
		length = 0;
		domain = "";
		kingdom = "";
		phylum = "";
		order = "";
		family = "";
		clas = "";
		genus = "";
		species = "";
}

*/

Read::Read(std::string i,
		int l)
	: id(i), length(l){
		
		domain = "";
		kingdom = "";
		phylum = "";
		order = "";
		family = "";
		clas = "";
		genus = "";
		species = "";
	}

void Read::setTaxonomy(std::vector<std::string>& tmp){
	for (int i = 0; i < tmp.size(); ++i)
	{
		if(tmp[i].length() > 3){
			char tax = tmp[i].at(0);
			std::string val = tmp[i].substr(3);
			switch(tax){
				case 's': species = val; break;
				case 'g': genus = val; break;
				case 'f': family = val; break;
				case 'o': order = val; break;
				case 'c': clas = val; break;
				case 'p': phylum = val; break;
				case 'k': kingdom = val; break;
				case 'd': domain = val; break;
				default: break;
			}
		}
	}	
}

bool Read::parentLevel(std::string level){
	//Check each level down to the requested comparison level, or until classification runs out
	//Function only called if exact match on requested level not found
	if(level == ""){
		return true;
	}
	if(domain != "" && level.find(domain) == std::string::npos){
		return false;
	}
	if(kingdom != "" && level.find(kingdom) == std::string::npos){
		return false;
	}
	if(phylum != "" && level.find(phylum) == std::string::npos){
		return false;
	}
	if(clas != "" && level.find(clas) == std::string::npos){
		return false;
	}
	if(order != "" && level.find(order) == std::string::npos){
		return false;
	}
	if(family != "" && level.find(family) == std::string::npos){
		return false;
	}
	if(genus != "" && level.find(genus) == std::string::npos){
		return false;
	}
	return true;
}



// Returns the read's full classification string from desired level up
// If read doesn't have classification to that level it returns an empty string
std::string Read::getClassification(char level){

	switch(level){
		case 's': {if(species != "") return getClassification('g') + " | " + species;} break;
		case 'g': {if(genus != "") return getClassification('f') + " | " + genus;} break;
		case 'f': {if(family != "") return getClassification('o') + " | " + family;} break;
		case 'o': {if(order != "") return getClassification('c') + " | " + order;} break;
		case 'c': {if(clas != "") return getClassification('p') + " | " + clas;} break;
		case 'p': {if(phylum != "") return getClassification('k') + " | " + phylum;} break;
		case 'k': {return getClassification('d') + " | " + kingdom;} break; // Kingdom can be empty for Bacteria
		case 'd': return domain;
		default: return "";
	}

	return "";
}

std::string Read::getClassificationForFileName(char level){

	switch(level){
		case 's': {if(species != "") return getClassificationForFileName('g') + "_" + species;} break;
		case 'g': {if(genus != "") return getClassificationForFileName('f') + "_" + genus;} break;
		case 'f': {if(family != "") return getClassificationForFileName('o') + "_" + family;} break;
		case 'o': {if(order != "") return getClassificationForFileName('c') + "_" + order;} break;
		case 'c': {if(clas != "") return getClassificationForFileName('p') + "_" + clas;} break;
		case 'p': {if(phylum != "") return getClassificationForFileName('k') + "_" + phylum;} break;
		case 'k': {return getClassificationForFileName('d') + "_" + kingdom;} break; // Kingdom can be empty for Bacteria
		case 'd': return domain;
		default: return "";
	}

	return "";
}