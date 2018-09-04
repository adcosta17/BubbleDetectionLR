all: BubbleDetect PafSubset ReadClassificationGen

BubbleDetect: BubbleDetect.cpp Match.hpp Match.cpp MatchUtils.cpp MatchUtils.hpp Read.cpp Read.hpp
	g++ -g -std=c++11 -lboost_iostreams BubbleDetect.cpp Match.cpp MatchUtils.cpp Read.cpp -o BubbleDetect

PafSubset: PafSubset.cpp Match.hpp Match.cpp MatchUtils.cpp MatchUtils.hpp Read.cpp Read.hpp
	g++ -g -std=c++11 -lboost_iostreams PafSubset.cpp Match.cpp MatchUtils.cpp Read.cpp -o PafSubset

ReadClassificationGen: ReadClassificationGen.cpp
	g++ -g -std=c++11 ReadClassificationGen.cpp -o ReadClassificationGen

clean:
	rm BubbleDetect PafSubset ReadClassificationGen
