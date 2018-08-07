all: BubbleDetect PafSubset

BubbleDetect: BubbleDetect.cpp Match.hpp Match.cpp MatchUtils.cpp MatchUtils.hpp Read.cpp Read.hpp
	g++ -g -std=c++11 -lboost_iostreams BubbleDetect.cpp Match.cpp MatchUtils.cpp Read.cpp -o BubbleDetect

PafSubset: PafSubset.cpp Match.hpp Match.cpp MatchUtils.cpp MatchUtils.hpp Read.cpp Read.hpp
	g++ -g -std=c++11 -lboost_iostreams PafSubset.cpp Match.cpp MatchUtils.cpp Read.cpp -o PafSubset

clean:
	rm BubbleDetect PafSubset
