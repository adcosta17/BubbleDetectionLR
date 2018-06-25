all: BubbleDetect PafSubset

BubbleDetect: BubbleDetect.cpp Match.hpp Match.cpp MatchUtils.cpp MatchUtils.hpp
	g++ -g -std=c++11 BubbleDetect.cpp Match.cpp MatchUtils.cpp -o BubbleDetect

PafSubset: PafSubset.cpp Match.hpp Match.cpp MatchUtils.cpp MatchUtils.hpp
	g++ -g -std=c++11 PafSubset.cpp Match.cpp MatchUtils.cpp -o PafSubset

clean:
	rm BubbleDetect PafSubset