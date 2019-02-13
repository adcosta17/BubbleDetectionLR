all: BubbleDetect PafSubset ReadClassificationGen SplitFastq CombinePaf ComputeContained ProcessBubbleList ReadCoverageCombined FilterContained

BubbleDetect: BubbleDetect.cpp Match.hpp Match.cpp MatchUtils.cpp MatchUtils.hpp Read.cpp Read.hpp
	g++ -g -std=c++11 -lboost_iostreams BubbleDetect.cpp Match.cpp MatchUtils.cpp Read.cpp -o BubbleDetect

PafSubset: PafSubset.cpp Match.hpp Match.cpp MatchUtils.cpp MatchUtils.hpp Read.cpp Read.hpp
	g++ -g -std=c++11 -lboost_iostreams PafSubset.cpp Match.cpp MatchUtils.cpp Read.cpp -o PafSubset

ReadClassificationGen: ReadClassificationGen.cpp
	g++ -g -std=c++11 ReadClassificationGen.cpp -o ReadClassificationGen

SplitFastq: SplitFastq.cpp Read.cpp
	g++ -g -std=c++11 -lboost_iostreams SplitFastq.cpp Read.cpp -lm -lz -o SplitFastq

CombinePaf: CombinePaf.cpp
	g++ -g -std=c++11 -lboost_iostreams CombinePaf.cpp -lm -lz -o CombinePaf

ComputeContained: ComputeContained.cpp
	g++ -g -std=c++11 -lboost_iostreams ComputeContained.cpp Match.cpp -lm -lz -o ComputeContained 

ProcessBubbleList: ProcessBubbleList.cpp Read_Alignment.cpp Read_Alignment.hpp Alignment.cpp Alignment.hpp
	g++ -g -std=c++11 ProcessBubbleList.cpp Read_Alignment.cpp Alignment.cpp -o ProcessBubbleList

ReadCoverageCombined: ReadCoverageCombined.cpp
	g++ -g -std=c++11 -lboost_iostreams ReadCoverageCombined.cpp -lm -lz -o ReadCoverageCombined

FilterContained: FilterContained.cpp
	g++ -g -std=c++11 -lboost_iostreams FilterContained.cpp -lm -lz -o FilterContained

clean:
	rm BubbleDetect PafSubset ReadClassificationGen SplitFastq CombinePaf ComputeContained ProcessBubbleList ReadCoverageCombined FilterContained
