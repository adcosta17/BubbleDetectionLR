all: BubbleDetect PafSubset ReadClassificationGen SplitFastq CombinePaf ComputeContained ProcessBubbleList ReadCoverageCombined FilterContained ProcessBubbleListModel

BubbleDetect: BubbleDetect.cpp Match.hpp Match.cpp MatchUtils.cpp MatchUtils.hpp Read.cpp Read.hpp
	g++ -g -std=c++11 -Wall -lboost_iostreams BubbleDetect.cpp Match.cpp MatchUtils.cpp Read.cpp -o BubbleDetect

PafSubset: PafSubset.cpp Match.hpp Match.cpp MatchUtils.cpp MatchUtils.hpp Read.cpp Read.hpp
	g++ -g -std=c++11 -Wall -lboost_iostreams PafSubset.cpp Match.cpp MatchUtils.cpp Read.cpp -o PafSubset

ReadClassificationGen: ReadClassificationGen.cpp
	g++ -g -std=c++11 -Wall ReadClassificationGen.cpp -o ReadClassificationGen

SplitFastq: SplitFastq.cpp Read.cpp
	g++ -g -std=c++11 -Wall -lboost_iostreams SplitFastq.cpp Read.cpp -lm -lz -o SplitFastq

CombinePaf: CombinePaf.cpp
	g++ -g -std=c++11 -Wall -lboost_iostreams CombinePaf.cpp -lm -lz -o CombinePaf

ComputeContained: ComputeContained.cpp
	g++ -g -std=c++11 -Wall -lboost_iostreams ComputeContained.cpp Match.cpp -lm -lz -o ComputeContained 

ProcessBubbleList: ProcessBubbleList.cpp Read_Alignment.cpp Read_Alignment.hpp Alignment.cpp Alignment.hpp
	g++ -g -std=c++11 -Wall ProcessBubbleList.cpp Read_Alignment.cpp Alignment.cpp -o ProcessBubbleList

ProcessBubbleListModel: ProcessBubbleListModel.cpp Read_Alignment.cpp Read_Alignment.hpp Alignment.cpp Alignment.hpp
	g++ -g -std=c++11 -Wall ProcessBubbleListModel.cpp Read_Alignment.cpp Alignment.cpp -o ProcessBubbleListModel

ReadCoverageCombined: ReadCoverageCombined.cpp
	g++ -g -std=c++11 -Wall -lboost_iostreams ReadCoverageCombined.cpp -lm -lz -o ReadCoverageCombined

FilterContained: FilterContained.cpp
	g++ -g -std=c++11 -Wall -lboost_iostreams FilterContained.cpp -lm -lz -o FilterContained

clean:
	rm BubbleDetect PafSubset ReadClassificationGen SplitFastq CombinePaf ComputeContained ProcessBubbleList ReadCoverageCombined FilterContained ProcessBubbleListModel
