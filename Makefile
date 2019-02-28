all: BubbleDetect PafSubset ReadClassificationGen SplitFastq CombinePaf ComputeContained ProcessBubbleList ReadCoverageCombined FilterContained ProcessBubbleListModel

BubbleDetect: BubbleDetect.cpp Match.hpp Match.cpp MatchUtils.cpp MatchUtils.hpp Read.cpp Read.hpp
	g++ -g -std=c++11 -Wall BubbleDetect.cpp Match.cpp MatchUtils.cpp Read.cpp -lboost_iostreams -lm -lz -o BubbleDetect

PafSubset: PafSubset.cpp Match.hpp Match.cpp MatchUtils.cpp MatchUtils.hpp Read.cpp Read.hpp
	g++ -g -std=c++11 -Wall PafSubset.cpp Match.cpp MatchUtils.cpp Read.cpp -lboost_iostreams -lm -lz -o PafSubset

ReadClassificationGen: ReadClassificationGen.cpp
	g++ -g -std=c++11 -Wall ReadClassificationGen.cpp -o ReadClassificationGen

SplitFastq: SplitFastq.cpp Read.cpp
	g++ -g -std=c++11 -Wall SplitFastq.cpp Read.cpp -lboost_iostreams -lm -lz -o SplitFastq

CombinePaf: CombinePaf.cpp
	g++ -g -std=c++11 -Wall CombinePaf.cpp -lboost_iostreams -lm -lz -o CombinePaf

ComputeContained: ComputeContained.cpp
	g++ -g -std=c++11 -Wall ComputeContained.cpp Match.cpp -lboost_iostreams -lm -lz -o ComputeContained 

ProcessBubbleList: ProcessBubbleList.cpp Read_Alignment.cpp Read_Alignment.hpp Alignment.cpp Alignment.hpp
	g++ -g -std=c++11 -Wall ProcessBubbleList.cpp Read_Alignment.cpp Alignment.cpp -o ProcessBubbleList

ProcessBubbleListModel: ProcessBubbleListModel.cpp Read_Alignment.cpp Read_Alignment.hpp Alignment.cpp Alignment.hpp
	g++ -g -std=c++11 -Wall ProcessBubbleListModel.cpp Read_Alignment.cpp Alignment.cpp -o ProcessBubbleListModel

ReadCoverageCombined: ReadCoverageCombined.cpp
	g++ -g -std=c++11 -Wall ReadCoverageCombined.cpp -lboost_iostreams -lm -lz -o ReadCoverageCombined

FilterContained: FilterContained.cpp
	g++ -g -std=c++11 -Wall FilterContained.cpp -lboost_iostreams -lm -lz -o FilterContained

clean:
	rm BubbleDetect PafSubset ReadClassificationGen SplitFastq CombinePaf ComputeContained ProcessBubbleList ReadCoverageCombined FilterContained ProcessBubbleListModel
