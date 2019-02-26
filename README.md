# BubbleDetectionLR

Genome assembly and Bubble Detection for Nanopore long reads

BubbleDetect is designed to take in a set of all v all overlaps, lay them out in an assembly graph and compute the bubbles that occur within them.


## Usage
```
BubbleDetect
 -p .paf Input File or directory containing paf files
 -o Output Path and Prefix
 (optional) -i Iterations# [10]
 (optional) -f fuzz [1000]
 (optional) -t threshold [5]
 (optional) -r read classification file
 (optional) -s species to colour map
 (optional) -h chimeric read map
 (optional) -c read coverage map
 (optional) -m kraken mpa classification file
 (optional) -g total estimated genome size (used for NG50 rather than N50 Calculation)
 (optional) -x skip bubble detection if needed [False]
 (optional) -z use_names_as_is, don't compress read ids. Useful if read ids aren't in normal Nanopre ID form [False]
 (optional) -b input is binned by species [False]
 (optional) -l compute a gfa file representing collapsed contigs [False]
 (optional) -a output assembly statistics [False]

```
## Arguments:
### -p
The main input to BubbleDetect. Either a paf file (can be gz compressed) or a directory of paf files. If a directory assumes reads have been split into subsets and an all vs all overlap done between sets.

### -o
The output path and prefix. All output files including the assembly stats, gfa files and bubble lists will be written here with the given prefix
