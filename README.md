# BubbleDetectionLR

Genome assembly and Bubble Detection for Nanopore long reads

BubbleDetect is designed to take in a set of all v all overlaps, lay them out in an assembly graph and compute the bubbles that occur within them. Each bubble is scored based on a number of characterisitcs including taxinomic classification and coverage if provided by the user. The assembly generated is outputted along with the bubble list to a user specified location.  


## Setup, Usage & Examples
```sh
# Setup:
git clone https://github.com/adcosta17/BubbleDetectionLR.git
cd BubbleDetectionLR && make

# Usage: 
# ./BubbleDetect -p INPUT -o OUT_PATH [OPTIONS]

# Examples:

# Basic Usage:
./BubbleDetect -p all_v_all.paf -o output/path/and/prefix

# With Taxinomic Classification and Coverage
./BubbleDetect -p all_v_all.paf -o output/path/and/prefix -r read_classification.txt -s species_to_colour.txt -c coverage_per_read.txt

# Input is a directory of paf files
./BubbleDetect -p directory/of/files/ -o output/path/and/prefix -r read_classification.txt -s species_to_colour.txt -c coverage_per_read.txt

# Input is binned by species
./BubbleDetect -p directory/of/files/binned/ -o output/path/and/prefix -r read_classification.txt -s species_to_colour.txt -c coverage_per_read.txt -b -m secondary_classification_file.txt

# Skip Bubble Detection, Assembly stats and output collapsed contig GFA
./BubbleDetect -p all_v_all.paf -o output/path/and/prefix -l -a -x

# Change default fuzz, iteration and threshold values
./BubbleDetect -p all_v_all.paf -o output/path/and/prefix -i 2 -f 1500 -t 3

```
## Arguments:
### Required
**-p** The main input to BubbleDetect. Either a paf file (can be gz compressed) or a directory of paf files. If a directory assumes reads have been split into subsets and an all vs all overlap done between sets.

**-o** The output path and prefix. All output files including the assembly stats, gfa files and bubble lists will be written here with the given prefix

### Optional
**-i** The number of times dead end and tip removal is done [10]

**-f** The fuzz value used in Myers transitive edge reduction [1000]

**-t** The multiple of the mean read length to use as a minimum size needed to retain a dead end or tip [5]

**-r** Read classification file, tab delimited. Reads as first column, taxinomic classification as second. Levels of classification separated by ;

**-s** Species to colour map. Provides a RGB hex clode for each species of interest for colouring the output GFA

**-h** Chimeric read map. Tab deliminted. Reads as first coloumn, Call as chimeric or not in second

**-c** Coverage per read. Tab deliminited. Reads as first column, Lengths as second. Coverage as third

**-m** Secondary classification file. Used for when reads are binned by a taxinomic group. Each read is classified to a set taxinomic level if possible (if using kraken it is in MPA format). Reads as first column, taxinomic classification as second. Levels of classification separated by |

**-g** Total estimated genome size, when passed in NG50 is used instead of N50

**-x** Flag to skip bubble detection. Useful if only a quick assembly is needed [False]

**-z** Use read names as is. Read names are assumed to be 32 hexchar nanopore read ids. These will normally be compressed. If flag is passed compression is skipped. [False]

**-b** Denotes input directory passed in represents a set of paf files where each file represents overlaps within a species or other taxinomic level. Reads are assembled within each species and the graphs combined. [False]

**-l** Compute GFA file of collapsed contigs. Normally GFA file output is at read level. Here a second GFA file is generated with read overlaps collapsed to contigs. [False]

**-a** Flag to skip computing assembly stats such as N50. Done if only a read level GFA is needed. [False]
