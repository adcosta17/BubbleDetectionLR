# BubbleDetectionLR

Genome assembly and Bubble Detection for Nanopore long reads

BubbleDetect is designed to take in a set of all v all overlaps, lay them out in an assembly graph and compute the bubbles that occur within them. Each bubble is scored based on a number of characterisitcs including taxinomic classification and coverage if provided by the user. The assembly generated is outputted along with the bubble list to a user specified location.  


## Usage
```sh
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
### Required
**-p** The main input to BubbleDetect. Either a paf file (can be gz compressed) or a directory of paf files. If a directory assumes reads have been split into subsets and an all vs all overlap done between sets.

**-o** The output path and prefix. All output files including the assembly stats, gfa files and bubble lists will be written here with the given prefix

### Optional
**-i** The number of times dead end and tip removal is done

**-f** The fuzz value used in Myers transitive edge reduction

**-t** The multiple of the mean read length to use as a minimum size needed to retain a dead end or tip

**-r** Read classification file, tab delimited. Reads as first column, taxinomic classification as second. Levels of classification separated by ;

**-s** Species to colour map. Provides a RGB hex clode for each species of interest for colouring the output GFA

**-h** Chimeric read map. Tab deliminted. Reads as first coloumn, Call as chimeric or not in second

**-c** Coverage per read. Tab deliminited. Reads as first column, Lengths as second. Coverage as third

**-m** Secondary classification file. Used for when reads are binned by a taxinomic group. Each read is classified to a set taxinomic level if possible. Reads as first column, taxinomic classification as second. Levels of classification separated by |

**-g** Total estimated genome size, when passed in NG50 is used instead of N50

**-x** Flag to skip bubble detection. Useful if only a quick assembly is needed

**-z** Use read names as is. Read names are assumed to be 32 hexchar nanopore read ids. These will normally be compressed. If flag is passed compression is skipped.

**-b** Denotes input directory passed in represents a set of paf files where each file represents overlaps within a species or other taxinomic level. Reads are assembled within each species and the graphs combined.

**-l** Compute GFA file of collapsed contigs. Normally GFA file output is at read level. Here a second GFA file is generated with read overlaps collapsed to contigs.

**-a** Flag to skip computing assembly stats such as N50. Done if only a read level GFA is needed. 
