# MetaCore
MetaCore: Unique core sequences from metagenomics samples for microbiome wide association studies

Both MetaCore.sh and remove_redundant_Cords.ST.py had to be in working folder

Positive.list : text file with the name of every positive metagenome to analyze (one name per line)

Negative.list : text file with the name of every negative metagenome to analyze (one name per line)


USAGE:  MetaCore.sh  -i Positive.list -n Negative.list -o <basename>  [options]

 OPTIONS:

    -i <input>      Input list of metagenomes one in each line(mandatory)
    -o <output>     Output name for metagenome core (mandatory)
    -t <int>        Number of threads to use (Default: 10)
    -s <int>        Size of the segment length to align (Default: 1000)
    -p <int>        Percentage of identity to used in the aligments (Default: 85)
    -n <input>      Input list of negative metagenomes one in each line(mandatory)
    -k <int>        Kmer size (Default: 19)
    -f  <string>    Results folder name (Default: Results)


