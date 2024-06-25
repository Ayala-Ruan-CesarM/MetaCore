#!/bin/bash 
#Codigo desarrollado para la ejecución y obtención de MetaCore
# Linea 4 a 59 permite ejecutar el codigo desde una terminal compatible con linux
# Tomando argumentos de entrada
#                                                                                                           #
  if [ "$1" = "-?" ] || [ "$1" = "-h" ] || [ $# -le 1 ]                                                     #
  then                                                                                                      #
    cat << EOF                                                       


 USAGE:  Extract_Core.sh  -i <list>  -o <basename>  [options]

 OPTIONS:
 
    -i <input>      Input list of metagenomes one in each line(mandatory)
    -o <output>     Output name for metagenome core (mandatory)
    -t <int>        Number of threads to use (Default: 10)
    -s <int>        Size of the segment length to align (Default: 1000)
    -p <int>        Percentage of identity to used in the aligments (Default: 85)
    -n <input>      Input list of negative metagenomes one in each line(mandatory)
    -k <int>        Kmer size (Default: 19)
    -f  <string>    Results folder name (Default: Results) 

EOF
    if [ "$1" = "-?" ] || [ "$1" = "-h" ] || [ $# -le 0 ]; then exit 0 ; fi                                 #
    exit 1 ;                                                                                                #
  fi                                                                                                        #
#                                                                                                           #

export LC_ALL=C ;

INPUT_LIST="N.O_L.I.S.T";       # -i (mandatory)
OUTPUT_FILE="N.O_O.U.T";        # -o (mandatory)
NPROC=10;                       # -t (10)
SEG_LENGTH=1000;                # -s (1000)
PERCENT_IDENT=85;               # -p (85)
INPUT_NEG_LIST="N.O_L.I.S.T";    # -n (mandatory)
KMER=19;                        # -k 19
FOLDER="Results";               # -f (Results)

while getopts :i:o:t:s:p:n:k::f:nfx option
do
  case $option in
    i) INPUT_LIST="$OPTARG"                                                            ;;
    o) OUTPUT_FILE="$OPTARG"                                                           ;;
    t) NPROC=$OPTARG                                                                 ;;
    s) SEG_LENGTH=$OPTARG                                                               ;;
    p) PERCENT_IDENT=$OPTARG                                                               ;;
    n) INPUT_NEG_LIST="$OPTARG"                                                     ;;
    k) KMER=$OPTARG                                                                  ;;  
    f) FOLDER=$OPTARG                                                                 ;;
    :) echo "missing argument ($OPTARG)" >&2 ; exit 1                               ;;
   \?) echo "invalid option ($OPTARG)"   >&2 ; exit 1                               ;;
  esac
done

if [ "$INPUT_LIST" == "N.O_L.I.S.T" ]; then echo "genome directory is not specified (option -i)"       >&2 ; exit 1 ; fi
if [ "$INPUT_NEG_LIST" == "N.O_L.I.S.T" ]; then echo "Negative metagenome list is not specified (option -n)"       >&2 ; exit 1 ; fi
if [ "$OUTPUT_FILE" == "N.O_O.U.T" ];  then echo "basename is not specified (option -o)"               >&2 ; exit 1 ; fi

# The menu is copied and modified to my needs from original Criscuolo A (2020), JolyTree bash script 
# available at https://gitlab.pasteur.fr/GIPhy/JolyTree
# Inicia variable que almacena la lista de metagenomas positivos
metagenome_list=($(cat $INPUT_LIST))
# 
mkdir -p Int_Files_S1
mkdir -p Int_Files_S2
mkdir -p "$FOLDER"    
ulimit -S -s unlimited

# Loop through the metagenomes in the list
for ((i=0; i<${#metagenome_list[@]}; i++)); do

    # Determine the first and second metagenomes
    if [ $i -eq 0 ]; then
        first_metagenome=${metagenome_list[i]}
        second_metagenome=${metagenome_list[i+1]}
        ref_metagenome=Ref_with_${second_metagenome%.fasta}.fasta
        output_prefix=MAP_1
        output_tmp_RN=Rename.fasta
    elif [[ $i -eq ${#metagenome_list[@]}-1 ]]; then
        break
    else
        first_metagenome=$ref_metagenome
        second_metagenome=${metagenome_list[i+1]}
        ref_metagenome=Ref_with_${second_metagenome%.fasta}.fasta
        output_prefix=MAP_$((i+1))
    fi

    file_size_bytes=$(stat -c %s "$first_metagenome")

    echo "Step 1: Running mashmap"

    if [ "$file_size_bytes" -lt 2600000000 ]; then
      echo "Este archivo pesa menos de 2.6GB"
      mashmap -r "$first_metagenome" -q "$second_metagenome" -s "$SEG_LENGTH" -t "$NPROC" -k "$KMER" -o "$output_prefix.txt" --pi "$PERCENT_IDENT"
      
    else
      echo "este archivo pesa más de 2.6G. Dividiendo el archivo "
      out=tmp_split
      folder=tmpsplit
      # This seqkit split2 line requires extensive testing as is the solution to RAM constraints
      # Maybe make it as a funcion of an user input RAM parameter. To Be implemented 25/06/24
      seqkit split2 -p 2 "$first_metagenome" -O "$folder" -o "$out"
      mashmap -r "$folder/$out.part_001.fasta" -q "$second_metagenome" -s "$SEG_LENGTH" -t "$NPROC" -k "$KMER" -o part1.split --pi "$PERCENT_IDENT"
      mashmap -r "$folder/$out.part_002.fasta" -q "$second_metagenome" -s "$SEG_LENGTH" -t "$NPROC" -k "$KMER" -o part2.split --pi "$PERCENT_IDENT"
      cat *.split > "$output_prefix.txt"
      if [ -d "$folder" ]; then
          rm -r "$folder"
        fi
        # Remove the split files
      rm -f *.split
    fi 
    sleep 5
    echo "Ending mashmap"
    echo "Step 2: Extracting coordinates"
    # Extract the coordinates
    awk '{print $6":"$8"-"$9}' $output_prefix.txt > ST_I$i.txt
    # Sort and overwrite ST_I$
    sort -u ST_I$i.txt -o sST_I$i.txt
	  sleep 5	
	  python remove_redundant_Cords_ST.py sST_I$i.txt rST_I$i.txt
    sed -i 's/:0/:1/g' rST_I$i.txt #change
    echo "Step 3: Running Samtools"
	  xargs samtools faidx $first_metagenome < rST_I$i.txt > MAP$i.fasta
	  sleep 5	
    echo "Sequences saved to "MAP$i.fasta" " ## new	
	  echo "Step 4: Creating New Reference"
	  # Create the reference for the next cycle
    cat $first_metagenome $second_metagenome > $ref_metagenome
    sleep 5
    # Rename all sequences in ref_metagenome
    awk '/>.*/{sub(/[^;]*/,">Sequence_" ++i )}1' $ref_metagenome > $output_tmp_RN 
    sleep 10
    mv $output_tmp_RN $ref_metagenome
    sleep 5	  
    mv *.txt Int_Files_S1/

done 

sleep 30

cat MAP* > "$OUTPUT_FILE.fasta"

sleep 30

awk '/>.*/{sub(/[^;]*/,">Sequence_" ++i )}1' $OUTPUT_FILE.fasta > $output_tmp_RN 
mv $output_tmp_RN $OUTPUT_FILE.fasta

echo "#################################### Creation of Metagenome Core complete #########################################"
echo "#################################### Removing intermidiate references #############################################"
rm -f Ref_with_*
rm -f *.fai
mv *MAP* "$FOLDER/"
echo "#################################### Individual Positive MAPS saved to "$FOLDER" folder ##############################"
echo "#################################### Extracting sequences share by negative metagenomes ############################"
# Read the input file
input_core=$OUTPUT_FILE.fasta
metagenome_neg_list=($(cat $INPUT_NEG_LIST))
seqkitfolder=tmp_split
# unlimited 
ulimit -S -s unlimited
# Loop through the metagenomes in the list
for ((i=0; i<${#metagenome_neg_list[@]}; i++)); do
	output_prefix=MAPn_$((i+1))
	second_metagenome=${metagenome_neg_list[i]}
  echo "Step 1: Running mashmap"
  ####Aqui este file_size_bytes esta declarado de acuerdo a la primer parte, positivos. Cambiar para que sea sobre el CORE 
  file_size_core=$(stat -c %s "$input_core")
  if [ "$file_size_core" -lt 2600000000 ]; then
    echo "Este archivo pesa menos de 2.6GB"
    mashmap -r "$input_core" -q "$second_metagenome" -s "$SEG_LENGTH" -t "$NPROC" -k "$KMER" -o "$output_prefix.txt" --pi "$PERCENT_IDENT"
  else
    echo "este archivo pesa más de 2.6G. Dividiendo el archivo "   
    if [ -d "$seqkitfolder" ]; then
      echo "The "$seqkitfolder" directory already exists. Skipping seqkit."
    else 
      seqkit split2 -p 2 "$input_core" -O "$seqkitfolder"
      echo "Splitting Core reference sequences completed"
    fi
    mashmap -r "$seqkitfolder"/"$OUTPUT_FILE.part_001.fasta" -q "$second_metagenome" -s "$SEG_LENGTH" -t "$NPROC" -k "$KMER" -o part1.split --pi "$PERCENT_IDENT"
    mashmap -r "$seqkitfolder"/"$OUTPUT_FILE.part_002.fasta" -q "$second_metagenome" -s "$SEG_LENGTH" -t "$NPROC" -k "$KMER" -o part2.split --pi "$PERCENT_IDENT"
    cat *.split > "$output_prefix.txt"
  fi
  sleep 5
  echo "Ending mashmap"
  echo "Step 2: Extracting coordinates"
  sleep 5
  awk '{print $6":"$8"-"$9}' $output_prefix.txt > ST_I$i.txt
  sort -u ST_I$i.txt -o sST_I$i.txt
  sleep 5       
  python remove_redundant_Cords_ST.py sST_I$i.txt rST_I$i.txt
  sed -i 's/:0/:1/g' rST_I$i.txt
  sleep 5
  echo "Step 3: Running Samtools"
  xargs samtools faidx $input_core < rST_I$i.txt > MAPNeg$i.fasta
  sleep 5
  mv *.txt Int_Files_S2/
  sleep 5
  rm -f *.split
done

mv *MAPNeg* "$FOLDER/"
rm -f *.fai
rm -r "$seqkitfolder"
echo "######################################### Individual Negative MAPS saved to folder "$FOLDER" #########################################"
