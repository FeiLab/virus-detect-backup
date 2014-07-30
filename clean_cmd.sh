#CTGTAGGCACCATCAAT
#ATCTCGTATGCCGTCTTCTGCTTG
#TGGAATTCTCGGGTGCC
#
#./tools/smallRNA_clean.pl --file_list list_fq --adapter TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC --adapterLen 11 --run_R --rRNA_removal --rRNA_reference rRNA.fasta

#./tools/smallRNA_clean.pl --file_list list_fq --adapter TCGTATGCCG --adapterLen 9 --run_R

#./tools/smallRNA_clean.pl --file_list list_fq --adapter CTGCTGGATCGT --adapterLen 11 --run_R 

./tools/smallRNA_clean.pl --file_list list_fq --adapter CTGTAGGCACCATCAAT --adapterLen 11 --run_R &

#--rRNA_removal --rRNA_reference rRNA_silva111.fasta
