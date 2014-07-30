##############################
## Hierarchy of subroutines ##
##############################   
virus_detect.pl
to produce known contigs:
	bwa-alignAndCorrect.pl
		bwa
		SAM_filter_out_unmapped_reads.pl
		samtools
		samFilter.pl
		pileup_filter()
		extractConsensus.clas
		renameFasta.pl
	removeRedundancy_batch.pl
		dust
		trim_XNseq1.pl
		removeRedundancy.pl
			formatdb
			megablast
		bowtie
		samtools			
		extractConsensus.class
		renameFasta.pl
		
to produce novel contigs:		
	bwa_remove.pl	
		bwa
		SAM_filter_out_unmapped_reads.pl
 	Velvet_Optimiser.pl
	runVelvet
	removeRedundancy_batch.pl
	
to produce combined contigs:		
	files_combine2.pl	
	removeRedundancy_batch.pl

virus_detect.pl	
#################################################
## clean还是unmapped来correct组装好的contigs   ##
#################################################
./bin/bwa-alignAndCorrect.pl --file_list file_names.list --reference ./databases/virus_genbank186.fasta --coverage 0.3

################
## 文件清单   ##
################
virus_detect
	databases
		virus_genbank186.fasta
		sweet_potato.fasta
		rRNA.fasta
		
		
column1: The identified virus ID from the virus database
column2: The length of the identified virus
column3: The total positions covered on the identified virus
column4: The coverage of the identified virus(col3/col2) 
column5: The average depth of the identified virus 

##########################
##   搜索报错的关键词   ##
##########################
Error
Aborted
cannot
No
core dumped