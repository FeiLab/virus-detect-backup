rm(list=ls());#����ڴ棬���¿�ʼ����
sub_routines<-paste(getwd(),"/bin/sRNA_clean.R",sep = "");
source(sub_routines);
library(ShortRead)
#ȥ�����С�N�������15bp������
#����Ĭ��ֵ
n_Cutoff=1;#ȥ������1����N����reads
read_Length=15; #ȥ������15bp��reads
Read_PerYield=5e5;#5Mreads*4*100=2G�ֽ�

#���봦�������в���
cmd_args = commandArgs(trailingOnly=TRUE);#��ȡ�ļ�����ȫ������
help_doc <- "
Usage: Rscript run_sRNA_clean.R filelist=<FILE> nCutoff=<INT> readLength=<INT> RdPerYield=<INT>
Required(1):
	filelist		The name of a txt file containing a list of input file names without any suffix
Options(7):
	nCutoff			reads with N number >= nCutoff after 5' and 3' end cleaned will be removed
	readLength		reads with length < readLength after trimming will be removed
	RdPerYield		how many reads will be processed at one time to control the memory usage
"; 
for (arg in cmd_args) {
	#cat("  ",arg, "\n", sep="");#������, �鿴��ȡ�������в���
	if ( grepl("^h(elp)?$", arg, ignore.case=TRUE, perl = TRUE, fixed = FALSE, useBytes = FALSE) ) {
		cat(help_doc); 
		stop("Stop for help.\n"); 
	} else if ( grepl("^filelist=", arg, ignore.case=TRUE, perl = TRUE, fixed = FALSE, useBytes = FALSE) ) {
		file_list <- unlist(strsplit(arg, "=", fixed=TRUE))[2];   
	} else if ( grepl("^nCutoff=", arg, ignore.case=TRUE, perl = TRUE, fixed = FALSE, useBytes = FALSE) ) {
		n_Cutoff <- as.numeric(unlist(strsplit(arg, "=", fixed=TRUE))[2]);  #argĬ����character����  
	} else if ( grepl("^readLength=", arg, ignore.case=TRUE, perl = TRUE, fixed = FALSE, useBytes = FALSE) ) {
		read_Length <- as.numeric(unlist(strsplit(arg, "=", fixed=TRUE))[2]);  #argĬ����character���� 
	} else if ( grepl("^RdPerYield=", arg, ignore.case=TRUE, perl = TRUE, fixed = FALSE, useBytes = FALSE) ) {
		Read_PerYield <- as.numeric(unlist(strsplit(arg, "=", fixed=TRUE))[2]);  #argĬ����character����  
	}
}

fastqfiles <- read.table(file_list)#�������е������ļ�����
sample_names <- as.matrix(fastqfiles)[,1]#��1�����ļ���������ֻ��1��
inputfiles <- paste(sample_names,".trimmed",sep="")
outputfiles <- paste(sample_names,".clean",sep="")
title <- c("sample_file","Trimmed_reads","Trimmed_length","Cleaned_reads","Cleaned_length")
write(title,file = "cleaned.report", ncolumns =5,append = T, sep = "\t")
for(i in 1:length(inputfiles))
{
 cleanRead(fastqfile=inputfiles[i], cleaned_file=outputfiles[i], nCutoff=n_Cutoff, readLength=read_Length, RdPerYield=Read_PerYield);
 #system(paste("rm", inputfiles[i]))#����clean��ɾ��ԭʼ�����ļ����������nuhup����������ֹ�ɾ��
}