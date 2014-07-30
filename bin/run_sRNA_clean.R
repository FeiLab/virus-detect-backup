rm(list=ls());#清空内存，重新开始计算
sub_routines<-paste(getwd(),"/bin/sRNA_clean.R",sep = "");
source(sub_routines);
library(ShortRead)
#去除含有“N”或短于15bp的序列
#变量默认值
n_Cutoff=1;#去除含有1个“N”的reads
read_Length=15; #去除短于15bp的reads
Read_PerYield=5e5;#5Mreads*4*100=2G字节

#读入处理命令行参数
cmd_args = commandArgs(trailingOnly=TRUE);#读取文件名后全部参数
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
	#cat("  ",arg, "\n", sep="");#调试用, 查看读取的命令行参数
	if ( grepl("^h(elp)?$", arg, ignore.case=TRUE, perl = TRUE, fixed = FALSE, useBytes = FALSE) ) {
		cat(help_doc); 
		stop("Stop for help.\n"); 
	} else if ( grepl("^filelist=", arg, ignore.case=TRUE, perl = TRUE, fixed = FALSE, useBytes = FALSE) ) {
		file_list <- unlist(strsplit(arg, "=", fixed=TRUE))[2];   
	} else if ( grepl("^nCutoff=", arg, ignore.case=TRUE, perl = TRUE, fixed = FALSE, useBytes = FALSE) ) {
		n_Cutoff <- as.numeric(unlist(strsplit(arg, "=", fixed=TRUE))[2]);  #arg默认是character类型  
	} else if ( grepl("^readLength=", arg, ignore.case=TRUE, perl = TRUE, fixed = FALSE, useBytes = FALSE) ) {
		read_Length <- as.numeric(unlist(strsplit(arg, "=", fixed=TRUE))[2]);  #arg默认是character类型 
	} else if ( grepl("^RdPerYield=", arg, ignore.case=TRUE, perl = TRUE, fixed = FALSE, useBytes = FALSE) ) {
		Read_PerYield <- as.numeric(unlist(strsplit(arg, "=", fixed=TRUE))[2]);  #arg默认是character类型  
	}
}

fastqfiles <- read.table(file_list)#读入所有的输入文件名称
sample_names <- as.matrix(fastqfiles)[,1]#第1列是文件名，可以只有1列
inputfiles <- paste(sample_names,".trimmed",sep="")
outputfiles <- paste(sample_names,".clean",sep="")
title <- c("sample_file","Trimmed_reads","Trimmed_length","Cleaned_reads","Cleaned_length")
write(title,file = "cleaned.report", ncolumns =5,append = T, sep = "\t")
for(i in 1:length(inputfiles))
{
 cleanRead(fastqfile=inputfiles[i], cleaned_file=outputfiles[i], nCutoff=n_Cutoff, readLength=read_Length, RdPerYield=Read_PerYield);
 #system(paste("rm", inputfiles[i]))#数据clean后，删除原始数据文件，如果不是nuhup，可以最后手工删除
}
