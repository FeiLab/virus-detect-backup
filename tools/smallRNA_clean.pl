#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;
my $usage = <<_EOUSAGE_;

#########################################################################################
# smallRNA_clean.pl --file_list <FILE>  --suffix [String] --adapter <String> --adapterLen [String]
# 		    --run_R --n_Cutoff [INT] --read_Length [INT] --read_PerYield [FLOAT] 
#                   --rRNA_removal --rRNA_reference [FILE] --virus_removal --virus_reference [FILE]
#                   --max_dist[INT] --max_open [INT] --max_extension [INT] --len_seed [INT] --dist_seed [INT] --thread_num [INT]
#                                  
# Required(2):
#  --file_list  The name of a txt file containing a list of input file names without any suffix 
#  --suffix     The suffix of input files
#
# Options(7):
#  --run_R            run R script to remove low quality, ambiguous and adapter bases [Not selected]
#  --rRNA_removal     Remove rRNA sequences [Not selected]
#  --rRNA_reference   a fasta file containing the rRNA reference sequences
#  --virus_removal    Remove virus sequences [Not selected]
#  --virus_reference  a fasta file containing the virus reference sequences
#  --adapter          the reverse complimentary of 3' end adapter(5'->3')
#  --adapterLen       the adapter offset from the 5' end [11]
#
# R script-related options(3): 
#  --n_Cutoff         after low quality bases trimmed, the reads containing "N" with a number above this parameter will be removed [1]  
#  --read_Length      reads with the length shorter than this parameter will be removed [15] 
#  --read_PerYield    Number of reads in each iteration to limit the memory usage [5e5] 
#
# BWA-related options(6):
#  --max_dist      Maximum edit distance [1]  
#  --max_open      Maximum number of gap opens [1]  
#  --max_extension Maximum number of gap extensions [1] 
#  --len_seed      Take the first INT subsequence as seed for the BWA alignment [15] 
#  --dist_seed     Maximum edit distance in the seed [1]  
#  --thread_num    Number of threads (multi-threading mode) [8] 
###########################################################################################

_EOUSAGE_

	;

#################
##   全局变量  ##
#################
our $file_list;
our $suffix="fastq";;          #输入文件（fastq格式）的后缀名
our $run_R;
our $rRNA_removal;
our $virus_removal;
our $rRNA_reference="rRNA.fasta";
our $virus_reference= "virus_genbank186.fasta";#包括全部参考序列的文件名称（FASTA格式）

our $adapter;         #adapter序列,所有的数据都是同一个adapter
our $adapterLen=11;  #原始adapter序列的前adapter_offset个bp，否则就用adapter全长用来匹配
our $n_Cutoff=1;
our $read_Length=15; 
our $read_PerYield=5e5;#5Mreads*4*100=2G字节
 
our $max_dist = 1;  #bwa允许的最大编辑距离 
our $max_open = 1;  #bwa允许的最大gap数量
our $max_extension = 1; #bwa允许的最大gap长度,-1表示不允许长gap
our $len_seed=15; #bwa中的种子区长度,100bp测序用50，50bp测序用36
our $dist_seed = 1; #bwa种子区允许的最大编辑距离
our $thread_num = 8; #bwa程序调用的线程数量 

################################
##   设置所有目录和文件的路径 ##
################################
our $WORKING_DIR=cwd();#工作目录就是当前目录
our $DATABASE_DIR=$WORKING_DIR."/databases";#所有数据库文件所在的目录
our $BIN_DIR=$WORKING_DIR."/bin";#所有可执行文件所在的目录
our $Tools=$WORKING_DIR."/tools";#所有作为工具的可执行文件所在的目录

##################
## 程序参数处理 ##
##################
&GetOptions( 'file_list=s' => \$file_list,#包括所有待处理的样本文件名称（无后缀）
	'suffix=s' => \$suffix,
	'adapter=s' => \$adapter,
	'adapterLen=i' => \$adapterLen,
	'run_R!' => \$run_R,	
	'n_Cutoff=i' => \$n_Cutoff,
	'read_Length=i' => \$read_Length,
	'read_PerYield=f' => \$read_PerYield,	
	'rRNA_removal!' => \$rRNA_removal,
	'rRNA_reference=s' => \$rRNA_reference,
	'virus_removal!' => \$virus_removal,	
	'virus_reference=s' => \$virus_reference,
	'max_dist=i' => \$max_dist,
	'max_open=i' => \$max_open,
	'max_extension=i' => \$max_extension,
	'len_seed=i' => \$len_seed,
	'dist_seed=i' => \$dist_seed,			 
	'thread_num=i' => \$thread_num
			 );

unless ($file_list) {
	die $usage;
}

#################
##  主程序开始 ##
################
main: {
    #这部分是clipper过程，如果不需要可以注释掉
	system("$Tools/fastq_clipper.pl --filelist $file_list --suffix $suffix --adapter $adapter --adapterLen $adapterLen --distance 1");
	system("$Tools/fastq_clipper.pl --filelist $file_list --suffix unmatched1 --adapter $adapter --adapterLen $adapterLen --distance 2");
	system("$BIN_DIR/files_combine1.pl --filelist $file_list --suffix trimmed --startNumber 1 --endNumber 2");
	system("rm *.trimmed1");
	system("rm *.trimmed2");
	system("rm *.null1");
	system("rm *.null2");
	system("rm *.unmatched1");
	system("rm *.unmatched2");
    #先去除含有1个"N"以上的reads，最后去除短reads（<15bp）
	if($run_R){
		system("Rscript $BIN_DIR/run_sRNA_clean.R filelist=$file_list nCutoff=$n_Cutoff readLength=$read_Length RdPerYield=$read_PerYield");
 	}       
	#如果需要，就去除rRNA污染
	if($rRNA_removal){	
		system("$BIN_DIR/bwa_remove.pl --file_list $file_list --reference $DATABASE_DIR/$rRNA_reference --max_dist $max_dist --max_open $max_open --max_extension $max_extension --len_seed $len_seed --dist_seed $dist_seed --thread_num $thread_num");
		system("$BIN_DIR/files_name_change.pl --file_list $file_list --suffix1 unmapped --suffix2 clean");
	}
	#如果需要，就去除virus污染
	if($virus_removal){	
		system("$BIN_DIR/bwa_remove.pl --file_list $file_list --reference $DATABASE_DIR/$virus_reference --max_dist $max_dist --max_open $max_open --max_extension $max_extension --len_seed $len_seed --dist_seed $dist_seed --thread_num $thread_num");
		system("$BIN_DIR/files_name_change.pl --file_list $file_list --suffix1 unmapped --suffix2 clean");
	}
}