#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;
my $usage = <<_EOUSAGE_;

#########################################################################################
# $0 --file_list <FILE> --file_type [String] --reference [FILE] --coverage [Float]
#                 --host_removal --host_reference [FILE] --objective_type [String]
#                 --max_dist[INT] --max_open [INT] --max_extension [INT] --len_seed [INT] --dist_seed [INT] --thread_num [INT]
#                 --strand_specific --min_overlap [INT] --max_end_clip [INT] --cpu_num [INT] --mis_penalty [INT] --gap_cost [INT] --gap_extension [INT]
#                                  
# Required(1):
#  --file_list The name of a txt file containing a list of input file names without any suffix
#  
# Options(6):
#  --file_type The format of input files(fastq or fasta)  [fastq]
#  --reference The name of a fasta file containing all of the virus reference sequences  [vrl_genbank.fasta]
#  --coverage  The reference sequence with at least a pencetage of its length covered by reads is detected as a aligned virus contig  [0.3]
#  --host_removal         Remove host sequences before assembly or not  [Not selected]
#  --host_reference       If using host_removal, the name of a host reference file should be provided  [Null]
#  --objective_type       The objective type of optimization(maxLen, n50 or avgLen)  [maxLen]
# 
# BWA-related options(6):
#  --max_dist      Maximum edit distance [1]  
#  --max_open      Maximum number of gap opens [1]  
#  --max_extension Maximum number of gap extensions [1]  
#  --len_seed      Take the first INT subsequence as seed [15] 
#  --dist_seed     Maximum edit distance in the seed [1]  
#  --thread_num    Number of threads (multi-threading mode) [8] 
# 
# Megablast-related options(7):
#  --strand_specific Only for sequences assembled from strand specific RNA-seq  [Not selected]
#  --min_overlap     The minimum overlap length between two contigs to be combined [30]
#  --max_end_clip    The maximum length of end clips [6]
#  --cpu_num         Number of processors to use [8] 
#  --mis_penalty     Penalty for a nucleotide mismatch [-1]
#  --gap_cost        Cost to open a gap [2] 
#  --gap_extension   Cost to extend a gap [1] 
###########################################################################################

_EOUSAGE_

	;

################################
##   设置所有目录和文件的路径 ##
################################
our $WORKING_DIR=cwd();#工作目录就是当前目录
our $DATABASE_DIR=$WORKING_DIR."/databases";#所有数据库文件所在的目录
our $BIN_DIR=$WORKING_DIR."/bin";#所有可执行文件所在的目录

###############################
##   全局变量（含默认设置）  ##
###############################
our $file_list;#包括所有待处理的输入文件（数据）的列表文件（无后缀）
our $file_type= "fastq";#输入文件的类型，当前只支持fastq和fasta格式
our $reference= "vrl_genbank.fasta";#包括全部参考序列的文件名称（FASTA格式）
our $coverage=0.3;  #每条参考序列如果被reads覆盖的部分占全长比例的阈值
our $host_removal;         #是否需要remove host
our $host_reference;       #如果要remove host，必须提供一个文件包括host的全部参考序列
our $objective_type='maxLen';#优化目标值的类型，只有n50、maxLen和avgLen三种

our $max_dist = 1;  #bwa允许的最大编辑距离 
our $max_open = 1;  #bwa允许的最大gap数量
our $max_extension = 1; #bwa允许的最大gap长度,-1表示不允许长gap
our $len_seed = 15; #bwa中的种子区长度
our $dist_seed = 1; #bwa种子区允许的最大编辑距离
our $thread_num = 8; #bwa程序调用的线程数量 
	 
our $strand_specific;  #专门用于strand specific转录组
our $min_overlap = 30; #hsp合并时，最短的overlap
our $max_end_clip = 6; #hsp合并时，两端允许的最小clip
our $cpu_num = 8;          #megablast使用的cpu数目
our $mis_penalty = -1;     #megablast中，对错配的罚分，必须是负整数
our $gap_cost = 2;         #megablast中，对gap open的罚分，必须是正整数
our $gap_extension = 1;    #megablast中，对gap open的罚分，必须是正整数

our $input_suffix='clean'; #数据文件后缀（clean或unmapped）
#################
## 输入参数处理##
#################
&GetOptions( 'file_list=s' => \$file_list,
		'file_type=s' => \$file_type,
		'reference=s' => \$reference,
		'coverage=f' => \$coverage,
		'host_removal!' => \$host_removal,
		'host_reference=s' => \$host_reference,
		'objective_type=s' => \$objective_type,
		
		'max_dist=i' => \$max_dist,
		'max_open=i' => \$max_open,
		'max_extension=i' => \$max_extension,
		'len_seed=i' => \$len_seed,
		'dist_seed=i' => \$dist_seed,			 
		'thread_num=i' => \$thread_num,
		
		'strand_specific!' => \$strand_specific,
		'min_overlap=i' => \$min_overlap,
		'max_end_clip=i' => \$max_end_clip,
		'cpu_num=i' => \$cpu_num ,
		'mis_penalty=i' => \$mis_penalty,
		'gap_cost=i' => \$gap_cost,
		'gap_extension=i' => \$gap_extension
			 );

unless ($file_list) {#至少需要1个参数
	die $usage;
}

#################
##  主程序开始 ##
#################
main: {
	#删除前一次的结果文件

	opendir(DIR, ".") || die "Can't open directory $WORKING_DIR";
	my @files = readdir(DIR);
	foreach (@files){#所有的文件或者子目录
		if (/^aligned$/) {system("rm -r aligned");}
		if (/^assembled$/) {system("rm -r assembled");}
		if (/^combined$/) {system("rm -r combined");}
		if (/\.finished$/) {system("rm $_");}
		if (/\.result$/) {system("rm $_");}
		if (/\.output$/) {system("rm $_");}
		if (/\.out$/) {system("rm $_");}
	}
	closedir DIR; 
	
	#detect 已知病毒
	system("mkdir aligned");
	system("$BIN_DIR/bwa-alignAndCorrect.pl --file_list $file_list --reference $DATABASE_DIR/$reference --coverage $coverage --max_dist $max_dist --max_open $max_open --max_extension $max_extension --len_seed $len_seed --dist_seed $dist_seed --thread_num $thread_num");  

	if($strand_specific){
	
		system("$BIN_DIR/removeRedundancy_batch.pl --file_list $file_list --file_type $file_type --input_suffix contigs1.fa --contig_type aligned --contig_prefix KNOWN --min_overlap $min_overlap --max_end_clip $max_end_clip --cpu_num $cpu_num --mis_penalty $mis_penalty --gap_cost $gap_cost --gap_extension $gap_extension");
	}
	else{
		system("$BIN_DIR/removeRedundancy_batch.pl --file_list $file_list --file_type $file_type --input_suffix contigs1.fa --contig_type aligned --contig_prefix KNOWN --strand_specific --min_overlap $min_overlap --max_end_clip $max_end_clip --cpu_num $cpu_num --mis_penalty $mis_penalty --gap_cost $gap_cost --gap_extension $gap_extension");
	}
	
	#detect 未知病毒
	system("mkdir assembled");
	if($host_removal){	
		system("$BIN_DIR/bwa_remove.pl --file_list $file_list --reference $DATABASE_DIR/$host_reference --max_dist $max_dist --max_open $max_open --max_extension $max_extension --len_seed $len_seed --dist_seed $dist_seed --thread_num $thread_num");
		$input_suffix='unmapped';#以后处理的文件后缀名称由clean变为unmapped
		#这里要注意，从sam文件中提取的数据都是fastq文件，无论原始文件是fastq还是fasta
		system("$BIN_DIR/Velvet_Optimiser.pl --file_list $file_list --input_suffix $input_suffix  --file_type fastq --objective_type $objective_type --hash_end 19 --coverage_end 25");
		system("$BIN_DIR/runVelvet.pl --parameters optimization.result --input_suffix $input_suffix --file_type fastq --output_suffix contigs1.fa");
		system("rm *.unmapped");#得到contig1.fa，就不需要unmapped文件了
	}	
	else{
		system("$BIN_DIR/Velvet_Optimiser.pl --file_list $file_list --input_suffix $input_suffix  --file_type $file_type --objective_type $objective_type --hash_end 19 --coverage_end 25");
		system("$BIN_DIR/runVelvet.pl --parameters optimization.result --input_suffix $input_suffix --file_type $file_type --output_suffix contigs1.fa");
	}
	if($strand_specific){
	
		system("$BIN_DIR/removeRedundancy_batch.pl --file_list $file_list --file_type $file_type --input_suffix contigs1.fa --contig_type assembled --contig_prefix NOVEL --min_overlap $min_overlap --max_end_clip $max_end_clip --cpu_num $cpu_num --mis_penalty $mis_penalty --gap_cost $gap_cost --gap_extension $gap_extension");
	}
	else{
		system("$BIN_DIR/removeRedundancy_batch.pl --file_list $file_list --file_type $file_type --input_suffix contigs1.fa --contig_type assembled --contig_prefix NOVEL --strand_specific --min_overlap $min_overlap --max_end_clip $max_end_clip --cpu_num $cpu_num --mis_penalty $mis_penalty --gap_cost $gap_cost --gap_extension $gap_extension");
	}

	#病毒合并
	system("mkdir combined");
	system("$BIN_DIR/files_combine2.pl --filelist $file_list --folder1 aligned --folder2 assembled");
	system("$BIN_DIR/removeRedundancy_batch.pl --file_list $file_list --file_type $file_type --input_suffix contigs1.fa --contig_type combined --contig_prefix CONTIG --strand_specific --min_overlap $min_overlap --max_end_clip $max_end_clip --cpu_num $cpu_num --mis_penalty $mis_penalty --gap_cost $gap_cost --gap_extension $gap_extension");
	print "All the samples have been processed by $0\n";
	#这个输出到log，表示运行结束
	print "###############################\n";
	print "###############################\n\n";
	#建立这个文件，也表示运行结束
	system("touch virus_detect.run.finished");
}