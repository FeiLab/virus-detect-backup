#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;
use File::Basename;
my $usage = <<_EOUSAGE_;

#########################################################################################
# paired-bowtie2pileup.pl --file_list <FILE> --reference [FILE]  --out_type <String> --mismatchs [INT]
#                  --len_seed [INT] --dist_seed [INT] --quality_sum [INT] --thread_num [INT]                
# Required:
#  --file_list The name of a txt file containing a list of input file names without any suffix
#  --reference The name of a fasta file containing the reference sequences
#  --out_type  output mapped, sam, bam or pileup files
#
# Bowtie-related options(6): 
#  --mismatchs     alignments may have no more than v mismatches(v<=3), ignoring qualities
#  --len_seed      Take the first INT subsequence as seed [36] 
#  --dist_seed     Maximum edit distance in the seed [1] 
#  --quality_sum   Max sum of mismatch quals across alignment in the seed region [70]  
#  --thread_num    Number of threads (multi-threading mode) [8]  
##########################################################################################

_EOUSAGE_
	;
###############################
##   全局变量（含默认设置）  ##
###############################
our $file_list;#包括所有待处理的输入文件（数据）的列表文件（无后缀）
our $reference;#包括全部参考序列的文件名称（FASTA格式）
our $index_name;#参考序列的索引名称
our $out_type;  #输出文件类型
our $mismatchs; #全局错配
our $len_seed = 36; #种子区长度
our $dist_seed = 1; #种子区允许的最大编辑距离
our $quality_sum = 70; #种子区所有错配位点的质量Phred quality值之和上限
our $thread_num = 8; #程序调用的线程数量 

################################
##   设置所有目录和文件的路径 ##
################################
our $WORKING_DIR=cwd();#工作目录就是当前目录
our $BIN_DIR=$WORKING_DIR."/bin";#所有可执行文件所在的目录

#################
## 输入参数处理##
#################
&GetOptions( 'file_list=s' => \$file_list,#包括所有待处理的样本的文本文件名称（无后缀）
		'reference=s' => \$reference,#参考转录组的文件名称（FASTA格式）
		'out_type=s' => \$out_type,
		'mismatchs=i' => \$mismatchs,
		'len_seed=i' => \$len_seed,
		'dist_seed=i' => \$dist_seed,	
		'quality_sum=i' => \$quality_sum,		
		'thread_num=i' => \$thread_num
			 );

unless ($file_list&&$reference&&$out_type) {#至少需要2个参数
	die $usage;
}
$index_name = basename($reference);#去掉目录名称，只保留文件名称
$index_name =~ s/\.\S*$//;#去掉文件后缀名

#################
##  主程序开始 ##
#################
main: {
    #调用bowtie为参考序列建立索引,用于下一步的alignment
    &process_cmd("$BIN_DIR/bowtie-build -f $reference $index_name") unless (-e "$index_name.ebwt");#建立索引文件，aligment用，运行完程序也不删除   
    &process_cmd("$BIN_DIR/samtools faidx $reference") unless (-e "$reference.fai");#建立索引文件，后面bam格式转换用，运行完程序也不删除      
    my $sample;
    my $i=0;
    open(IN, "$file_list");

    while (<IN>) {
		$i=$i+1;
		chomp;
		my @paired_seq = split(/\t/, $_);
		$sample="a"; #每次循环读入一行，后续代码都是处理该样本文件（名称无后缀）。
		print "#processing sample $i by $0: $sample\n";
		
		#aligment -> sam -> bam -> sorted bam -> pileup
		if($mismatchs){
			&process_cmd("$BIN_DIR/bowtie $index_name -v $mismatchs -p $thread_num -I 50 -X 600 -1 paired_seq[1].clean -2 paired_seq[2].clean -S -a $sample.pre.sam") unless (-s "$sample.pre.sam");		

		}
		else{
			&process_cmd("$BIN_DIR/bowtie $index_name -n $dist_seed -l $len_seed -e $quality_sum -p $thread_num -I 50 -X 600 -1 paired_seq[1].clean -2 paired_seq[2].clean -S -a $sample.pre.sam") unless (-s "$sample.pre.sam");	
		}
		&process_cmd("$BIN_DIR/SAM_filter_out_unmapped_reads.pl $sample.pre.sam $sample.unmapped $sample.mapped > $sample.sam") unless (-s "$sample.sam");	
		if($out_type eq "mapped"){next;}
		&process_cmd("samtools view -bt $reference.fai $sample.sam > $sample.bam") unless (-s "$sample.bam");
		if($out_type eq "sam"){next;}
		&process_cmd("samtools sort $sample.bam $sample.sorted") unless (-s "$sample.sorted.bam");
		if($out_type eq "bam"){next;}
		&process_cmd("$BIN_DIR/samtools sort $sample.bam $sample.sorted") unless (-s "$sample.sorted.bam");
		&process_cmd("$BIN_DIR/samtools mpileup -f $reference $sample.sorted.bam > $sample.pileup") unless (-s "$sample.pileup");
		if($out_type eq "pileup"){next;}
		my $file_size= -s "$sample.pileup";#根据pileup文件大小是不是0，进入下面处理流程
		if($file_size!=0){#如果文件大小不是0，
			&process_cmd("java -cp $BIN_DIR extractConsensus $sample 0 40 1");#提取连续片段（depth>=1,length>=40），文件名含有1			
		}
	}
	close(IN);
	print "###############################\n";
	print "All the input files have been processed by the $0\n";
}
system("rm *.ebwt");
####
sub process_cmd {
	my ($cmd) = @_;	
	print "CMD: $cmd\n";
	my $ret = system($cmd);	
	if ($ret) {
		die "Error, cmd: $cmd died with ret $ret";
	}
	return($ret);
}
