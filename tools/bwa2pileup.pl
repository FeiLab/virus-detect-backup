#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;
use File::Basename;
my $usage = <<_EOUSAGE_;

#########################################################################################
# bwa2pileup.pl --file_list <FILE> --reference <FILE> --index_type <String> --out_type <String> 
#                 --max_dist[INT] --max_open [INT] --max_extension [INT] --len_seed [INT] --dist_seed [INT] --thread_num [INT]
# 
# Required(3):
#  --file_list A txt file containing a list of input file names without any suffix
#  --reference A fasta file containing all the reference sequences
#  --out_type  output mapped, sam, bam, pileup or contig files
# 
# Required(1):
#  --index_type  
#
# BWA-related options(6):
#  --max_dist      Maximum edit distance [1]  
#  --max_open      Maximum number of gap opens [1]  
#  --max_extension Maximum number of gap extensions [1]  
#  --len_seed      Take the first INT subsequence as seed [15] 
#  --dist_seed     Maximum edit distance in the seed [1]  
#  --thread_num    Number of threads (multi-threading mode) [8]
###########################################################################################

_EOUSAGE_

	;
#################
##   全局变量  ##
#################
our $file_list;#包括所有待处理的样本的文本文件名称（无后缀）
our $reference;#包括全部参考序列的文件名称（FASTA格式）
our $index_name;#参考序列的索引名称
our $index_type = "is";  #bwa索引文件类型
our $out_type;  #输出文件类型
our $max_dist = 1; # 
our $max_open = 1; #
our $max_extension = 1; # 
our $len_seed = 15; #
our $dist_seed = 1; # 
our $thread_num = 8; # 

################################
##   设置所有目录和文件的路径 ##
################################
our $WORKING_DIR=cwd();#工作目录就是当前目录
our $BIN_DIR=$WORKING_DIR."/bin";#所有可执行文件所在的目录

##################
## 程序参数处理 ##
##################
&GetOptions( 'file_list=s' => \$file_list,
	'reference=s' => \$reference,
	'out_type=s' => \$out_type,
	'max_dist=i' => \$max_dist,
	'max_open=i' => \$max_open,
	'max_extension=i' => \$max_extension,
	'len_seed=i' => \$len_seed,
	'dist_seed=i' => \$dist_seed,			 
	'thread_num=i' => \$thread_num 
			 );

unless ($file_list&&$reference&&$out_type) {#这3个参数必须通过输入得到
	die $usage;}
$index_name = basename($reference);#去掉目录名称，只保留文件名称
$index_name =~ s/\.\S*$//;#去掉文件后缀名

#################
##  主程序开始 ##
#################
main: {
    #调用bwa为参考序列建立索引,用于下一步的alignment
    &process_cmd("$BIN_DIR/bwa index -p $index_name -a $index_type $reference 2> bwa.log") unless (-e "$index_name.amb");#建立索引文件，aligment用，运行完程序也不删除   
    &process_cmd("$BIN_DIR/samtools faidx $reference") unless (-e "$reference.fai");#建立索引文件，后面bam格式转换用，运行完程序也不删除   
    my $sample;
    my $i=0;
    open(IN, "$file_list");

    while (<IN>) {
		$i=$i+1;
		chomp;
		$sample=$_; #每次循环读入一行，后续代码都是处理该样本文件（名称无后缀）。
		print "#processing sample $i by $0: $sample\n";
		
		#aligment -> sam -> bam -> sorted bam -> pileup
		&process_cmd("$BIN_DIR/bwa aln -n $max_dist -o $max_open -e $max_extension -i 0 -l $len_seed -k $dist_seed -t $thread_num $index_name $sample.clean 1> $sample.sai 2>> bwa.log") unless (-s "$sample.sai");
		&process_cmd("$BIN_DIR/bwa samse -n 10000 -s $index_name $sample.sai $sample.clean 1> $sample.pre.sam 2>> bwa.log") unless (-s "$sample.pre.sam");			
		system("rm $sample.sai");
		&process_cmd("$BIN_DIR/SAM_filter_out_unmapped_reads.pl $sample.pre.sam $sample.unmapped $sample.mapped > $sample.sam1") unless (-s "$sample.sam1");
		if($out_type eq "mapped"){next;}		
		&process_cmd("$BIN_DIR/samFilter.pl $sample.sam1 > $sample.sam") unless (-s "$sample.sam");#只保留最好一个级别的的hits
		if($out_type eq "sam"){next;}
		&process_cmd("$BIN_DIR/samtools view -bt $reference.fai $sample.sam > $sample.bam") unless (-s "$sample.bam");
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
	print "All the input files have been processed by $0\n";
}
system("rm *.sam1");	
#################
##    子程序   ##
#################
sub process_cmd {
	my ($cmd) = @_;	
	print "CMD: $cmd\n";
	my $ret = system($cmd);	#成功就返回0，否则就失败
	if ($ret) {
		print "Error, cmd: $cmd died with ret $ret";
	}
	return($ret);
}