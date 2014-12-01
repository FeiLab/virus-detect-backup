#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;
#先通过多次循环去冗余，然后再做base correction
my $usage = <<_EOUSAGE_;

#########################################################################################################
# removeRedundancy_batch.pl --file_list <FILE> --file_type <String> --input_suffix <String> --contig_type <String> --contig_prefix <String> 
#                 --strand_specific --min_overlap [INT] --max_end_clip [INT] --cpu_num [INT] --mis_penalty [INT] --gap_cost [INT] --gap_extension [INT]
#
# Required(4):
#  --file_list       A txt file containing a list of input file names without any suffix
#  --file_type       The format of input files(fastq or fasta)
#  --input_suffix    A file name suffix of input data 
#  --contig_type     The contig type for redundancy removal (aligned or assembled)
#  --contig_prefix   A prefix of each contig name in the data (fasta format)
#
# Megablast-related options(7):
#  --strand_specific Only for sequences assembled from strand specific RNA-seq [Not selected]
#  --min_overlap     The minimum overlap length between two contigs to be combined [30]
#  --max_end_clip    The maximum length of end clips [4]
#  --cpu_num         Number of processors to use [8] 
#  --mis_penalty     Penalty for a nucleotide mismatch [-1]
#  --gap_cost        Cost to open a gap [2] 
#  --gap_extension   Cost to extend a gap [1]  
#########################################################################################################

_EOUSAGE_
	;
	
#################
##   全局变量  ##
#################
our $file_list; #包括所有待处理的样本的文本文件名称（无后缀）
our $file_type;#输入文件的类型，当前只支持fastq和fasta格式
our $input_suffix;#输入文件的后缀名称
our $contig_type; #contig类型，决定了输出文件夹是aligned还是assembled
our $contig_prefix;#fasta文件改名时，每条记录用到的前缀

our $strand_specific;  #专门用于strand specific转录组
our $min_overlap = 30; #hsp合并时，最短的overlap
our $max_end_clip = 6; #hsp合并时，两端允许的最小clip
our $cpu_num = 8;          #megablast使用的cpu数目
our $mis_penalty = -1;     #megablast中，对错配的罚分，必须是负整数
our $gap_cost = 2;         #megablast中，对gap open的罚分，必须是正整数
our $gap_extension = 1;    #megablast中，对gap open的罚分，必须是正整数
################################
##   设置所有目录和文件的路径 ##
################################
our $WORKING_DIR=cwd();#工作目录就是当前目录
our $DATABASE_DIR=$WORKING_DIR."/databases";#所有数据库文件所在的目录
our $BIN_DIR=$WORKING_DIR."/bin";#所有可执行文件所在的目录
our $CONTIG_DIR;#指定的病毒contig所在的目录

#################
## 程序参数处理##
#################
&GetOptions( 'file_list=s' => \$file_list,
	'file_type=s' => \$file_type,
	'input_suffix=s' => \$input_suffix,
	'contig_type=s' => \$contig_type,
	'contig_prefix=s' => \$contig_prefix,
	'strand_specific!' => \$strand_specific,
	'min_overlap=i' => \$min_overlap,
	'max_end_clip=i' => \$max_end_clip,
	'cpu_num=i' => \$cpu_num,
	'mis_penalty=i' => \$mis_penalty,
	'gap_cost=i' => \$gap_cost,
	'gap_extension=i' => \$gap_extension			 			 
			 );
			 
unless ($file_list&&$input_suffix&&$contig_type&&$contig_prefix) {#这4个参数通过输入得到
die $usage;}

$CONTIG_DIR=$WORKING_DIR."/$contig_type";
#################
##  主程序开始 ##
#################
open(IN1,$file_list) || die "Can't open the file $file_list\n";
my $sample;
my $j=0;
while(<IN1>){
	$j=$j+1;
	chomp;
	$sample=$_; #每次循环读入一行，后续代码都是处理该样本文件（名称无后缀）。
	print "#processing sample $j by $0: $sample\n";
		
	my $file_size= -s "$CONTIG_DIR/$sample.$input_suffix";#根据该样本文件是不是0，进入下面处理流程
	if($file_size==0){#如果文件大小是0，跳过下面的指令直接进入下一次循环
		system ("mv $CONTIG_DIR/$sample.$input_suffix $CONTIG_DIR/$sample.$contig_type.fa");
		next;
	}
	#否则，移动到当前目录来去冗余
	system ("$BIN_DIR/dust $CONTIG_DIR/$sample.$input_suffix 1> $sample.masked 2> dust.log");#首先mask简单重复序列
	#去掉简单重复序列，可以大大节省序列合并的时间
	system ("$BIN_DIR/trim_XNseq1.pl $sample.masked $CONTIG_DIR/$sample.$input_suffix 0.8 40 > $sample.$input_suffix");
	system ("rm $sample.masked");
	system ("rm $CONTIG_DIR/$sample.$input_suffix");
	my $contig_num1 =  `grep \'>\' $sample.$input_suffix | wc -l `;#首先得到去冗余前的序列总数
	chomp($contig_num1);
	my $contig_num2 = 0;#此变量保存每轮去除冗余后的序列总数	
	my @a= split(/\./,$input_suffix);
	my ($i)= $a[0] =~ /(\d+)/;#提取contigs开始的编号,默认是1
	
	
	while($contig_num2!=$contig_num1)#如果新contigs数量不等于旧数量，继续计算，直到不变为止
	{
		process_cmd ("$BIN_DIR/removeRedundancy.pl --input $sample.contigs$i.fa --min_overlap $min_overlap --max_end_clip $max_end_clip --cpu_num $cpu_num --mis_penalty $mis_penalty --gap_cost $gap_cost --gap_extension $gap_extension");		
		process_cmd("rm $sample.contigs$i.fa");#删除旧文件
		$i++;
		$contig_num1=$contig_num2;#更新contig_num1
		#更新contig_num2
		$contig_num2 =  `grep \'>\' inset | wc -l `;#校正后生成的contigs
		chomp($contig_num2);	
		system ("mv inset $sample.contigs$i.fa");	
	}
	#去冗余完毕，开始进行base correction
	my $sample_reference= $sample.".contigs$i.fa";
	
	#这里需要检测结果文件是否为空，处理无法组装出新序列的情况
	system "wc -l $sample.contigs$i.fa";
	my $content_sig =  `grep \'>\' $sample.contigs$i.fa | wc -l `;#校正后生成的contigs		
	chomp($content_sig);
	# 如果为空，则跳过该样本
	if (!$content_sig){
		print "Skip sample:$sample\n";
		# 要确保一下两个文件必须存在
		system "touch ./aligned/$sample.aligned.fa" unless (-e "./aligned/$sample.aligned.fa");
		system "touch ./assembled/$sample.assembled.fa" unless (-e "./assembled/$sample.assembled.fa");
		system("rm $sample.contigs$i.fa");
		next;
	}
		
	#去除冗余的结果文件，产生新contigs，将作为参考，进行校正
	my $sample_reads= $sample.".clean";#read文件，需要align到新contigs上	
	#aligment -> sam -> bam -> sorted bam -> pileup
	my $format="-q";
	if($file_type eq "fasta"){$format="-f"};
	&process_cmd("$BIN_DIR/bowtie-build --quiet -f $sample_reference $sample") unless (-e "$sample.1.amb");
	&process_cmd("$BIN_DIR/samtools faidx $sample_reference") unless (-e "$sample_reference.fai");
	&process_cmd("$BIN_DIR/bowtie --quiet $sample -v 1 -p $cpu_num $format $sample.clean -S -a --best $sample.sam") unless (-s "$sample.sam");
	&process_cmd("$BIN_DIR/samtools view -bt $sample_reference.fai $sample.sam > $sample.bam") unless (-s "$sample.bam");
	&process_cmd("$BIN_DIR/samtools sort $sample.bam $sample.sorted") unless (-s "$sample.sorted.bam");
	&process_cmd("$BIN_DIR/samtools mpileup -f $sample_reference $sample.sorted.bam > $sample.pileup") unless (-s "$sample.pileup");	

	$file_size= -s "$sample.pileup";#根据pileup文件大小是不是0，进入下面处理流程
	if($file_size==0){#如果文件大小是0，需要建立一个空文件，然后跳出循环
		system ("touch $CONTIG_DIR/$sample.$contig_type.fa");
		next;
	}
	#否则继续下面处理
	$i++;
	&process_cmd("java -cp $BIN_DIR extractConsensus $sample 1 40 $i");
	system("$BIN_DIR/renameFasta.pl --inputfile $sample.contigs$i.fa --outputfile contigs --prefix $contig_prefix");#每条序列的名字要统一命名以去除重复
	system ("mv contigs $CONTIG_DIR/$sample.$contig_type.fa");#去冗余得到的最终文件转回到$CONTIG_DIR目录
	
	#删除本次循环中产生的中间文件
	system("rm $sample.sam");
	system("rm $sample.bam");
	system("rm $sample.sorted.bam");
	system("rm $sample.pileup");#每次必须删除，否则下次不能继续
	system("rm $sample_reference");
	system("rm $sample_reference.fai");
	system("rm *.ebwt");
	system("rm $sample.contigs$i.fa");	
}
close(IN1);
print "All the samples have been processed by $0\n";
print "###############################\n";
system("touch removeRedundancy_batch.$contig_type.finished");#建立这个文件，表示结束标志
system("rm tem.*");
system("rm tem");
system("rm *.log");
system("rm restset");
system("rm query");


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