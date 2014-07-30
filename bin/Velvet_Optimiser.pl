#!/usr/bin/env perl
use strict;
use Getopt::Long qw(:config no_ignore_case bundling);
use Bio::SeqIO;
use Cwd;
my $usage = <<_EOUSAGE_;

#########################################################################################
# Velvet_Optimiser.pl --file_list <File> --input_suffix <String> --file_type <String> --velvet_dir [String] 
                      --objective_type [String] --hash_end <INT> --coverage_end <INT> 
# Required[4]:
#  --file_list             The name of a txt file containing a list of input file names without any suffix
#  --input_suffix          The suffix name of all the input files
#  --hash_end              The max k-mer length for the velvet to try to optimize the objective_type
#  --coverage_end          The max coverage cutoff for the velvet to try to optimize the objective_type
#  
# Options[2]:
#  --file_type             The format of input files(fastq or fasta)  [fastq]
#  --velvet_dir            The directory to install the program velvet [workingdir/bin]
#  --objective_type        The objective type of optimization(maxLen, n50 or avgLen)  [maxLen]
###########################################################################################

_EOUSAGE_
	;
	
###############################
##   全局变量（含默认设置）  ##
###############################	
our $file_list;              #包括所有待处理的样本的文本文件名称（无后缀）
our $input_suffix;           #输入数据文件的后缀名（clean或unmapped），注意没有"."
our $file_type="fastq";#输入文件的类型，当前只支持fastq和fasta格式
our $velvet_dir;             #velvet的安装目录
our $objective_type='maxLen';#优化目标值的类型，只有n50、maxLen和avgLen三种
our $hash_end;               #优化所尝试的k-mer长度范围的最大值
our $coverage_end;           #优化所尝试的coverage cutoff范围的最大值

our $hash_start=9;
our $coverage_start=5;

################################
##   设置所有目录和文件的路径 ##
################################
our $WORKING_DIR=cwd();#工作目录就是当前目录
$velvet_dir=$WORKING_DIR."/bin";#设置默认velvet路径; 

#################
## 输入参数处理##
#################
&GetOptions( 'file_list=s' => \$file_list,#包括所有待处理的样本的文本文件名称（无后缀）			 
			'input_suffix=s' => \$input_suffix,
			'file_type=s' => \$file_type,
			'velvet_dir=s' => \$velvet_dir,
			'objective_type=s' => \$objective_type,			 			 
			'hash_end=i' => \$hash_end,
			'coverage_end=i' => \$coverage_end);
			 
unless ($file_list&&$hash_end&&$coverage_end) {
	die $usage;
}

#################
##  主程序开始 ##
#################
main: {
    my $sample;
    my $sampleNum=0;
    my $current_folder;
    my $statfile;
    my $objective;                #每次计算得到的目标值
    my $max_objective;            #目标值的最大值
    my $opt_hash_length=$hash_start; #目标值取得最大值时（最优的）的hash_length
    my $opt_coverage=$coverage_start;#目标值取得最大值时（最优的）的coverage
    my $opt_avgLen=0;#目标值取得最大值时的平均长度
    open(IN, "$file_list");
    open(OUT1, ">optimization.log") or die "fhgfhg\n";#保存计算的中间结果
    open(OUT2, ">optimization.result");#保存最终结果
    while (<IN>) {
		$sampleNum=$sampleNum+1;
		chomp;
		$sample=$_; #每次循环读入一行，后续代码都是处理该样本文件（名称无后缀）。		
		print "#processing sample $sampleNum: $sample\n";		
		
		$max_objective=0;#每次优化开始必须置0
		$opt_hash_length=$hash_start;#每次优化开始必须重置
		$opt_coverage=$coverage_start;#每次优化开始必须重置
		
		#开始优化k-mer length
		for(my $i=$hash_start; $i<= $hash_end; $i=$i+2) {
		    runVelvet($sample,$i,$coverage_start);
			$current_folder =$sample."_".$i."_".$coverage_start;
			$statfile=$current_folder."/contigs.fa";
			my $aa=contigStats($statfile);        #注意返回值是hash表的引用，而不是hash表
			$objective=$aa->{$objective_type};
			print OUT1 $i."\t".$coverage_start."\t".$objective."\t".$aa->{avgLen}."\t".$aa->{numSeqs}."\n";
			#输出每轮的hash长度，coverage(固定的),优化目标值，平均长度，contigs的数目
			if($objective>$max_objective){
			    print OUT1 "yes"."\n";#如果上一列的优化目标值提高了，就写"yes"
				$max_objective=$objective;
				$opt_hash_length=$i;
				$opt_avgLen=$aa->{avgLen};
			}
			&process_cmd("rm $current_folder -r");		
		} 
		#开始优化k-mer length
		for(my $j=$coverage_start+2; $j<=$coverage_end; $j=$j+1) {#注意这里从7开始，因为5上面都算过了
			runVelvet($sample,$opt_hash_length,$j);
			$current_folder=$sample."_".$opt_hash_length."_".$j;
			$statfile=$current_folder. "/contigs.fa";
			my $aa=contigStats($statfile);
			$objective=$aa->{$objective_type};
			print OUT1 $opt_hash_length."\t".$j."\t".$objective."\t".$aa->{avgLen}."\t".$aa->{numSeqs}."\n";
            #输出最优的hash长度（从上面得到），每轮的coverage ,优化目标值，平均长度，contigs的数目
			if($objective>$max_objective){
				print OUT1 "yes"."\n";#如果上一列的优化目标值提高了，就写"yes"
				$max_objective=$objective;
				$opt_coverage=$j;
				$opt_avgLen=$aa->{avgLen};
			}
			&process_cmd("rm $current_folder -r");
		}       
		print OUT2 $sample."\t".$opt_hash_length."\t".$opt_coverage."\t".$max_objective."\t".$opt_avgLen."\n";
	}
	close(IN);
	close(OUT1);
	close(OUT2);
	system("rm velvet.log");
	print "###############################\n";
        print "All the samples have been processed by $0\n";
	system("touch Velvet_Optimiser.run.finished");#建立这个文件，表示结束标志
}

#################
##    子程序   ##
#################
sub runVelvet {
	my $sample1=shift;
	my $hash_length=shift;
	my $cov_cutoff=shift;
	my $outputDir=$sample1."_".$hash_length."_".$cov_cutoff;
	#下面执行command lines
	&process_cmd($velvet_dir."/velveth $outputDir $hash_length -$file_type $sample1.$input_suffix >> velvet.log");
	&process_cmd($velvet_dir."/velvetg $outputDir -cov_cutoff $cov_cutoff -min_contig_lgth 30 >> velvet.log");	
}
sub process_cmd {
	my ($cmd) = @_;	
	print "CMD: $cmd\n";
	my $ret = system($cmd);	
	if ($ret) {
		print "Error, cmd: $cmd died with ret $ret";
	}
	return($ret);
}
sub contigStats {
	
	my $file = shift;
	#my $minsize = shift;
	
	#print "In contigStats with $file, $minsize\n" if $interested;
	
	my $numseq=0;
	my $avglen=0;
	my $minlen=1E9;
	my $maxlen=0;
	my @len;
	my $toosmall=0;
	my $nn=0;
	
	my $in = Bio::SeqIO->new(-file => $file, -format => 'Fasta');
	while(my $seq = $in->next_seq()){
		my $L = $seq->length;
		#if($L < $minsize){
			#$toosmall ++;
			#next;
		#}
		#count Ns
		my $s = $seq->seq;
		my $n = $s =~ s/N/N/gi;
		$n ||= 0;
		$nn += $n;
		#count seqs and other stats
		$numseq ++;
		$avglen += $L;
		$maxlen = $L if $L > $maxlen;
		$minlen = $L if $L < $minlen;
		push @len, $L;
	}
	@len = sort { $a <=> $b } @len;
	my $cum = 0;
	my $n50 = 0;
	for my $i (0 .. $#len){
		$cum += $len[$i];
		if($cum >= $avglen/2) {
			$n50 = $len[$i];
			last;
		}
	}
	
	my %out;
	if($numseq > 0){
		$out{numSeqs} = $numseq;
		$out{numBases} = $avglen;
		$out{numOK} = ($avglen - $nn);
		$out{numNs} = $nn;
		$out{minLen} = $minlen;
		$out{avgLen} = $avglen/$numseq;
		$out{maxLen} = $maxlen;
		$out{n50} = $n50;
		#$out{minsize} = $minsize;
		$out{numTooSmall} = $toosmall;
	}
	else {
		$out{$numseq} = 0;
	}
	
	#print "Leaving contigstats!\n" if $interested;
	return (\%out);#返回hash表的引用
}