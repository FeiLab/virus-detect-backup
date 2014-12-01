#!/usr/bin/perl -w 
use strict; 
use IO::File; 
use Getopt::Long; #这个模块可以接受完整参数
use Cwd;
my $usage = <<_EOUSAGE_;

#########################################################################################
# hit_filter2.pl --input <FILE> --diff_ratio <FLOAT> --output <FILE>
#
# Required[1]:
#  --input              A name of an input file containing sequences in fasta format
#
# options[1]:
#   
##########################################################################################

_EOUSAGE_
	;
#################
##   全局变量  ##
#################	
our $input;            #需要处理的文件
our $diff_ratio;     
our $output;            

################################
##   设置所有目录和文件的路径 ##
################################
our $WORKING_DIR=cwd();#工作目录就是当前目录
our $BIN_DIR=$WORKING_DIR."/bin";#所有可执行文件所在的目录

##################
## 程序参数处理 ##
##################
&GetOptions( 'input=s' => \$input, 
	'diff_ratio=f' => \$diff_ratio,	
	'output=s' => \$output 	
			 );
unless ($input) {#这个参数通过输入得到
die $usage;}
			 
#################
##  主程序开始 ##
#################
			 
my (@all_data, $name, $sequence); 
$name = ''; $sequence = ''; 

open(IN, "$input");
while (<IN>) {
	chomp; 
	my @each_line = split(/\t/, $_);
    my @contigs= split(/,/, $each_line[4]);#第5列是所有contig的名称
    my $covered_bp=$each_line[7];          #第8列是read covered bp总数，不是contig覆盖的bp总数
    if ($each_line[7] eq ""){$covered_bp=0;}	
	push(@all_data, [scalar(@contigs), $covered_bp, $_]);#contig数量，covered bp总数和整行数据
}
close(IN);

my @inset;  #用于存储非冗余序列，这样会再次存储@all_data中的序列，应该改进为只保存index
my @restset = ''; #用于存储冗余序列，这些数据的保存仅为了校对用
@all_data = sort { -1*($a->[0] <=> $b->[0]) || -1*($a->[1] <=> $b->[1])} @all_data; #@all_data中数据按照contig数量，然后covered bp总数降序排序
=head;#调试时，输出看一下
for my $tr (@all_data) {
 print $tr->[2]."\n";
}
=cut;
my $contig_count=1; 
for my $tr (@all_data) {
	if (scalar(@inset) == 0) {#第一次
		push(@inset, $tr->[2]); #把contig数量最多的一条hit先放进来
	}else{
		my @aa = split(/\t/, $tr->[2]);
		#print $tr->[2]."\n";
		my $return_string = &ifRedundant(\@inset, \$aa[4]);
		if ($return_string eq "n") {#如果query是非冗余，包括进来
			push(@inset, $tr->[2]); 
		}else{#如果query是冗余，扔掉
			push(@restset, $tr->[2]); 
		}
	}
}

open(OUT1, ">$output");#将所有非冗余序列输出
for my $tr (@inset) {
	print OUT1 $tr."\n";
}
close(OUT1);

=head;#没用，就不输出了
open(OUT2, ">restset");#将所有冗余序列输出
for my $tr (@restset) {
	print OUT2 $tr."\n";
}
close(OUT2);
=cut;
#######################
##     子程序开始    ##
#######################
sub ifRedundant {
	my ($inset, $query) = @_; 
	my @query_contigs= split(/,/, $$query);
	my $ratio;
	for my $tr (@$inset) {
		my @aa = split(/\t/, $tr);
		my @contigs= split(/,/, $aa[4]);
		my %contigs;
		for my $each_contig (@contigs) {#首先把本次所有的
			$contigs{$each_contig}=1;
		}
		my $total_contigs=0;#注意每次必须初始化
		my $diff_contigs=0; #注意每次必须初始化
		for my $each_contig (@query_contigs) {
			$total_contigs++;
			if(not defined $contigs{$each_contig}){
				$diff_contigs++;			
			}
		}
		$ratio=$diff_contigs*1.0/$total_contigs;
		if($ratio<=$diff_ratio){return "r";} #一旦发现有相似contig，就认为是冗余，立刻返回
	}
	return "n";#如果inset中的所有数据比较一遍后，没有发现相似的，就认为不是冗余的记录
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