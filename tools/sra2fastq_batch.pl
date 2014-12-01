#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;
#先通过多次循环去冗余，然后再做base correction
my $usage = <<_EOUSAGE_;

#########################################################################################################
# sra2fastq_batch.pl --file_list <FILE>
#
# Required(1):
#  --file_list       A txt file containing a list of input file names without any suffix
#########################################################################################################

_EOUSAGE_
	;
	
#################
##   全局变量  ##
#################
our $file_list; #包括所有待处理的样本的文本文件名称（无后缀）

#################
## 程序参数处理##
#################
&GetOptions( 'file_list=s' => \$file_list		 			 
			 );
			 
unless ($file_list) {#这个参数通过输入得到
die $usage;}

################################
##   设置所有目录和文件的路径 ##
################################
our $WORKING_DIR=cwd();#工作目录就是当前目录
our $BIN_DIR=$WORKING_DIR."/tools";#所有可执行文件所在的目录

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
	process_cmd ("$BIN_DIR/fastq-dump $_.sra");		
}
close(IN1);
print "###############################\n";
print "All the samples have been processed by $0\n";

#################
##    子程序   ##
#################
sub process_cmd {
	my ($cmd) = @_;	
	print "CMD: $cmd\n";
	my $ret = system($cmd);	#成功就返回0，否则就失败
	if ($ret) {
		die "Error, cmd: $cmd died with ret $ret";
	}
	return($ret);
}