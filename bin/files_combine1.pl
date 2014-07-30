#!/usr/bin/perl
#use strict;   #变量必须定义才可以使用
#use warnings; #对未初始化变量是给出警告
use Getopt::Long; #这个模块可以接受完整参数

#参数处理
my $usage = <<_EOUSAGE_;

#########################################################################################
# files_combine.pl --filelist <FILE> --suffix <String> --startNumber <INT> --endNumber <INT>
# Required:
#  --filelist       a fastq file name without suffix for processing 
#  --suffix
#  --number
##########################################################################################

_EOUSAGE_
	;
#定义变量接受命令行参数	
our $filelist;        #包括需要处理的所有fastq文件名称列表
our $suffix;        #
our $startNumber = 1;  
our $endNumber; 

#从命令行参数向变量传值
&GetOptions( 'filelist=s' => \$filelist,
			'suffix=s' => \$suffix,
			'startNumber=i' => \$startNumber,
			'endNumber=i' => \$endNumber
			 );

#主程序开始
open(IN1,$filelist) || die "Can't open the $filelist file\n";
my $files1="";
while(<IN1>){
	chomp;
	my $files="";
	
	for(my $i = $startNumber; $i<= $endNumber; $i++){
		$files=$files.$_.".$suffix$i ";#所有需要合并的文件名称连在一起
	}
	
	my $finalFile=$_.".$suffix";
	print "cat $files > $finalFile\n";
	my $result = `cat $files> $finalFile`;		
}
close(IN1);

