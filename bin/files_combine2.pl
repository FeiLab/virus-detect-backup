#!/usr/bin/perl
#use strict;   #变量必须定义才可以使用
#use warnings; #对未初始化变量是给出警告
use Getopt::Long; #这个模块可以接受完整参数

#参数处理
my $usage = <<_EOUSAGE_;

#########################################################################################
# files_combine2.pl --filelist <FILE> --folder1 <String> --folder2 <String>
# Required:
#  --filelist       a fastq file name without suffix for processing 

##########################################################################################

_EOUSAGE_
	;
#定义变量接受命令行参数	
our $filelist;        #包括需要处理的所有fastq文件名称列表
our $folder1;  
our $folder2; 
our $output_folder="combined";

#从命令行参数向变量传值
&GetOptions( 'filelist=s' => \$filelist,
	     'folder1=s' => \$folder1,
	     'folder2=s' => \$folder2
			 );

#主程序开始
open(IN1,$filelist) || die "Can't open the $filelist file\n";
while(<IN1>){
	chomp;
	$sample=$_; #每次循环读入一行，后续代码都是处理该样本文件（名称无后缀）。
	$file1="./".$folder1."/".$sample.".".$folder1.".fa";
	$file2="./".$folder2."/".$sample.".".$folder2.".fa";
	my $result_file = "./".$output_folder."/".$sample.".contigs1.fa";
	&process_cmd("cat $file1 $file2 > $result_file");
}
close(IN1);
####
sub process_cmd {
	my ($cmd) = @_;	
	print "CMD: $cmd\n";
	my $ret = system($cmd);	
	if ($ret) {#如果不为0，表示不成功
		print "Error, cmd: $cmd died with ret $ret";
	}
	return($ret);
}
