#!/usr/bin/perl -w 
use strict;
# 首先调换列次数，然后排序，最后加表头输出
# 下列column需要互相调换,先根据Genus，后根据read_cov(bp)排序


if (@ARGV < 1)
{
  print "usage: arrange_col2.pl input1 input2 > output\n";
  exit(0);
}

our $input1 = $ARGV[0]; #输入文件
our $input2 = $ARGV[1]; #输入文件，包括数据的总read数

my @all_data;
open(IN2, "$input2");
my $first_line =<IN2>;
my ($total_reads)=($first_line=~m/# reads processed: (\d+)/);
close(IN2);

open(IN1, "$input1");
while (<IN1>) {
	chomp; 
	my @ta = split(/\t/, $_);
    my $coverage= 1.0*$ta[7]/$ta[1];	
	my $depth_norm= 1e6*$ta[6]*$ta[7]/($ta[1]*$total_reads); # 把depth编程标准化后的depth
	push(@all_data, [@ta[0,1,7],$coverage,@ta[4,5],$depth_norm,@ta[8,9,3]]);#选择需要的列，并重新排列顺序 
}
close(IN1);

@all_data = sort { ($a->[7] cmp $b->[7]) || ($b->[2] <=> $a->[2])} @all_data; #根据Genus和read_cov(bp)排序

foreach my $each (@all_data){
	print join("\t", @$each)."\n";#就输出到标准输出	
}
