#!/usr/bin/perl -w 
use strict;
# 首先调换列次数，然后排序，最后加表头输出
# 下列column需要互相调换,先根据decription，后根据Hit_ID排序
# Contig_ID	Contig_length	Hit_ID	Hit_length	strand	Contig_start	Contig_end	Hit_start	Hit_end	identity	genus	 description	contig_sequence

if (@ARGV < 1)
{
  print "usage: arrange_col2.pl input > output\n";
  exit(0);
}

our $input = $ARGV[0]; #输入文件

my @all_data;
open(IN, "$input");
while (<IN>) {
	chomp; 
	my @ta = split(/\t/, $_);
    my $coverage= 1.0*$ta[7]/$ta[1];	
	push(@all_data, [@ta[0,1,7],$coverage,@ta[4,5,6,8,9]]);#选择需要的列，并重新排列顺序 
}
close(IN);
@all_data = sort { ($a->[7] cmp $b->[7]) || ($b->[2] <=> $a->[2])} @all_data; #根据Genus和hit_covered(bp)排序

foreach my $each (@all_data){
	print join("\t", @$each)."\n";#就输出到标准输出	
}
