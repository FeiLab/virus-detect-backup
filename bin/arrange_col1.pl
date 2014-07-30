#!/usr/bin/perl -w 
use strict;
# 首先调换列次数，然后排序，最后加表头输出
# 下列column需要互相调换,先根据decription，后根据Hit_ID排序
# Contig_ID	Contig_length	Hit_ID	Hit_length	strand	Contig_start	Contig_end	Hit_start	Hit_end	identity	genus	 description	contig_sequence

if (@ARGV < 1)
{
  print "usage: arrange_col1.pl input > output\n";
  exit(0);
}

our $input = $ARGV[0]; #输入文件

my @all_data;
open(IN, "$input");
while (<IN>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	push(@all_data, [@ta[0,19,1,2,3,17,18,9,10,11,12,13,6,8]]);#选择需要的列，并重新排列顺序 
}
close(IN);
@all_data = sort { ($a->[5] cmp $b->[5]) || ($a->[3] cmp $b->[3])} @all_data; #根据Genus和hit排序

print "Contig_ID\tContig_Seq\tContig_Len\tHit_ID\tHit_Len\tGenus\tDescription\tContig_start\tContig_end\tHit_start\tHit_end\tHsp_identity\tE_value\tHsp_strand\n";#就输出到标准输出
foreach my $each (@all_data){
	print join("\t", @$each)."\n";#就输出到标准输出	
}
