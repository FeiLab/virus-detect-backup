#!/usr/bin/perl -w 
use strict;
# 首先调换列次数，然后排序，最后加表头输出
# 下列column需要互相调换,先根据decription，后根据Hit_ID排序
# Contig_ID	Contig_length	Hit_ID	Hit_length	strand	Contig_start	Contig_end	Hit_start	Hit_end	identity	genus	 description	contig_sequence

if (@ARGV < 2)
{
  print "usage: arrange_col.pl input output\n";
  exit(0);
}

our $input = $ARGV[0]; #这个是作为index的文件
our $output = $ARGV[1]; #这是要得到的输出文件，包括了所有提取的序列

my @all_data;
open(IN, "$input");
while (<IN>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	push(@all_data, [@ta[0,13,1,2,3,11,12,6,7,8,9,10,4,5]]);#contig数量，hit长度和整行 
}
close(IN);
@all_data = sort { ($a->[6] cmp $b->[6]) || ($a->[4] cmp $b->[4])} @all_data; #排序
open(OUT, ">$output");
	print OUT "Contig_ID\tContig_Seq\tContig_Len\tHit_ID\tHit_Len\tGenus\tDescription\tContig_start\tContig_end\tHit_start\tHit_end\tHsp_identity\tE_value\tHsp_strand\n";#就输出到标准输出
	foreach my $each (@all_data){
		print OUT join("\t", @$each)."\n";#就输出到标准输出	
	}
close(OUT);
