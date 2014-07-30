#!/usr/bin/perl -w 
use strict; 
# usage: ./blast_fileter2.pl inputfile outputfile
# 来自blast_parse_table2.pl程序的输出
# query_name	query_length	hit_name	hit_length	hsp_length	identity	evalue	score	strand	query_start	query_end	hit_start	hit_end

# 只保留每个<query,hit>中 evalue最低的(即第一个)

if (@ARGV < 2)
{
  print "usage: hit_cov.pl input output1 query_coverage\n";
  exit(0);
}

our $input = $ARGV[0]; #去掉文件后缀名称
our $output = $ARGV[1];

open(IN, "$input");
open(OUT, ">$output");

my $last_query="";
my $last_hit="";
while (<IN>) {
	my @cols = split(/\t/, $_); 
	if($cols[0] ne $last_query || $cols[2] ne $last_hit){#如果当前的与上一个不一致
		print OUT $_;
	}
	$last_query=$cols[0];
	$last_hit=$cols[2];
}
close(IN);
close(OUT);


