#!/usr/bin/perl -w 
use strict; 

# 来自blast_parse_table2.pl程序的输出
# query_name	query_length	hit_name	hit_length	hsp_length	identity	evalue	score	strand	query_start	query_end	hit_start	hit_end
# query或hit至少要被hsp覆盖一定比例
# 相同的<query,hit>中，只保留evalue最低的(即第一个)
# 一个query对应的hits中 只保留比最低evalue那个hsp的identity少于drop_off(5%)以内的
if (@ARGV < 2)
{
  print "usage: blast_filter5.pl input output coverage identity_dropoff\n";
  exit(0);
}

our $input = $ARGV[0]; #去掉文件后缀名称
our $output = $ARGV[1];
our $coverage=$ARGV[2];
our $identity_dropoff = $ARGV[3];

open(IN, "$input");
open(OUT, ">$output");

my $last_query="";
my $last_hit="";
my $current_query;
my $current_hit;

my $high_identity;#保存每个query最高的identity
my $current_identity;
while (<IN>) {
	my @cols = split(/\t/, $_);
	if(3.0*$cols[4]/$cols[1] >= $coverage || 3.0*$cols[4]/$cols[3] >= $coverage){#query或hit至少要被hsp覆盖一定比例
		$current_query=$cols[0];
		$current_hit=$cols[2];	
		$current_identity=$cols[5];
		if($current_query ne $last_query){#如果遇到一个新的query
			$high_identity=$cols[5];#记录最高的identity
			print OUT $_;#query的第一个hit肯定要，因为是evalue最低的
		}else{
			if ($current_hit ne $last_hit && ($current_identity >= $high_identity-$identity_dropoff))#相同的<query,hit>中，只保留evalue最低的(即第一个)
			{
				print OUT $_;	
			}	
		}
		$last_query=$cols[0];
		$last_hit=$cols[2];
	}
}
close(IN);
close(OUT);


