#!/usr/bin/perl
use strict;
use warnings;
if (@ARGV < 2)
{
  print "usage: fasta2tab input output\n";
  exit(0);
}
our $input = $ARGV[0]; #
our $output = $ARGV[1]; #


#读fasta file文件，
open(IN, "$input");
open(OUT, ">$output");
my $last_header = "";
my $last_sequence;
while(<IN>) {
	if($_ =~ m/^>/) {
		if($last_header ne ""){			
			print OUT $last_header."\t".$last_sequence."\n";
		}
		$last_header = $_;#输出后立刻更新为当前的header
		chomp($last_header);
		$last_header=~ s/>//;
		$last_sequence = "";#同时，开始记录一条新序列
	}
	else {
		chomp;
		$last_sequence .= $_;
	}
}
print OUT $last_header."\t".$last_sequence."\n"; #一定不要忘记输出最后一个记录
close(IN);
close(OUT);
