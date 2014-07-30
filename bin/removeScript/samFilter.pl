#!/usr/bin/perl -w 
#
use strict; 

-t and !@ARGV and die "perl $0 2.sam\n"; 

our $max_distance = 2;			# set $max_distance for all selected hits
my @lls; 
my $pre_query_name; 
my $nmv; 				# the minimum edit distance for each query
$pre_query_name = ''; 
$nmv = -1; 
my ($query_col, $opt_col) = (0, 11); 	# query and option column number for sam
while (<>) {
	chomp; s/\s+$//; 
	my @ta = split(/\t/, $_);				# split one line to array ta
	my $query_name = $ta[$query_col];			# get query read name

	if ($pre_query_name eq '') {				# check pre query name
		$pre_query_name = $query_name; 
	}elsif ($query_name ne $pre_query_name) {		# 
		$nmv < $max_distance and print STDOUT join("\n", @lls)."\n";#nmv必须小于$max_distance，才保留，注意bwa有时会输出超过预定$max_distance的hits
		$nmv=-1; 
		$pre_query_name = $query_name; 
	}

	for my $each_opt (@ta[$opt_col .. $#ta]) {		# scan for column from 11 to end
		if ($each_opt =~ /^NM:i:(\d+)$/) {		# 'NM:i:X'
			my $distance = $1; 			# gent the distance between query and reference
			if ($nmv == -1 or $distance < $nmv) {	# find the minimum edit distance for each query
				$nmv = $distance;
				@lls = ($_); 
			}elsif ($distance == $nmv) {
				push(@lls, $_); 
			}else{
				; 
			}
			last; 					# each line must have one 'NM:i:<X>'
		}
	}
}

if ($nmv > -1) {
	$nmv < $max_distance and print STDOUT join("\n", @lls)."\n"; # for MAX NM:i:? 
	$nmv = -1; 
	$pre_query_name = ''; 
}
