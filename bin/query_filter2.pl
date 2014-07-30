#!/usr/bin/perl -w 

# annotation of this script (by kentnf)
# 1. find the best hit (virus) for each query (contigs from sRNA), than output it.
# 2. if the best hit has more than one record, output them
# 3. if the best hit not list in the input1 file, output it

# the input and output file: 
# input:  $sample.known.identified   
# input:  $sample.known.table
# output: $sample.contigs.table

# for format of know.identified:
# HQ593108        7458    96      0.0113971574148565      CONTIG282       1       1.31764705882353        Badnavirus      Banana streak UI virus, complete genome.

# the format of known table (output of blast_parse_talble2.pl);
# query_name \t query_length \t hit_name \t hit_length \t hsp_length \t identity \t evalue \t score \t strand
# \t query_start \t query_end \t hit_start \t hit_end \n
# * the query_name [0], query_length [1], query_start [9], query_end [10] are required 

use strict;

if (@ARGV < 2)
{
	print "usage: query_filter2.pl inputfile1 inputfile2 output1\n";
	exit(0);
}

my $input1 = $ARGV[0]; # input as index
my $input2 = $ARGV[1]; # input as blast table
my $output = $ARGV[2]; # output

# put the refID to hash
# key: refID (The reference is virus reference sequence)
# value: 1
my %hit_index; 
open(IN1, "$input1") || die $!;
while (<IN1>) {
	chomp; 
	my @cols = split(/\t/, $_); 
	defined $hit_index{$cols[0]} or $hit_index{$cols[0]} = 1;
}
close(IN1);

# main
my $last_query="";
my $current_query;

my $high_evalue;	# high evalue for each query. (Should Be The first one)
my $current_evalue;

my $high_identity;	# high identity for each query
my $current_identity;

my $high_record;	# best record for each query
my $if_has_ref=0;	# reference status

open(IN2, "$input2") || die $!;
open(OUT, ">$output") || die $!;

while (<IN2>) {
	my @cols = split(/\t/, $_);
	$current_query=$cols[0];
	$current_identity=$cols[5];
	$current_evalue=$cols[6];

	if($current_query ne $last_query)		# new query record
	{
		# parse pre query
		# if the previous query (contig) do not have reference, ouput the best one.
		if( $if_has_ref==0 && $last_query ne "") {
			print OUT $high_record;
		}

		# init this vars for new record
		$if_has_ref=0;  

		# check if the query (contigs) could aligned to ref (virus)
		# then ouput the first result (should be best one according to evalue)
		if (defined($hit_index{$cols[2]}))
		{
			$high_identity=$cols[5];
			$high_evalue=$cols[6];
			print OUT $_;
			$if_has_ref=1;		
		}
		else
		{
			$high_record = $_; # the best record for each query, output it when meet next query			
		}
	}
	else						# same as pre-query
	{
		if($if_has_ref==1) 
		{
			# output result better than/equal to the first one
			# * that is the reason why two result with same evalue could be output
			#if (defined($hit_index{$cols[2]}) && $current_evalue >= $high_evalue && $current_identity >= $high_identity)
			#{
			#	print OUT $_;
			#}

			# output all the result (changed by kentnf)
			if (defined ($hit_index{$cols[2]})) { print OUT $_; }
		}
		else
		{	
			# check if the query (contigs) could aligned to ref (virus)
			# then ouput the first result (should be best one according to evalue)
			if (defined($hit_index{$cols[2]}))
			{
				$high_identity=$cols[5];
				$high_evalue=$cols[6];		
				print OUT $_;
				$if_has_ref=1;	# yes, the contig has a reference
			}
		}	
	}
	$last_query=$cols[0];
}

# changed by kentnf
# print the last one result if best hit is not list in input1 file
if ($if_has_ref == 0) {
	print OUT $high_record;
}

close(IN2);
close(OUT);

