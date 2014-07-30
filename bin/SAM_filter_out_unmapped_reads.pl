#!/usr/bin/env perl

use strict;
use FindBin;
use lib ("$FindBin::Bin/PerlLib");
use SAM_reader;
use SAM_entry;

my $usage = "usage: SAM_filter_out_unmapped_reads.pl <file.sam> out.unmapped out.mapped > result.sam\n";

my $sam_file	= $ARGV[0] or die $usage;	# input sam file; required
my $outputfile 	= $ARGV[1];			# output unmapped reads; fastq format
my $outputfile2 = $ARGV[2];			# output mapped reads; fastq format

open(OUT,">$outputfile");
open(OUT2,">$outputfile2");
main: {

	my $sam_reader = new SAM_reader($sam_file);
        my $filtered_count = 0;			# count of unmapped read for input sam file
        my $total_count = 0;			# count of read for input sam file 
           
	while ($sam_reader->has_next()) 
	{
		my $sam_entry = $sam_reader->get_next();
		$total_count++;       
	        if ($sam_entry->is_query_unmapped()) {
				$filtered_count++;
				my @current_entry=$sam_entry->get_fields();
				print OUT "@".$current_entry[0]."\n".$current_entry[9]."\n"."+"."\n".$current_entry[10]."\n";	 
		}
		else {
			if($outputfile2){	# output mapped reads to file 
				my @current_entry=$sam_entry->get_fields();
				print OUT2 "@".$current_entry[0]."\n".$current_entry[9]."\n"."+"."\n".$current_entry[10]."\n";
			}
			print $sam_entry->toString() . "\n";#同时输出该sam记录到标准输出
		}
	}
    
	# report unmapped reads for BWA
	print STDERR "this program filtered $filtered_count out of $total_count reads (" . sprintf("%.2f", $filtered_count / $total_count * 100) . ") as unmapped reads, only for BWA\n";
    	exit(0);
}
close(OUT);
close(OUT2);
