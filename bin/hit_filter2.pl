#!/usr/bin/perl -w 

use strict; 
use IO::File; 
use FindBin;
use Getopt::Long;
use Cwd;
my $usage = <<_EOUSAGE_;

#########################################################################################
# hit_filter2.pl --input1 <FILE> input2 <FILE> --diff_ratio <FLOAT> 
# 		 --diff_contig_cover <FLOAT> --diff_contig_length <INT> --output <FILE>
#
# Required[3]:
#  --input1	A name of an input file sample.known.identified1
#  --input2	A name of an input file sample.known.table
#  --output	A name of a output file sample.known.identified
#
# options[3]:
#  --diff_ratio	[0.25]
#  --diff_contig_cover  [0.5]
#  --diff_contig_length [100]
#
##########################################################################################

_EOUSAGE_
	;
#################
# set vars  	#
#################	
my ($input1, $input2, $diff_ratio, $diff_contig_cover, $diff_contig_length, $output);            
$diff_ratio = 0.25;
$diff_contig_cover = 0.5;
$diff_contig_length = 100;

################################
# set folder and file path     #
################################
my $WORKING_DIR = cwd();			# current folder is working folder
my $BIN_DIR = ${FindBin::RealBin}."/bin";	# script folder

#########################
# input parameter	#
#########################
GetOptions( 
	'input1=s' 		=> \$input1,
	'input2=s'		=> \$input2,
	'diff_ratio=f' 		=> \$diff_ratio,
	'diff_contig_cover=f'	=> \$diff_contig_cover,
	'diff_contig_length=i'	=> \$diff_contig_length,
	'output=s' 		=> \$output
);

die $usage unless $input1 && $input2 && $output;# required parameter
			 
#################
##   main      ##
#################

# load all data to array: @all_data
# structure:
# [contig num, covered bp, line],[contig num, covered bp, line], ... ,[contig num, covered bp, line]
my (@all_data, $name, $sequence); 
$name = ''; $sequence = ''; 

open(IN, "$input1") || die "Can not open input file $input1 $!\n";
while (<IN>) {
	chomp; 
	my @each_line = split(/\t/, $_);	
    	my @contigs= split(/,/, $each_line[4]);	# column 5 has all contig name for each hit
    	my $covered_bp = $each_line[7];
    	if ( $each_line[7] eq "") { $covered_bp = 0; }	
	push(@all_data, [scalar(@contigs), $covered_bp, $_]);
}
close(IN);

# sort all dataset by 1) number of contigs for each hit
# 		      2) covered bp
# 		      in descending order
@all_data = sort { -1*($a->[0] <=> $b->[0]) || -1*($a->[1] <=> $b->[1])} @all_data;

# load contig hit covered start, end to hash
# hash : contig_position
# key : contig \t hit
# value : hit_start \t hit_end
my %contig_position; 
my $fh = IO::File->new($input2) || die "Can not open input file 2: $input2 $!\n";
while(<$fh>)
{
	chomp;
	next if $_ =~ m/^#/;
	my @a = split(/\t/, $_);
	$contig_position{$a[0]."\t".$a[2]} = "$a[11]\t$a[12]";
}
$fh->close;

=head
for my $tr (@all_data) { print $tr->[2]."\n"; }	# debug for sorted data
=cut

# story non-redundancy and redundancy to array
# @inset : lines of non-redundancy contigs
# @restset : lines of redundancy contigs 
my @inset;
my @restset = ();
my $contig_count = 1; 

foreach my $tr (@all_data) 
{
	if (scalar(@inset)  == 0) {
		# the first line contain most contigs than others
		push(@inset, $tr->[2]);
	}else{
		my @aa = split(/\t/, $tr->[2]);

		# get return string for non-redundancy and redundancy status
		# return_string: n -- non-redundancy 
		# return_string: r -- redundancy
		# put non-redundancy and redundancy data to array inset, restset
		my $return_string = ifRedundant(\@inset, \$aa[4], $aa[0], $aa[1]);
		if ($return_string eq "n") {	
			push(@inset, $tr->[2]); 
		}else{
			push(@restset, $tr->[2]); 
		}
	}
}

# output all non-redundancy data
open(OUT1, ">$output") || die "Can not open output file $output $!";
foreach my $tr (@inset) { print OUT1 $tr."\n"; }
close(OUT1);

=head2 output the redundancy data for debug
open(OUT2, ">restset");
for my $tr (@restset) { print OUT2 $tr."\n"; }
close(OUT2);
=cut

#################
# subroutine	#
#################
sub ifRedundant 
{
	my ($inset, $query, $hit_name, $hit_length) = @_; 

	my @query_contigs = split(/,/, $$query);
	my $total_contigs = scalar(@query_contigs);
	my $ratio;

	foreach my $tr (@$inset) 
	{
		# put contigs to hash for non-redundancy data 
		my @aa = split(/\t/, $tr);
		my @contigs= split(/,/, $aa[4]);
		my %contigs;
		foreach my $each_contig ( @contigs ) { $contigs{$each_contig} = 1; }

		# identify the new non-redundancy data
		my $diff_contigs = 0;

		my @hit_range;
		foreach my $each_contig ( @query_contigs )
		{
			unless (defined $contigs{$each_contig})
			{
				$diff_contigs++;
				if ( defined $contig_position{$each_contig."\t".$hit_name} )
				{
					push(@hit_range, $contig_position{$each_contig."\t".$hit_name});
				}
				else
				{
					die "Error for pair of $each_contig $hit_name\n";
				}		
			}
		}

		$ratio = $diff_contigs * 1.0 / $total_contigs;
		my $length_of_diff_contigs_cover = find_length(@hit_range);
		my $cover_of_diff_contigs = $length_of_diff_contigs_cover / $hit_length;
		
		if( $ratio <= $diff_ratio && $cover_of_diff_contigs <= $diff_contig_cover && $length_of_diff_contigs_cover <= $diff_contig_length ) 
		{ return "r"; } 
	}
	
	return "n";
}

=head2 find_length 

=cut
sub find_length
{
	my @hit_range = @_;
	my %pos;
	foreach my $range (@hit_range)
	{
		my ($start, $end) = split(/\t/, $range);
		for(my $i=$start; $i<=$end; $i++)
		{
			$pos{$i}++;
		}
	}
	my $len = scalar(keys(%pos));
	return $len;
}

sub process_cmd {
	my ($cmd) = @_;	
	print "CMD: $cmd\n";
	my $ret = system($cmd);	
	if ($ret) {
		die "Error, cmd: $cmd died with ret $ret";
	}
	return($ret);
}
