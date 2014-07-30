#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
my $usage = <<_EOUSAGE_;

#########################################################################################
# extractFromFasta.pl --inputfile <FILE> --type <String> --query <String> --output1 <FILE> --output2 <FILE>
# Required:
#  --inputfile        a txt file containing a list of input sample names
#  --type             "single" or "list"
#  --query            fasta record list
#  --output1          a file to contain all the fasta record in query
#
###########################################################################################

_EOUSAGE_

	;


my $fasta_file;
my $type;
my $query;
my $output1;
my $output2;

&GetOptions( 'inputfile=s' => \$fasta_file,#
             'type=s' => \$type,
             'query=s' => \$query,
	     'output1=s' => \$output1,
	     'output2=s' => \$output2,
			 );

unless ($fasta_file&&$type&&$query) {
	die $usage;
}



my %list=();
if ($type eq "list") {
	open(IN, $query);
	while(<IN>) {
		chomp($_);
		$list{$_} = 1;
	}
	close(IN);
}
elsif($type eq "single") { 
	$list{$query}=1;
}
else {
	die &usage();
}


open(IN, $fasta_file);
open(OUT1, ">$output1");
open(OUT2, ">$output2");
my $flag = "off";
while(<IN>) {
	if($_ =~ /^>(\S+)/) {
		my $head = $1;
		if(defined $list{$head}) {
			print OUT1 $_;
			$flag = "on";
		}
		else {
			if($type eq "single" && $flag eq "on") {
				exit;
			}
			print OUT2 $_;
			$flag = "off";
		}
	}
	else {
		if($flag eq "on") {
			print OUT1 $_;
		}
		else {
			print OUT2 $_;
		}
	}
}
close(IN);
close(OUT1);
close(OUT2);