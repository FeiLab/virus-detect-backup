#!/usr/bin/perl
#use strict;
#use warnings;
use Getopt::Long;

# usage info
my $usage = <<_EOUSAGE_;

#########################################################################################
# files_combine2.pl --filelist <FILE> --folder1 <String> --folder2 <String>
# Required:
#  --filelist       a fastq file name without suffix for processing 

##########################################################################################

_EOUSAGE_
;

# define parameters	
our $filelist;
our $folder1;  
our $folder2; 
our $output_folder="combined";

# get input parameters
GetOptions( 
	'filelist=s' => \$filelist,
	'folder1=s' => \$folder1,
	'folder2=s' => \$folder2
);

# main
open(IN1,$filelist) || die "Can't open the $filelist file\n";
while(<IN1>){
	chomp;
	$sample = $_;
	$file1	= "./".$folder1."/".$sample.".".$folder1.".fa";
	$file2	= "./".$folder2."/".$sample.".".$folder2.".fa";
	my $result_file = "./".$output_folder."/".$sample.".contigs1.fa";
	process_cmd("cat $file1 $file2 > $result_file");
}
close(IN1);

# subroutine
sub process_cmd {
	my ($cmd) = @_;	
	print "CMD: $cmd\n";
	my $ret = system($cmd);	
	if ($ret) {
		print "Error, cmd: $cmd died with ret $ret";
	}
	return($ret);
}
