#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;
use FindBin;
use File::Basename;
my $usage = <<_EOUSAGE_;

#########################################################################################
# bwa_remove.pl --file_list <FILE> --reference <FILE>
#		--max_dist[INT] --max_open [INT] --max_extension [INT] --len_seed [INT] --dist_seed [INT] --thread_num [INT]
#
# Required(2):
#  --file_list	The name of a txt file containing a list of input file names without any suffix
#  --reference  a fasta file containing the host reference genome or transcriptom
#
# BWA-related options(6):
#  --max_dist      Maximum edit distance [1]  
#  --max_open      Maximum number of gap opens [1]  
#  --max_extension Maximum number of gap extensions [1]  
#  --len_seed      Take the first INT subsequence as seed [15] 
#  --dist_seed     Maximum edit distance in the seed [1]  
#  --thread_num    Number of threads (multi-threading mode) [8]  
###########################################################################################

_EOUSAGE_
;

#################
# global vars   #
#################
our $file_list;
our $reference;
our $index_name;
our $max_dist = 1;	# bwa allowed max edit distance 
our $max_open = 1;	# bwa allowed max gap number
our $max_extension = 1;	# bwa allowed max gap length, -1 means disallow long gap
our $len_seed = 15;	# bwa seed length
our $dist_seed = 1;	# bwa seed allowed max edit distance
our $thread_num = 8;
 
################################
# set folder and file path     #
################################
our $WORKING_DIR	= cwd();				# current folder
our $DATABASE_DIR 	= ${FindBin::RealBin}."/../databases";	# database folder
our $BIN_DIR 		= ${FindBin::RealBin};			# programs folder

##################
# get input para #
##################
GetOptions( 
	'file_list=s'	=> \$file_list,	# sample list
	'reference=s'	=> \$reference,	# reference file (fasta format)
	'max_dist=i'	=> \$max_dist,
	'max_open=i'	=> \$max_open,
	'max_extension=i' => \$max_extension,
	'len_seed=i' 	=> \$len_seed,
	'dist_seed=i' 	=> \$dist_seed,			 
	'thread_num=i' 	=> \$thread_num
);

die $usage unless ($file_list && $reference);	# required parameters

$index_name = basename($reference);		# remove the path 
#$index_name =~ s/\.\S*$//;			# remove the suffix

# main
main: {

	# creat bwa index for reference
	process_cmd("$BIN_DIR/bwa index -p $DATABASE_DIR/$index_name -a bwtsw $reference 2> bwa.log") unless (-e "$DATABASE_DIR/$index_name.amb");

	my $sample;
	my $i=0;
        open(IN, "$file_list") || die $!;
        while (<IN>) {
		$i=$i+1;
		chomp;
		$sample=$_;
		print "#processing sample $i by $0: $sample\n";

		# command lines		
		process_cmd("$BIN_DIR/bwa aln -n $max_dist -o $max_open -e $max_extension -i 0 -l $len_seed -k $dist_seed -t $thread_num $DATABASE_DIR/$index_name $sample 1> $sample.sai 2>bwa.log") unless (-s "$sample.sai");
		process_cmd("$BIN_DIR/bwa samse -n 1 $DATABASE_DIR/$index_name $sample.sai $sample 1> $sample.pre.sam 2>bwa.log") unless (-s "$sample.pre.sam");				
		process_cmd("$BIN_DIR/SAM_filter_out_unmapped_reads.pl $sample.pre.sam $sample.unmapped $sample.mapped > $sample.sam") unless (-s "$sample.sam");
		system("rm $sample.sai");
		system("rm $sample.pre.sam");
		system("rm $sample.sam");
		system("rm $sample.mapped");
	}
        close(IN);
	print "###############################\n";
	print "All the input files have been processed by $0\n";
	#system("touch $index_name.remove.run.finished");	# create finished sign
}

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
