#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;
use FindBin;
use lib "$FindBin::RealBin/bin/PerlLib";
use IO::File;
use Util;

my $usage = <<_EOUSAGE_;

#########################################################################################
# virus_detect.pl --file_type [String] --reference [FILE] --coverage [Float]
#                 --host_removal --host_reference [FILE] --objective_type [String]
#                 --max_dist[INT] --max_open [INT] --max_extension [INT] --len_seed [INT] --dist_seed [INT] --thread_num [INT]
#                 --strand_specific --min_overlap [INT] --max_end_clip [INT] --cpu_num [INT] --mis_penalty [INT] --gap_cost [INT] --gap_extension [INT]
#		  input_file1 input_file2 ... input_fileN
#  
# Options(6):
#  --file_type The format of input files(fastq or fasta)  [fastq]
#  --reference The name of a fasta file containing all of the virus reference sequences  [vrl_genbank.fasta]
#  --coverage  The reference sequence with at least a pencetage of its length covered by reads is detected as a aligned virus contig  [0.3]
#  --host_removal         Remove host sequences before assembly or not  [Not selected]
#  --host_reference       If using host_removal, the name of a host reference file should be provided  [Null]
#  --objective_type       The objective type of optimization(maxLen, n50 or avgLen)  [maxLen]
# 
# BWA-related options(6):
#  --max_dist      Maximum edit distance [1]  
#  --max_open      Maximum number of gap opens [1]  
#  --max_extension Maximum number of gap extensions [1]  
#  --len_seed      Take the first INT subsequence as seed [15] 
#  --dist_seed     Maximum edit distance in the seed [1]  
#  --thread_num    Number of threads (multi-threading mode) [8] 
# 
# Megablast-related options(7): for virus detection
#  --strand_specific Only for sequences assembled from strand specific RNA-seq  [Not selected]
#  --min_overlap     The minimum overlap length between two contigs to be combined [30]
#  --max_end_clip    The maximum length of end clips [6]
#  --cpu_num         Number of processors to use [8] 
#  --mis_penalty     Penalty for a nucleotide mismatch [-1]
#  --gap_cost        Cost to open a gap [2] 
#  --gap_extension   Cost to extend a gap [1] 
#
# Megablast-related options(): for virus identification
#  --word_size	     [11] 
#  --exp_value	     [1e-5]
#  --identity_percen [80] 
#  --mis_penalty_b   Penalty for a nucleotide mismatch [-1]
#  --gap_cost_b      Cost to open a gap [2] 
#  --gap_extension_b Cost to extend a gap [1]
#
# Result filter options(4):
#  --hsp_cover		The coverage of hsp should be more than this cutoff for query or hit [0.75]
#  --diff_ratio 	The hits with distance less than 0.25 will be combined into one  [0.25]
#  --diff_contig_cover	The coverage for different contigs [0.5]
#  --diff_contig_length	The length of different contigs [100] 
#
###########################################################################################

_EOUSAGE_
;

################################
##  set file folder path      ##
################################
my $WORKING_DIR   = cwd();				# set current folder as working folder
my $DATABASE_DIR  = ${FindBin::RealBin}."/databases";	# set database folder
my $BIN_DIR	  = ${FindBin::RealBin}."/bin";		# set script folder 
my $TEMP_DIR	  = $WORKING_DIR."/temp";		# set temp folder

###############################
##   global vars	     ##
###############################
my $file_type= "fastq";				# input file type, fasta or fastq
my $reference= "vrl_genbank.fasta";		# virus sequence
my $coverage=0.3;  				# 每条参考序列如果被reads覆盖的部分占全长比例的阈值
my $host_removal;         			# switch for host removing
my $host_reference;       			# host reference
my $objective_type='maxLen';			# objective type for Velvet assembler: n50、maxLen, avgLen

# paras for BWA
my $max_dist = 1;  				# max edit distance 
my $max_open = 1;  				# max gap opening
my $max_extension = 1; 				# max gap extension (gap length)
my $len_seed = 15; 				# bwa seed length
my $dist_seed = 1; 				# bwa seed max edit distance
my $thread_num = 8; 				# bwa thread number

# paras for megablast detection (remove redundancy )
my $strand_specific;  				# switch for strand specific transcriptome data? 
my $min_overlap = 30; 				# minimum overlap for hsp combine
my $max_end_clip = 6; 				# max end clip for hsp combine
my $cpu_num = 8;          			# thread number
my $mis_penalty = -1;     			# megablast mismatch penlty, minus integer
my $gap_cost = 2;         			# megablast gap open cost, plus integer
my $gap_extension = 1;    			# megablast gap extension cost, plus integer 

# paras for blast && identification 
my $word_size = 11;
my $exp_value = 1e-5;				#
my $identity_percen = 25;			# tblastx 以蛋白质序列来比对时hsp的最小同一性
my $mis_penalty_b = -1;				# megablast mismatch penlty, minus integer
my $gap_cost_b = 2;				# megablast gap open cost, plus integer
my $gap_extension_b = 1;			# megablast gap extension cost, plus integer

my $filter_query = "F";				# megablast switch for remove simple sequence
my $hits_return = 500;				# megablast number of hit returns

# paras for result filter
my $hsp_cover = 0.75;
my $diff_ratio= 0.25;
my $diff_contig_cover = 0.5;
my $diff_contig_length= 100; 

# disabled parameters
my $input_suffix='clean'; 			# input_suffix, disabled

# get input paras #
GetOptions(
	'file_type=s'	=> 	\$file_type,
	'reference=s'	=> 	\$reference,
	'coverage=f'	=> 	\$coverage,
	'host_removal!' => 	\$host_removal,
	'host_reference=s' => 	\$host_reference,
	'objective_type=s' => 	\$objective_type,
		
	'max_dist=i' => 	\$max_dist,
	'max_open=i' => 	\$max_open,
	'max_extension=i' => 	\$max_extension,
	'len_seed=i' => 	\$len_seed,
	'dist_seed=i' => 	\$dist_seed,			 
	'thread_num=i' => 	\$thread_num,
		
	'strand_specific!' => 	\$strand_specific,
	'min_overlap=i' => 	\$min_overlap,
	'max_end_clip=i' => 	\$max_end_clip,
	'cpu_num=i' => 		\$cpu_num ,
	'mis_penalty=i' => 	\$mis_penalty,
	'gap_cost=i' => 	\$gap_cost,
	'gap_extension=i' => 	\$gap_extension,
	
	'word_size=i' =>  	\$word_size,
	'exp_value-s' =>  	\$exp_value,
	'identity_percen=s' => 	\$identity_percen,	# tblastx以蛋白质序列来比对时hsp的最小同一性
	'mis_penalty_b=i' => 	\$mis_penalty_b,
	'gap_cost_b=i' => 	\$gap_cost_b,
	'gap_extension_b=i' => 	\$gap_extension_b,

	'hsp_cover=s' =>	\$hsp_cover,
	'diff_ratio=s' => 	\$diff_ratio,
	'diff_contig_cover=s' =>\$diff_contig_cover,
	'diff_contig_length=s'=>\$diff_contig_length
);

# put input file parameters to filelist array
my @filelist;
if ( !@ARGV ) {
	die "\nPlease input files:\n\n$usage\n";
} else {
	foreach my $m (@ARGV) {
		if (-s $m) {
			push(@filelist, $m);
		} else {
			print "Error, can not parse input file $m\n";
		}
	}
}

# set paramsters
my $parameters_remove_redundancy = "--min_overlap $min_overlap --max_end_clip $max_end_clip --cpu_num $cpu_num ".
				   "--mis_penalty $mis_penalty --gap_cost $gap_cost --gap_extension $gap_extension";
if ($strand_specific) { $parameters_remove_redundancy.=" --strand_specific"; }

my $parameters_bwa_align = "--max_dist $max_dist --max_open $max_open --max_extension $max_extension ".
			   "--len_seed $len_seed --dist_seed $dist_seed --thread_num $thread_num";

# main
main: {

	# create temp folder
	# create softlink for input file
	# create temp file list for input file
	# *the file in file list has full path
	Util::process_cmd("mkdir $TEMP_DIR") unless -e $TEMP_DIR;
	my $file_list = "$TEMP_DIR/temp_file_list";
	my $fh_list = IO::File->new(">".$file_list) || "Can not open temp file list $file_list $!\n";
	foreach my $sample (@filelist) { 
		Util::process_cmd("ln -s $WORKING_DIR/$sample $TEMP_DIR/$sample") unless -e "$TEMP_DIR/$sample";
		print $fh_list $TEMP_DIR."/".$sample."\n"; 
	}
	$fh_list->close;

	# detect known virus
	# align small RNA reads to know plant virus sequence 
	print "alignend small RNA to know plant virus database: $reference\n";
	my $cmd_align = "$BIN_DIR/bwa-alignAndCorrect.pl --file_list $file_list --reference $DATABASE_DIR/$reference --coverage $coverage ".
			"$parameters_bwa_align --output_suffix aligned";
	Util::process_cmd($cmd_align);

	# remove redundancy
	my $cmd_removeRedundancy = "$BIN_DIR/removeRedundancy_batch.pl --file_list $file_list --file_type $file_type --input_suffix aligned ".
				   "--contig_prefix KNOWN $parameters_remove_redundancy";
	Util::process_cmd($cmd_removeRedundancy);

	#detect unknown virus, then assembly them 
	print "\nassembly small RNA to contigs\n\n";

	if($host_removal){	
		Util::process_cmd("$BIN_DIR/bwa_remove.pl --file_list $file_list --reference $DATABASE_DIR/$host_reference $parameters_bwa_align");
		# kentnf: this could be writen to a subroutine

		# the input suffix of unmapped reads is 'unmapped'
		# the seq from sam file is fastq format, no matter the format of input file
		Util::process_cmd("$BIN_DIR/Velvet_Optimiser.pl --file_list $file_list --input_suffix unmapped --file_type fastq --objective_type $objective_type --hash_end 19 --coverage_end 25");
		Util::process_cmd("$BIN_DIR/runVelvet.pl --parameters $TEMP_DIR/optimization.result --input_suffix unmapped --file_type fastq --output_suffix assemblied");
	}	
	else
	{
		Util::process_cmd("$BIN_DIR/Velvet_Optimiser.pl --file_list $file_list --file_type $file_type --objective_type $objective_type --hash_end 19 --coverage_end 25");
		Util::process_cmd("$BIN_DIR/runVelvet.pl --parameters $TEMP_DIR/optimization.result --file_type $file_type --output_suffix assemblied");
	}

	# remove redundancy of assembly results
	print "\nremove redundancy for assemblied unknown small RNA\n\n";
	$cmd_removeRedundancy = "$BIN_DIR/removeRedundancy_batch.pl --file_list $file_list --file_type $file_type --input_suffix assemblied ".
				"--contig_prefix NOVEL $parameters_remove_redundancy";
	Util::process_cmd($cmd_removeRedundancy);

	# combine the known and unknown virus
	print "combine the known and unkown contigs\n";
	files_combine2($file_list);

	# remove redundancy of combined results, it must be using strand_specific parameter
	print "remove redundancy\n";
	$cmd_removeRedundancy = "$BIN_DIR/removeRedundancy_batch.pl --file_list $file_list --file_type $file_type --input_suffix combined ".
       				"--contig_prefix CONTIG --strand_specific $parameters_remove_redundancy";
	Util::process_cmd($cmd_removeRedundancy);

	# identify the virus
	my $cmd_identify = "$BIN_DIR/virus_identify.pl ";
       	$cmd_identify .= "--file_list $file_list --file_type $file_type --reference $reference --contig_type combined ";
	$cmd_identify .= "--word_size $word_size --exp_value $exp_value --identity_percen $identity_percen ";
	$cmd_identify .= "--cpu_num $cpu_num --mis_penalty $mis_penalty_b --gap_cost $gap_cost_b --gap_extension $gap_extension_b ";
	$cmd_identify .= "--hsp_cover $hsp_cover --diff_ratio $diff_ratio --diff_contig_cover $diff_contig_cover --diff_contig_length $diff_contig_length ";
	Util::process_cmd($cmd_identify);

	# delete temp files and log files 
	system("rm -r *.log");
	system("rm -r temp");
}


#################################################################
# subroutine							#
#################################################################

=head2
 files_combine2 : combine known and unknown virus
=cut
sub files_combine2
{
	my $file_list = shift;
	my $fh = IO::File->new($file_list) || die "Can not open innput file list $file_list $!\n";
	while(<$fh>)
	{
		chomp;
		my $file = $_;
		my $file_aligned	= $file.".aligned";
		my $file_assemblied	= $file.".assemblied";
		my $file_combined	= $file.".combined";
		Util::process_cmd("cat $file_aligned $file_assemblied > $file_combined");
	}
	$fh->close;
}
