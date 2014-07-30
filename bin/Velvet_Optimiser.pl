#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case bundling);
use Bio::SeqIO;
use FindBin;
use Cwd;

my $usage = <<_EOUSAGE_;

#########################################################################################
# Velvet_Optimiser.pl --file_list <File> --input_suffix <String> --file_type <String> --velvet_dir [String] 
                      --objective_type [String] --hash_end <INT> --coverage_end <INT> 
# Required[4]:
#  --file_list             The name of a txt file containing a list of input file names without any suffix
#  --input_suffix          The suffix name of all the input files
#  --hash_end              The max k-mer length for the velvet to try to optimize the objective_type
#  --coverage_end          The max coverage cutoff for the velvet to try to optimize the objective_type
#  
# Options[2]:
#  --file_type             The format of input files(fastq or fasta)  [fastq]
#  --velvet_dir            The directory to install the program velvet [workingdir/bin]
#  --objective_type        The objective type of optimization(maxLen, n50 or avgLen)  [maxLen]
###########################################################################################

_EOUSAGE_
;
	
my $file_list;              
my $input_suffix = '';		# clean or unmapped, do not have "."
my $file_type = "fastq";
my $velvet_dir;			# velvet bin dir
my $objective_type='maxLen';	# objective type (n50, maxLen, avgLen) for optimization
my $hash_end;			# the max range for k-mer
my $coverage_end;		# the max range for coverage cutoff

my $hash_start = 9;
my $coverage_start = 5;

#################################
# set folder and file path	#
#################################
my $WORKING_DIR = cwd();		# current folder
$velvet_dir = ${FindBin::RealBin};	# bin folder
my $tf = $WORKING_DIR."/temp";		# temp folder

###################3#############
# get input parameters		#
#################################
GetOptions(
	'file_list=s'		=> \$file_list,		 
	'input_suffix=s'	=> \$input_suffix,
	'file_type=s'		=> \$file_type,
	'velvet_dir=s'		=> \$velvet_dir,
	'objective_type=s'	=> \$objective_type,			 			 
	'hash_end=i'		=> \$hash_end,
	'coverage_end=i'	=> \$coverage_end
);
			 
die $usage unless ($file_list && $hash_end && $coverage_end);	# required parameters

# main
main: {
    my $sample;
    my $sampleNum=0;
    my $current_folder;
    my $statfile;
    my $objective;						# the objective for each run
    my $max_objective;						# the maximum objective
    my $opt_hash_length=$hash_start; 				# the optimization hash_length when objective is max
    my $opt_coverage=$coverage_start;				# the optimization coverage when objective is max
    my $opt_avgLen=0;                				# the avg Length when objective is max
    open(IN, "$file_list");
    open(OUT1, ">$tf/optimization.log") or die "$!\n"; 		# save optimization information
    open(OUT2, ">$tf/optimization.result") || die "$!\n";		# save final optimization result
    while (<IN>) {
		$sampleNum = $sampleNum+1;
		chomp;
		$sample = $_;		
		print "#processing sample $sampleNum: $sample\n";		
		
		# reset parameters
		$max_objective = 0;
		$opt_hash_length = $hash_start;
		$opt_coverage = $coverage_start;
		
		# optimize k-mer length using fixed coverage
		for(my $i=$hash_start; $i<= $hash_end; $i=$i+2) {
			runVelvet($sample, $i, $coverage_start);
			$current_folder = $sample."_".$i."_".$coverage_start;
			$statfile = $current_folder."/contigs.fa";
			my $aa = contigStats($statfile);        # return hash reference
			if ( defined $aa->{$objective_type} ) {
				$objective = $aa->{$objective_type};
				print OUT1 $i."\t".$coverage_start."\t".$objective."\t".$aa->{avgLen}."\t".$aa->{numSeqs}."\n";
			
				# output best hash length, coverage, objective, avgLength, contigs num
				if ( $objective > $max_objective ) {
					print OUT1 "yes"."\n"; #print yes if the optimization value improved(higher than before)
					$max_objective = $objective;
					$opt_hash_length = $i;
					$opt_avgLen=$aa->{avgLen};
				}
			}
			process_cmd("rm $current_folder -r");		
		} 
		# optimize coverage using fixed k-mer length
		for(my $j=$coverage_start+2; $j<=$coverage_end; $j=$j+1) {  # start from 7, 注意这里从7开始，因为5上面都算过了
			runVelvet($sample,$opt_hash_length,$j);
			$current_folder=$sample."_".$opt_hash_length."_".$j;
			$statfile=$current_folder. "/contigs.fa";
			my $aa=contigStats($statfile);
			if ( defined $aa->{$objective_type} ) {
				$objective=$aa->{$objective_type};
				print OUT1 $opt_hash_length."\t".$j."\t".$objective."\t".$aa->{avgLen}."\t".$aa->{numSeqs}."\n";

				# output the best hast length, coverage, objective, avglength, contigs num
				if($objective>$max_objective){
					print OUT1 "yes"."\n";# print yes of the optimization value improved
					$max_objective=$objective;
					$opt_coverage=$j;
					$opt_avgLen=$aa->{avgLen};
				}
			}
			process_cmd("rm $current_folder -r");
		}       
		print OUT2 $sample."\t".$opt_hash_length."\t".$opt_coverage."\t".$max_objective."\t".$opt_avgLen."\n";
	}
	close(IN);
	close(OUT1);
	close(OUT2);
	print "###############################\n";
    	print "All the samples have been processed by $0\n";
}

# subroutine
sub runVelvet {
	my ($sample1, $hash_length, $cov_cutoff) = @_;
	my $outputDir=$sample1."_".$hash_length."_".$cov_cutoff;

	my $file;
	if ($input_suffix)	{ $file = "$sample1.$input_suffix"; } 
	else 			{ $file = $sample1; }

	process_cmd($velvet_dir."/velveth $outputDir $hash_length -$file_type $file >> $tf/velvet.log");
	process_cmd($velvet_dir."/velvetg $outputDir -cov_cutoff $cov_cutoff -min_contig_lgth 30 >> $tf/velvet.log");	
}
sub process_cmd {
	my ($cmd) = @_;	
	print "CMD: $cmd\n";
	my $ret = system($cmd);	
	if ($ret) {
		print "Error, cmd: $cmd died with ret $ret";
	}
	return($ret);
}
sub contigStats {
	
	my $file = shift;
	#my $minsize = shift;
	
	#print "In contigStats with $file, $minsize\n" if $interested;
	
	my $numseq=0;
	my $avglen=0;
	my $minlen=1E9;
	my $maxlen=0;
	my @len;
	my $toosmall=0;
	my $nn=0;
	
	my $in = Bio::SeqIO->new(-file => $file, -format => 'Fasta');
	while(my $seq = $in->next_seq()){
		my $L = $seq->length;
		#if($L < $minsize){
			#$toosmall ++;
			#next;
		#}
		#count Ns
		my $s = $seq->seq;
		my $n = $s =~ s/N/N/gi;
		$n ||= 0;
		$nn += $n;
		#count seqs and other stats
		$numseq ++;
		$avglen += $L;
		$maxlen = $L if $L > $maxlen;
		$minlen = $L if $L < $minlen;
		push @len, $L;
	}
	@len = sort { $a <=> $b } @len;
	my $cum = 0;
	my $n50 = 0;
	for my $i (0 .. $#len){
		$cum += $len[$i];
		if($cum >= $avglen/2) {
			$n50 = $len[$i];
			last;
		}
	}
	
	my %out;
	if($numseq > 0)
	{
		$out{numSeqs} = $numseq;
		$out{numBases} = $avglen;
		$out{numOK} = ($avglen - $nn);
		$out{numNs} = $nn;
		$out{minLen} = $minlen;
		$out{avgLen} = $avglen/$numseq;
		$out{maxLen} = $maxlen;
		$out{n50} = $n50;
		#$out{minsize} = $minsize;
		$out{numTooSmall} = $toosmall;
	}
	else 
	{
		$out{$numseq} = 0;
	}
	
	#print "Leaving contigstats!\n" if $interested;
	return (\%out);
}
