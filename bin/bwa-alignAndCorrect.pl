#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Cwd;
use FindBin;
use lib "$FindBin::RealBin/PerlLib";
use Util;
use File::Basename;
use IO::File;
use Bio::SeqIO;

my $usage = <<_EOUSAGE_;

#########################################################################################
# bwa-alignAndCorrect.pl --file_list <FILE> --reference <FILE> --coverage <Float> 
#                 		 --max_dist[INT] --max_open [INT] --max_extension [INT] --len_seed [INT] --dist_seed [INT] --thread_num [INT]
# 
# Required(3):
#  --file_list A txt file containing a list of input file names without any suffix
#  --reference A fasta file containing all the reference sequences
#  --coverage  A reference sequence covered by reads will be output as results
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
#  global vars  #
#################
our $file_list;		# list of sample name
our $reference;		# reference sequence, fasta format, full path
our $coverage;  	# cutoff of mapping reads / reference
our $max_dist = 1; 	# BWA allowed max distance 
our $max_open = 1; 	# BWA allowed max gap opens
our $max_extension = 1; # BWA allowed max gap extension
our $len_seed = 15; 	# length of seed
our $dist_seed = 1; 	# BWA allowed seed max distance
our $thread_num = 8; 	# thread number
our $output_suffix;	# output suffix

################################
# set folder and file path     #
################################
our $WORKING_DIR	= cwd();				# current folder
our $TEMP_DIR		= $WORKING_DIR."/temp";			# temp foldr
our $DATABASE_DIR	= ${FindBin::RealBin}."/../databases";	# database folder
our $BIN_DIR		= ${FindBin::RealBin};			# programs folder
our $seq_info		= $DATABASE_DIR."/vrl_genbank.info";	# genbank info

my $tf = $TEMP_DIR; # short name of temp folder for easy to do

###################
# get input paras #
###################
GetOptions( 
	'file_list=s'		=> \$file_list,
	'reference=s' 		=> \$reference,
	'coverage=f' 		=> \$coverage,
	'max_dist=i' 		=> \$max_dist,
	'max_open=i' 		=> \$max_open,
	'max_extension=i' 	=> \$max_extension,
	'len_seed=i' 		=> \$len_seed,
	'dist_seed=i' 		=> \$dist_seed,			 
	'thread_num=i' 		=> \$thread_num,
	'output_suffix=s'	=> \$output_suffix
);

die $usage unless ($file_list && $reference && $coverage && $output_suffix);	# required parameters

#################
# main          #
#################
main: {

    # create bwa and fasta index file for reference file (plant virus)   
    Util::process_cmd("$BIN_DIR/bwa index -p $reference -a bwtsw $reference 2> $tf/bwa.log") unless (-e "$reference.amb"); 
    Util::process_cmd("$BIN_DIR/samtools faidx $reference") unless (-e "$reference.fai");
    
    # parse samples in file list 
    my ($i, $sample, $output_file);
    $i=0;
    open(IN, "$file_list") || die $!;
    while (<IN>) {
		$i=$i+1;
		chomp;
		$sample = $_;
		$output_file = $sample.".".$output_suffix;
		die "Error, file $sample does not exist\n" unless -s $sample;		
		#print "# processing the $i sample -- $sample using $0\n";

		#aligment -> sam -> bam -> sorted bam -> pileup
		Util::process_cmd("$BIN_DIR/bwa aln -n $max_dist -o $max_open -e $max_extension -i 0 -l $len_seed -k $dist_seed -t $thread_num $reference $sample 1> $sample.sai 2>> $tf/bwa.log") unless (-s "$sample.sai");
		Util::process_cmd("$BIN_DIR/bwa samse -n 10000 -s $reference $sample.sai $sample 1> $sample.sam 2>> $tf/bwa.log") unless (-s "$sample.sam");

		# filter out unmapped reads
		# filter out 2nd hits, and only keep the best hits of reads alignment to reference		
		filter_SAM_unmapped($sample.".sam");
		filter_SAM_2nd_hits($sample.".sam");

		# sort sam to bam and generate pileup file
		Util::process_cmd("$BIN_DIR/samtools view -bt $reference.fai $sample.sam > $sample.bam 2> $tf/samtools.log") unless (-s "$sample.bam");
		Util::process_cmd("$BIN_DIR/samtools sort $sample.bam $sample.sorted 2> $tf/samtools.log") unless (-s "$sample.sorted.bam");
		Util::process_cmd("$BIN_DIR/samtools mpileup -f $reference $sample.sorted.bam > $sample.pre.pileup 2> $tf/samtools.log") unless (-s "$sample.pre.pileup");

		# if coverage of any reference is lower than ratio cutoff (default : 0.3 ), remove the alignment for this reference
		my $depth_cutoff;
		$depth_cutoff = pileup_filter("$sample.pre.pileup", "$seq_info", "$coverage", "$sample.pileup") unless (-s "$sample.pileup");
		# *** notice, could rewrite this function *** 

		# get the size of pileup file	
		my $file_size = -s "$sample.pileup";
		if($file_size != 0)
		{
			# get continous fragment for depth>=1 and length>=40, output file name has '1'
			Util::process_cmd("java -cp $BIN_DIR extractConsensus $sample 0 40 1");			
			my $result_contigs = $sample.".contigs1.fa";			# this is a corrected sequences
			my $count = `grep \'>\' $result_contigs | wc -l`; chomp($count);# get seq number
			Util::print_user_submessage("$count of contigs were retrieved from aligned reads");

			if($count!=0){ 							# if seq number > 0, move seq file to align folder
				#system("$BIN_DIR/renameFasta.pl --inputfile $result_contigs --outputfile $tsample.contigs2.fa --prefix ALIGNED");# the seq name should be formatted
				renameFasta($result_contigs, "$sample.contigs2.fa", "ALIGNED");
				Util::process_cmd("mv $sample.contigs2.fa $output_file");# move the seq file to folder
				Util::process_cmd("rm $result_contigs");
			}
			else
			{
				Util::process_cmd("mv $result_contigs $output_file");	# if seq number = 0, move file to align folder
			}
		}	
		else{	
			Util::print_user_submessage("None of contigs were retrieved from aligned reads");# file size is 0 or sample.pileup is not exist
			Util::process_cmd("touch $output_file");
		}
		
		system("rm $sample.sai");
		system("rm $sample.sam");
		system("rm $sample.bam");
		system("rm $sample.sorted.bam");
		system("rm $sample.pre.pileup");
		system("rm $sample.pileup");	
		system("rm $tf/bwa.log");		
	}
	close(IN);
}

#################
# subroutine    #
#################

sub pileup_filter 
{
	my ($inputfile, $seq_info, $ratio, $outputfile) = @_;

	# get seq length 
	open(IN1, "$seq_info" ) || die $!;
	my %hash;
	my @a;
	my $id;
	while(my $line=<IN1>){
		chomp($line);
		@a= split(/\t/,$line);
		$id=$a[0];			# seq ID
		$hash{$id}=$a[1];	# length
	}
	close IN1;

	open(IN2, "$inputfile" ) || die $!;
	open(OUT, ">$outputfile" ) || die $!;

	my $pre_reference="a";
	my $current_reference;
	my $i=1;				#从1开始计数，得到当前reference sequence的长度
	my @data;
	my $len_threshold=0;
	my $line;
	my $totalDepth=0;
	my $totalPositions=0;
	while (<IN2>) {
	    $line = $_;    
	    chomp($line);
	    my @a = split( /\t/, $line);
	    $current_reference= $a[0];			#取当前的参考序列id
	    $totalDepth=$totalDepth+$a[3];		#累加这个文件中所有位点的depth
	    $totalPositions=$totalPositions+1;	#累加这个文件中所有位点的数量
	    if($current_reference eq $pre_reference){#如果当前的参考序列id与前面的相同，表示处于同一个参考序列内部
			$i++;#@a数组的指针加1，后面会加入对应的元素	
	    }
	    else{#id改变，开始判断并处理
		if($pre_reference ne "a"){#如果不是第1个query，第1个实际是空
			$len_threshold=$hash{$pre_reference}*$ratio;
			#print $pre_reference."\t".$len_threshold."\t".$i."\n";#输出上一个reference上覆盖碱基的实际长度，和长度阈值（总长*要求的覆盖度）
			if($i>=$len_threshold){#如果此序列的实际长度超过阈值
				foreach my $data (@data) {
					print OUT "$data\n";#输出到结果文件中
				}
			}
			$i=1;#重新开始计数
			@data=();  
		} 
	    }
	    push @data, $line;#每次循环都要存pileup文件的一行
	    $pre_reference=$current_reference;
	}
	
	#print "the average depth of $inputfile is \t", 1.0*$totalDepth/$totalPositions, "\n";

	#别忘了处理剩下的部分内容
	$len_threshold=$hash{$current_reference}*$ratio;
	#print $current_reference."\t".$len_threshold."\t".$i."\n";
	if($i>=$len_threshold){
		foreach my $data (@data) {
		    print OUT "$data\n";
		}
	}
	close(IN2);
	close(OUT);
	my $aveDepth=1.0*$totalDepth/$totalPositions;
	return	$aveDepth;
}

=head2
 renameFasta: rename the fasta file with spefic prefix
=cut
sub renameFasta
{
	my ($input_fasta_file, $output_fasta_file, $prefix) = @_;

	my $seq_num = 0;
	my $out = IO::File->new(">".$output_fasta_file) || die $!;
	my $in = Bio::SeqIO->new(-format=>'fasta', -file=>$input_fasta_file);
       	while(my $inseq = $in->next_seq)
	{
		$seq_num++;
		print $out ">".$prefix.$seq_num."\n".$inseq->seq."\n";
	}
	$out->close;
}


=head2
 
 filter_SAM_unmapped: filter out unmapped hit base on cigar (cigar == 4)

=cut
sub filter_SAM_unmapped
{
	my $input_SAM = shift;
	
	my $temp_SAM = $input_SAM.".temp";
	my ($total_count, $filtered_count) = (0, 0);

	my $in  = IO::File->new($input_SAM) || die $!;
	my $out = IO::File->new(">".$temp_SAM) || die $!;
	while(<$in>)
	{
		chomp;
		if ($_ =~ m/^@/) { print $out $_."\n"; next; }
		my @a = split(/\t/, $_);
		if ( $a[1] == 4 ) { $filtered_count++; }
		else {	print $out $_."\n"; }
		$total_count++;	
	}
	$in->close;
	$out->close;
	Util::process_cmd("mv $temp_SAM $input_SAM");
	#print STDERR "This program filtered $filtered_count out of $total_count reads (" . sprintf("%.2f", $filtered_count / $total_count * 100) . ") as unmapped reads, only for BWA\n";
}

=head2

 filter_SAM_2nd_hits: filter out 2nd hits, only kept 1st hits alignment for each query id

=cut
sub filter_SAM_2nd_hits
{
	my $input_SAM = shift;
	my $temp_SAM = $input_SAM.".temp";

	my ($query_col, $opt_col) = (0, 11); 	# query and option column number for sam
	my $max_distance = 2;			# set $max_distance for all selected hits
	my $bestEditDist = -1;			# set best edit distance
	my @alignment = ();			# alignment to array
	my $pre_query_name = '';		# previous query name
	my ($total_count, $kept_align) = (0,0);

	my $in  = IO::File->new($input_SAM) || die $!;
	my $out = IO::File->new(">".$temp_SAM) || die $!;
	while(<$in>)
	{
		chomp;
		if ($_ =~ m/^@/) { next; #print $out $_."\n"; next;
	       	}
		my @a = split(/\t/, $_);

		my $query_name = $a[$query_col];

		if ($query_name ne $pre_query_name) 
		{
			# parse the pre results
			foreach my $align (@alignment)
			{
				my $editDistance;
				if ($align =~ m/\tNM:i:(\d+)/) { $editDistance = $1; }	
				else { die "Error, this alignment info do not have edit distance : $align\n"; }
				if ($editDistance == $bestEditDist) { print $out $align."\n"; $kept_align++; }
			}	

			# init vars;
			@alignment = ();
			$bestEditDist = -1;
			$pre_query_name = $query_name;
		}

		my $distance;
		if ($_ =~ /\tNM:i:(\d+)/) { $distance = $1; }
		else { die "Error, this alignment info do not have edit distance : $_\n"; }
		next if $distance >= $max_distance;
		if ($bestEditDist == -1) { $bestEditDist = $distance; }
		if ($distance < $bestEditDist) { $bestEditDist = $distance; }	
		push (@alignment, $_);

		$total_count++;	
	}
	$in->close;

	# parse final query recoed
	if (scalar(@alignment) > 0) 
	{
 		foreach my $align (@alignment)
		{
			my $editDistance;
			if ($align =~ m/\tNM:i:(\d+)/) { $editDistance = $1; }	
			else { die "Error, this alignment info do not have edit distance : $align\n"; }
			if ($editDistance == $bestEditDist) { print $out $align."\n"; $kept_align++; }
		}	
	}

	$out->close;

	Util::process_cmd("mv $temp_SAM $input_SAM");
	my $filtered_count = $total_count - $kept_align;
	#print STDERR "This program filtered $filtered_count out of $total_count reads (" . sprintf("%.2f", $filtered_count / $total_count * 100) . ") as 2ndhits reads, only for BWA\n";
}


