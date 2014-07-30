#!/usr/bin/perl

use strict;   
use warnings; 
use IO::File;
use Getopt::Long; 

my $usage = <<_EOUSAGE_;

fastq_clipper.pl -- clip the adapter using levenshtein or hamming algorithm

Author: Shan Gao

usage:
  fastq_clipper.pl --filelist <FILE> --suffix <String> --adapter <String> --adapterLen INT --distance <INT> --levenshtein

  Required:
  -i / --filelist       a file containing all fastq file name for processing without suffix
  -a / --adapter        the reverse complimentary of 3' end adapter(5'->3')
  -l / --adapterLen     the adapter offset from the 5' end
  -d / --distance       the largest edit distance can be allowed between a adapter and reads
  -e / --levenshtein    using levenshtein distance, default is hamming distance

  Options(1):
  -s / --suffix	    	the suffix of the input files

_EOUSAGE_
	;

our ($filelist, $suffix, $adapter, $adapter_offset, $distance, $levenshtein);

# description of vars
# $filelist;		# list of file name, just prefix
# $suffix="fastq";	# suffix of file name
# $adapter;		# adapter sequence
# $adapter_offset;  	# the top N number of adapter sequence, if not assign this number, full length of adapter will be use
# $distance;		# the distance between reads and adapter, default is 1
# $levenshtein=0;	# using levenshtein to cumpute distance

# get parameter  from command
&GetOptions( 
	'i|filelist=s'	=> \$filelist,
	'a|adapter=s'	=> \$adapter,
	's|suffix=s'	=> \$suffix,
	'l|adapterLen=i'=> \$adapter_offset,
	'd|distance=i'	=> \$distance,
	'e|levenshtein'	=> \$levenshtein 
);

die $usage unless $filelist;
die $usage unless $adapter;
$distance ||= 1;

print "distance is $distance\n";
print "levenshtein is $levenshtein\n";
our $adapter_substr = substr($adapter, 0, $adapter_offset);

check_file($filelist, $suffix);
die;

#################################################################
# main								#
#################################################################
my $fh = IO::File->new($filelist) || die "Error, can not open list file $filelist $!\n";
while(<$fh>)
{
	chomp;
	my $filebase = $_;
	trim_fastq($filebase, $suffix, $distance, $adapter, $adapter_offset);	
}
$fh->close;

print "all the files have been processed\n";

#################################################################
# kentnf: subroutine		                                #
#################################################################
=head1 check_file

 check_file: check if file exist

=cut
sub check_file
{
	my ($list, $suffix) = @_;
	my $fh = IO::File->new($list) || die "Can not open list file $list $!\n";
	while(<$fh>)
	{
		chomp;
		my $filebase = $_;
		my $filename = $filebase.".".$suffix;
		unless (-s $filename) { die "Error! file $filename do not exist!\n"; }
	}
	$fh->close;
}

=head1 trim_fastq

 trim_fastq: 

=cut

sub trim_fastq
{
	my ($basename, $suffix, $distance, $adapter, $adapter_offset) = @_;
	my $input_file	 =$basename.".$suffix";			# input read file                  
	my $output_file1 =$basename.".trimmed".$distance;	# output 1
	my $output_file2 =$basename.".unmatched".$distance;	# output 2
	my $output_file3 =$basename.".null".$distance;		# output 3
	
	my $match_position;	# second element of @match_result, the position of barcode matched to reads 
	my ($id1, $DNAseq, $id2, $qul, $format);

	my ($trimmed_reads, $untrimmed_reads, $null_reads) = (0,0,0);
	
	my $out1 = IO::File->new(">".$output_file1) || die "can't create $output_file1 $!\n";
	my $out2 = IO::File->new(">".$output_file2) || die "can't create $output_file2 $!\n";
	my $out3 = IO::File->new(">".$output_file3) || die "can't create $output_file3 $!\n";

	my $in2 = IO::File->new($input_file) || die "Can not open input sequence file $input_file $!\n";
	while(<$in2>)
	{
		chomp;		
		$id1 = $_;
		if	($id1 =~ m/^>/) { $format = "fasta"; }#$id =~ s/^>//;
		elsif	($id1 =~ m/^@/) { $format = "fastq"; }#$id =~ s/^@//;
		else	{ die "Error in read format for $id1\n"; }

		$DNAseq = <$in2>; chomp($DNAseq);
		$match_position = adapter_match($DNAseq, $adapter_substr);

		if ($format eq "fastq") 
		{ 
			$id2 = <$in2>; chomp($id2);
			$qul = <$in2>; chomp($qul);
		
			if ( $match_position == 0 ) 
			{
				# this reads is fully contaminated, output to null file
				print $out3 $id1."\n".$DNAseq."\n".$id2."\n".$qul."\n";
				$null_reads++;
			}	
			elsif($match_position == length($DNAseq))
			{
				# can not find matched reads, output to unmatched file
				print $out2 $id1."\n".$DNAseq."\n".$id2."\n".$qul."\n";
				$untrimmed_reads++;
			}	
			else						
			{
				# trim part of reads, output to trimmed file
				my $trimmed_DNAseq = substr($DNAseq, 0, $match_position); 
			        my $trimmed_qual = substr($qul, 0, $match_position);        
				print $out1 $id1."\n".$trimmed_DNAseq."\n".$id2."\n".$trimmed_qual."\n";
				$trimmed_reads++;
			}
		}
		else
		{
			if ( $match_position == 0 )
			{
				# this reads is fully contaminated, output to null file
				print $out3 $id1."\n".$DNAseq."\n";
				$null_reads++;
			}
			elsif($match_position == length($DNAseq))
			{
				# can not find matched reads, output to unmatched file
				print $out2 $id1."\n".$DNAseq."\n";
				$untrimmed_reads++;
			}
			else
			{
				# trim part of reads, output to trimmed file
				my $trimmed_DNAseq = substr($DNAseq, 0, $match_position);
				print $out1 $id1."\n".$trimmed_DNAseq."\n";
				$trimmed_reads++;
			}
		}		
	}
	$in2->close;
	$out1->close;
	$out2->close;
	$out3->close;

	print $basename."\t".$null_reads."\t".$trimmed_reads."\t".$untrimmed_reads."\n";
}

=head adapter_match

=cut

sub adapter_match
{
	my ($read_seq, $adapter_substr) = @_;    

	my $read_len = length($read_seq);	
	my $position = $read_len;
	
	my $min_offset=10;
 
	if($levenshtein==1)
	{
		for(my $i = 0; $i<$read_len-$adapter_offset+1; $i++)
		{ 
			my $read_substr = substr($read_seq, $i, $adapter_offset);
			my $editDistance = levenshtein($read_substr,$adapter_substr);
			if ($editDistance <= $distance){
				$position=$i;  #记下这个位置，返回为了下一步trim
				last;          #找到了，就不继续找了，跳出循环	
			}				
		}			
	}
	else
	{
		for(my $i = 0; $i<$read_len-$min_offset+1; $i++)
		{            
			my $read_substr = substr($read_seq, $i, $adapter_offset);
			#注意如果剩下的长度不够$adapter_offset，有多少截取多少而不报错
			my $editDistance = hamming($read_substr,$adapter_substr);
			#$read_substr有可能短于$adapter_offset，短的必须在前面
			if ($editDistance <= $distance){				
				$position=$i;  
				#记下这个位置，返回为了下一步trim
				last;          
				#找到了，就不继续找了，跳出循环	
			}			
		}
	}		
	return $position;
}

=head1 levenshtein

 levenshtein: compute the levenshtein between two reads

=cut
sub levenshtein
{
    my ($s1, $s2) = @_;
    my ($len1, $len2) = (length $s1, length $s2);
    return $len2 if ($len1 == 0);
    return $len1 if ($len2 == 0);

    my %mat;
    for (my $i = 0; $i <= $len1; ++$i)
    {
        for (my $j = 0; $j <= $len2; ++$j)
        {
            $mat{$i}{$j} = 0;
            $mat{0}{$j} = $j;
        }
        $mat{$i}{0} = $i;
    }

    my @ar1 = split(//, $s1);
    my @ar2 = split(//, $s2);

    for (my $i = 1; $i <= $len1; ++$i)
    {
        for (my $j = 1; $j <= $len2; ++$j)
        {
            my $cost = ($ar1[$i-1] eq $ar2[$j-1]) ? 0 : 1;
            $mat{$i}{$j} = min($mat{$i-1}{$j} + 1, $mat{$i}{$j-1} + 1, $mat{$i-1}{$j-1} + $cost);
        }
    }

    return $mat{$len1}{$len2};
}

sub min
{
    my @list = @_;
    my $min = $list[0];
    foreach my $i (@list)
    {
        $min = $i if ($i < $min);
    }
    return $min;
}

=head1 hamming

 hamming: compute hamming between two reads

=cut
sub hamming($$) { length( $_[ 0 ] ) - ( ( $_[ 0 ] ^ $_[ 1 ] ) =~ tr[\0][\0] ) }
