#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
my $usage = <<_EOUSAGE_;

#########################################################################################
# renameFasta.pl --inputfile <FILE> --outputfile <FILE> --prefix [String]
# Required:
#  --inputfile        a txt file containing a list of input sample names
#  --outputfile       a txt file containing a list of input sample names
#
###########################################################################################

_EOUSAGE_

	;
my $inputfile; #输入文件，必须是fasta格式
my $outputfile;#输出文件，必须是fasta格式
my $prefix;    #原文件中序列名称改为$prefix+数字（自0开始）

&GetOptions( 'inputfile=s' => \$inputfile,
             'outputfile=s' => \$outputfile,
			 'prefix=s' => \$prefix  
			 );

unless ($inputfile&&$outputfile) {
	die $usage;
}


#读fasta file文件，
open(IN, $inputfile);
open(OUT, ">$outputfile");
my $sequenceNum=0;
while(<IN>) {
	if($_ =~ m/^>/) {
		$sequenceNum++;
		print OUT ">".$prefix.$sequenceNum."\n";

	}
	else {
		print OUT $_;#否则，后面的序列需要继续向OUT2输出
	}
}
close(IN);
close(OUT);