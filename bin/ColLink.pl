#!/usr/bin/perl

#最终输出的文件中id与inputfile2中保持一致
#文件inputfile1中id与inputfile2中id一样的，对应的信息保留
#每个文件内部的各列不是必须选择的，而且顺序可以变换
use strict;
use Getopt::Long;
my %opts;
GetOptions(\%opts,
	"keyC1:s","keyC2:s",
	"col1:s","col2:s",
	"f1:s",
	"out:s",
	"add!",
	"sign:s",
	"link!", # used for join two files. 
	"help!");
my $info = <<INFO;
############################################################
./ColLink.pl <inputfile2> -keyC1 0 -Col1 0,1,.. -keyC2 0 -Col2 0,1,.. -f1 <inputfile1> 

-add        if given, -Col2 全部保留
-link       if given, -Col1 全部保留
-sign		\'yes,no\'
############################################################
INFO

$opts{help} and die "$info";
!@ARGV and -t and die "$info"; 

my @Col1 = split(/,/,$opts{col1});
my @Col2 = split(/,/,$opts{col2});
my @KC1 = (0); 
my @KC2 = (0); 
defined $opts{keyC1} and @KC1 = split(/,/, $opts{keyC1});
defined $opts{keyC2} and @KC2 = split(/,/, $opts{keyC2}); 
#my $kc1 = $opts{keyC1};
#my $kc2 = $opts{keyC2};
my $file1 = $opts{f1};
my $out = $opts{out};
my ($yS,$nS);
if (defined $opts{sign}) {
	($yS,$nS) = split(/,/,$opts{sign});
	$yS =~ s/^\s+//;
}

my $f1_col_fit = 0; 
my @f1_link_col; 

open (F1,'<',$file1) or die "f1?\n$info";
my %need;		# 
while (<F1>) {
	chomp;
	my @temp = split(/\t/,$_);
	my $tk = join("\t", @temp[@KC1]); 
	if ($opts{link}) {
#		if ($f1_col_fit == 0) {
#			F1_COL: 
#			for my $id1 (0..$#temp) {
#				for my $id2 (@KC1) {
#					$id1 == $id2 and next F1_COL; 
#				}
#				push(@f1_link_col, $id1); 
#			}# end of F1_COL; 
#			$f1_col_fit = 1; 
#		}
#		$need{$tk} = join("\t", @temp[@f1_link_col]); 
		my @tt; 
		F1_COL:
		for my $id1 (0..$#temp) {
			for my $id2 (@KC1) {
				$id1 == $id2 and next F1_COL;
			}
			push(@tt, $id1);
		}# end of F1_COL;
		$need{$tk} = join("\t", @temp[@tt]);
	}else{
		$need{$tk} = join("\t",@temp[@Col1]);
	}
}
close F1;
#open(OUT,'>',$out) or die $info;
while (<>) {
	chomp;
	my @temp = split(/\t/,$_);
	my ($line,$add);
	my $tk = join("\t", @temp[@KC2]); 
	if (defined $need{$tk}) {
		$add = $need{$tk};
	}else{
		my @tt;
		$#tt = $#Col1;
		$add = join("\t",@tt);
	}
	if ($opts{add}) {
		$line = $_;
	}else{
		$line = join("\t",@temp[@Col2]);
	}
	if (defined $opts{sign}) {
		$add = (defined $need{$tk})?$yS:$nS;
	}
	print STDOUT "$line\t$add\n";
}
#close OUT;



