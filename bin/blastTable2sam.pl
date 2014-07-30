#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;


#把blast得到的paired结果文件转为sam文件类型
if (@ARGV ==0)
{
  print "usage: blast2sam.pl input > output\n";
  exit(0);
}

our $input = $ARGV[0]; 
my $hit_name;
my $hit_len;
my %name_len;
#先通过一次循环，提取所有的hit，保存到hash表里
open(IN, "$input");
while (<IN>) {
	my @cols= split(/\t/,$_);
	$hit_name = $cols[2];
	$hit_len = $cols[3];
	defined $name_len{$hit_name} or $name_len{$hit_name} = $hit_len;	
}
close(IN);
#输出这些hit的信息，作为sam文件的header
foreach my $key (keys(%name_len))
{
	print "\@SQ\tSN:".$key."\tLN:".$name_len{$key}."\n";#输出hit名称，平均深度
}

my ($qlen, $slen, $q, $s, $qbeg, $qend, @sam, @cigar, @cmaux);
@sam = (); @sam[4,6,7,8,10] = (255, '*', 0, 0, '*');#这些列对所有数据都是固定值
open(IN, "$input");
while (<IN>) {
	chomp;
	my @cols= split(/\t/,$_);
	$sam[0] = $cols[0]; #query的名称
	if($cols[8] eq "-1"){
		$sam[1] = 0x10;}
	else{$sam[1] = 0x00;}
	$sam[2] = $cols[2]; #subject的名称
	$qlen = $cols[1]; #query的长度
	$qend = $cols[10]; #query的最后一个碱基
	$sam[3] = $cols[11]; #subject的开始碱基，也就是query aligned的位置
	@cigar = ();
	@cmaux = (0, 0);
	my @q_seqs= split(/,/,$cols[14]);#解析所有query的分段序列
	my @s_seqs= split(/,/,$cols[15]);#解析所有query的分段序列
	my $l = scalar(@q_seqs);#取数组size,经核对最后一个","后面的空字符不在数组内
	$sam[9] = '';#每次必须初始化
	for (my $i = 0; $i < $l; ++$i)
	{
		$q_seqs[$i] =~ /Query\:\s(\d+)\s*(\S+)\s(\d+)/;
		$q = $2;
		my $x = $q;
		$x =~ s/-//g; 
		$sam[9] .= $x;#不断连接query序列片段
		
		$s_seqs[$i] =~ /Sbjct\:\s(\d+)\s*(\S+)\s(\d+)/;
		$s = $2;
		&aln2cm(\@cigar, \$q, \$s, \@cmaux);		
	}	
	#提取附加行信息
	my ($as, $ev) = (int($cols[7] + .499), $cols[6]);
	$ev = "1$ev" if ($ev =~ /^e/);
	@sam[11,12] = ("AS:i:$as", "EV:Z:$ev");
	&blast_print_sam(\@sam, \@cigar, \@cmaux, $qlen - $qend);
}
close(IN);


sub blast_print_sam {
  my ($sam, $cigar, $cmaux, $qrest) = @_;
  push(@$cigar, $cmaux->[1] . substr("MID", $cmaux->[0], 1));
  push(@$cigar, $qrest . 'H') if ($qrest);
  if ($sam->[1] & 0x10) {
	@$cigar = reverse(@$cigar);
	$sam->[9] = reverse($sam->[9]);
	$sam->[9] =~ tr/atgcrymkswATGCRYMKSW/tacgyrkmswTACGYRKMSW/;
  }
  $sam->[9] = '*' if (!$sam->[9]);
  $sam->[5] = join('', @$cigar);
  print join("\t", @$sam), "\n";
}

sub aln2cm {
  my ($cigar, $q, $s, $cmaux) = @_;
  my $l = length($$q);
  for (my $i = 0; $i < $l; ++$i) {
	my $op;
	# set $op
	if (substr($$q, $i, 1) eq '-') { $op = 2; }#如果query中发现一个'-'，应该是"D"
	elsif (substr($$s, $i, 1) eq '-') { $op = 1; }#如果subject中发现一个'-'，应该是"I"
	else { $op = 0; }#否则就是"M"
	# for CIGAR
	if ($cmaux->[0] == $op) {
	  ++$cmaux->[1];#记录"MID"每种类型的数量
	} else {
	  push(@$cigar, $cmaux->[1] . substr("MID", $cmaux->[0], 1));
	  $cmaux->[0] = $op; $cmaux->[1] = 1;
	}
  }
}
