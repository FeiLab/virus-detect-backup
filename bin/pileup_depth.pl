#!/usr/bin/perl -w 
use strict; 
# usage: pileup_depth.pl inputfile outputfile
# 从pileup中提取每个chr的平均depth信息


if (@ARGV < 1)
{
  print "usage: pileup_depth.pl inputfile outputfile\n";
  exit(0);
}

our $input = $ARGV[0]; #去掉文件后缀名称
our $output = $ARGV[1];

my $last_chr = ''; 
my %depth;#这个哈希表，映射chr名称与平均depth
my $mean_depth;
my $total_depth = 0; 
my $total_positions = 1; 
open(IN, "$input");
while (<IN>) {
	my @cols = split;
	if ($cols[0] ne $last_chr) { #如果当前chr与前一个不同
                if ($last_chr){      #并且不是第一个chr
			#$depth{$last_chr} = 1.0*$total_depth/$total_positions;#存储chr对应的平均depth
			$mean_depth = 1.0*$total_depth/$total_positions;
			push(@{$depth{$last_chr}}, ($mean_depth,$total_positions)); 
		}			
		$total_depth = $cols[3];#重新初始化
		$total_positions = 1;   #重新初始化
	}
	else{
		$total_depth = $total_depth + $cols[3];#累加第4列深度信息
		$total_positions = $total_positions + 1;#累加位点信息
	}
        $last_chr = $cols[0];	
}
#$depth{$last_chr} = 1.0*$total_depth/$total_positions;
$mean_depth = 1.0*$total_depth/$total_positions;
push(@{$depth{$last_chr}}, ($mean_depth,$total_positions));  
close(IN);

open(OUT, ">$output");
foreach my $chr (keys(%depth))
{
	print OUT $chr."\t".$depth{$chr}[0]."\t".$depth{$chr}[1]."\n";#输出hit名称，平均深度和hit被覆盖的位置总数
}
close(OUT);


