#!/usr/bin/perl -w 
use strict; 
# usage: ./hit_cov1.pl inputfile outputfile identity query_cov > output2
# 来自blast_parse_table2.pl程序的输出
# query_name	query_length	hit_name	hit_length	hsp_length	identity	evalue	score	strand	query_start	query_end	hit_start	hit_end
# 其中2、3、4、5、11、12列(0开始)是必须的
# 将blast结果的table中的每个hit被query覆盖的hsp全部加起来的长度，然后除以该hit总长度，得到一个ratio(结果)
# 还可以得到，所有这个组成这个hit的query
# 对于每个hit，只有满足相似度高于identity，query被覆盖query_cov的hsp才予以考虑
if (@ARGV < 2)
{
  print "usage: hit_cov1.pl input output1 identity query_cov > output2\n";
  exit(0);
}

our $input = $ARGV[0]; #去掉文件后缀名称
our $output = $ARGV[1];
our $identity = $ARGV[2];
our $query_cov = $ARGV[3];

open(IN, "$input");
open(OUT, ">$output");

my %blk; 
my %hit_len; #hit的长度
while (<IN>) {
	chomp; 
	my @ta = split(/\t/, $_);
	my $query_covered= 1.0*$ta[4]/$ta[1]; #hsp length除以query length	
	if ($ta[5]>=$identity && $query_covered>=$query_cov){#只保存符合identity要求的hsp
		push(@{$blk{$ta[2]}}, [@ta[11,12]]); #建立hit和(hit_start,hit_end)之间的映射
		$blk{$ta[2]}[-1][2]{$ta[0]}=1; 
		defined $hit_len{$ta[2]} or $hit_len{$ta[2]} = $ta[3];#建立hit和hit_length之间的映射 
	}
}
close(IN);

#print OUT join("\t", qw/hit_name hit_len total_cov cov_rate% query_names/)."\n";   #向结果文件输出的各列名称
print join("\t", qw/hit_name block_len block_start block_end covered_len query_names/)."\n"; #向标准输出输出的各列名称

for my $tk (sort keys %blk) {#先根据hit排序
	my @o; #存储没有overlap的hit上的block
	for my $ar (sort {$a->[0]<=>$b->[0] || $a->[1]<=>$b->[1] || $a->[2] cmp $b->[2];} @{$blk{$tk}}) {
		if (scalar(@o) == 0) {
			push(@o, [@$ar]); 
		}elsif ($o[-1][1] >= $ar->[0]-1) {#把一个hit上重叠的区域合并
			for my $qname (keys %{$ar->[2]}) {
				$o[-1][2]{$qname}=1; 
			}
			$o[-1][1] < $ar->[1] and $o[-1][1] = $ar->[1]; 
		}else{
			push(@o, [@$ar]); 
		}
	}

	my $total_cov = 0;#hit被覆盖的总长度 
	my %aa; 
	for my $ar (@o) {
		my @query_names = sort keys %{$ar->[2]}; 
		print join("\t", $tk, $hit_len{$tk}, $ar->[0], $ar->[1], $ar->[1]-$ar->[0]+1, join(",", @query_names))."\n"; #向标准输出输出各列
		#hit名称，hit长度，hit一个block的起点，hit一个block的终点，hit一个block的长度
		$total_cov += ($ar->[1]-$ar->[0]+1); #一个hit上所有非重叠block的长度之和
		@aa{@query_names} = (1) x scalar(@query_names); #将所有contigs的名称存入hash表
	}
	print OUT join("\t", $tk, $hit_len{$tk}, $total_cov, $total_cov/$hit_len{$tk}*1.0, join(",", sort keys %aa), scalar(keys %aa))."\n"; #向结果文件输出各列
	#hit名称，hit长度，hit被覆盖bp数，hit被覆盖的百分比，所有contigs名称，contigs数量
}
close(OUT);


