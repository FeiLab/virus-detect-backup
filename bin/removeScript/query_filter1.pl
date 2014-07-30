#!/usr/bin/perl -w 
use strict; 
# 来自blast_parse_table2.pl程序的输出
# query_name	query_length	hit_name	hit_length	hsp_length	identity	evalue	score	strand	query_start	query_end	hit_start	hit_end
# 其中0、1、9、10是必须的
# 将blast结果的table中的query被hit覆盖的全部长度加起来，然后除以query总长度，得到一个ratio
# 所有的hit必须要求大于一定的identity
# 把query被覆盖超过一定ratio的记录下来，然后从一个fasta文件中提取记录中不包括的相应序列，输出到output1
# 剩下的序列的名称输出到ouput2
# 如果hit都是病毒，可以用这个ratio来判断所得到contig是不是病毒
if (@ARGV < 2)
{
  print "usage: query_filter1.pl inputfile1 inputfile2 output1 identity min_ratio > output2\n";
  exit(0);
}

our $input1 = $ARGV[0]; #
our $input2 = $ARGV[1]; #
our $output = $ARGV[2]; #这是要得到的输出文件，包括了所有提取的序列
our $identity = $ARGV[3];#需要累积的hsp的identity最小值
our $min_ratio = $ARGV[4];#一个query累计被covered的最小比例

open(IN1, "$input1");
my %blk; 
my %query_len; 
my %query_filtered;
while (<IN1>) {
	chomp; 
	my @ta = split(/\t/, $_); 
	if ($ta[5]>=$identity){#只保存符合identity要求的记录
		push(@{$blk{$ta[0]}}, [@ta[9,10]]); #建立query和(query_start,query_end)之间的映射
		defined $query_len{$ta[0]} or $query_len{$ta[0]} = $ta[1];#建立query和query_length之间的映射 
	}
}
close(IN1);

#print join("\t", qw/query_name block_len block_start block_end covered_len/)."\n"; #向标准输出输出的各列名称
#print OUT join("\t", qw/query_name query_len total_cov cov_rate%/)."\n";   #向结果文件输出的各列名称

for my $tk (sort keys %blk) {#先根据query排序，然后提取所有需要filtered掉的query名称
	my @o; #存储没有overlap的query上的block
	for my $ar (sort { $a->[0]<=>$b->[0] || $a->[1]<=>$b->[1];} @{$blk{$tk}}) {
		if (scalar(@o) == 0) {
			push(@o, [@$ar]); 
		}elsif ($o[-1][1] >= $ar->[0]-1) {#把一个query上重叠的区域合并
			$o[-1][1] < $ar->[1] and $o[-1][1] = $ar->[1]; 
		}else{
			push(@o, [@$ar]); 
		}
	}

	my $total_cov = 0; 
	for my $ar (@o) {
		#print join("\t", $tk, $query_len{$tk}, $ar->[0], $ar->[1], $ar->[1]-$ar->[0]+1)."\n"; #向标准输出输出的各列
		#query名称，query长度，query一个block的起点，query一个block的终点，query一个block的长度
		$total_cov += ($ar->[1]-$ar->[0]+1); #一个query上所有非重叠block的长度之和
	}
	my $ratio=int($total_cov/$query_len{$tk}*10000+0.5)/100;
	if($ratio >= $min_ratio){#所有符合条件的query,都需要存起来
		defined $query_filtered{$tk} or $query_filtered{$tk} = 1;
	}
	#向结果文件输出的各列
	#query名称，query长度，query被覆盖的总长度，query被覆盖的百分比
}

#从input2中把不包括在%query_filtered中的序列提取出来，输出到OUT
open(IN2, $input2);
open(OUT, ">$output");
my $flag = "off";
while(<IN2>) {
	if($_ =~ m/^>/) {
		my $head = $_;
		chomp($head);
		$head=~s/>//;

		if(defined $query_filtered{$head}) {#如果包括这个name
			print $head."\t";#就输出到标准输出
			$flag = "on";#同时改变标志，表示后面的序列需要继续向OUT1输出
		}
		else {#如果不包括这个name
			print OUT $_;#就输出到OUT
			$flag = "off";#同时改变标志，表示后面的序列需要继续向OUT2输出
		}
	}
	else {
		if($flag eq "on") {#表示为"on"
			print $_;#后面的序列需要继续向标准输出输出
		}
		else {
			print OUT $_;#否则，后面的序列需要继续向OUT输出
		}
	}
}
close(IN2);
close(OUT);



