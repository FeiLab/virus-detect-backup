#!/usr/bin/perl -w 
#
use strict; 

-t and !@ARGV and die "perl $0 2.sam\n"; 

our $max_distance = 2; #要求选中的所有hit的编辑距离都不能大于$max_distance
my @lls; 
my $pre_query_name; 
my $nmv; #每个query都有一个最小编辑距离，符合这个距离的hits才保留
$pre_query_name = ''; 
$nmv = -1; 
my ($query_col, $opt_col) = (0, 11); #sam文件一行中的第1列和第12列，分别对应query和option的开始列
while (<>) {
	chomp; s/\s+$//; 
	my @ta = split(/\t/, $_);#ta数组包括sam文件一行中的所有列
	my $query_name = $ta[$query_col];#读入当前query的名称 

	if ($pre_query_name eq '') {#如果前一个query为空，则赋值
		$pre_query_name = $query_name; 
	}elsif ($query_name ne $pre_query_name) {#如果当前query与前一个不同
		$nmv < $max_distance and print STDOUT join("\n", @lls)."\n";#nmv必须小于$max_distance，才保留，注意bwa有时会输出超过预定$max_distance的hits
		$nmv=-1; 
		$pre_query_name = $query_name; 
	}

	for my $each_opt (@ta[$opt_col .. $#ta]) {#扫描11列到最后一列的内容
		if ($each_opt =~ /^NM:i:(\d+)$/) {#从当前行中找到'NM:i:X'模式
			my $distance = $1; #取出query与reference序列之间的distance
			if ($nmv == -1 or $distance < $nmv) {#如果找到更小的距离
				$nmv = $distance; #就更新最小距离
				@lls = ($_); 
			}elsif ($distance == $nmv) {
				push(@lls, $_); 
			}else{
				; 
			}
			last; #找到并处理过NM:i:<N>模式，就跳出来，要求每行的NM:i:<N>只出现一次
		}
	}
}

if ($nmv > -1) {
	$nmv < $max_distance and print STDOUT join("\n", @lls)."\n"; # for MAX NM:i:? 
	$nmv = -1; 
	$pre_query_name = ''; 
}


