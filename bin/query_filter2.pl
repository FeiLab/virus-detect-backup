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
  print "usage: query_filter2.pl inputfile1 inputfile2 output1\n";
  exit(0);
}

our $input1 = $ARGV[0]; #这个是作为index的文件
our $input2 = $ARGV[1]; #这个是blast table
our $output = $ARGV[2]; #这是要得到的输出文件，包括了所有提取的序列

open(IN1, "$input1");
my %hit_index; #保存在第一个文件中出现的hit，作为reference

while (<IN1>) {
	chomp; 
	my @cols = split(/\t/, $_); 
	defined $hit_index{$cols[0]} or $hit_index{$cols[0]} = 1; #把所有出现的hit用hash保存
}
close(IN1);

open(IN2, "$input2");
open(OUT, ">$output");
my $last_query="";
my $current_query;

my $high_evalue;#保存每个query最高的identity，其实就是第一次出现的(因为已排序）
my $current_evalue;

my $high_identity;#保存每个query最高的identity
my $current_identity;

my $high_record;#保存前一个query的最好hsp整行记录
my $if_has_ref=0;#表示前一个query或当前query没有reference

while (<IN2>) {
	my @cols = split(/\t/, $_);
	$current_query=$cols[0];
	$current_identity=$cols[5];
	$current_evalue=$cols[6];	
	if($current_query ne $last_query){#每次出现一个新的query
		if($if_has_ref==0 && $last_query ne ""){#先看前一个query的所有hsp有无ref，如果都没有ref
			print OUT $high_record;#就把前一个query中所有hsp中最好hsp输出
		}
		#处理完上个query，开始分析当前query的第一个hsp
		$if_has_ref=0;   #表示当前query还没有找到ref
		if (defined($hit_index{$cols[2]})){#如果第一个hsp就有ref
			$high_identity=$cols[5];#这个就是最高的identity
			$high_evalue=$cols[6];	#这个就是最高的evalue
			print OUT $_;           #输出这条记录
			$if_has_ref=1;	        #表示当前query找到了ref		
		}else{
			$high_record=$_; #如果第一个没有ref，保留为最高的，如果本次的所有hsp都没有ref，就在下次输出这个			
		}
	}else{#对于当前query后续的hsp
		if($if_has_ref==1){#如果当前query已经找到了ref
			if (defined($hit_index{$cols[2]}) && $current_evalue>=$high_evalue && $current_identity>=$high_identity)
			{
				print OUT $_;#只输出有ref同时和最好的hsp同样identity和evalue的结果(注意用>=,后面的有可能更高)			
			}
		}else{#如果当前query还没找到ref
			if (defined($hit_index{$cols[2]})){#一旦找到ref
				$high_identity=$cols[5];#这个就是最高的
				$high_evalue=$cols[6];		
				print OUT $_;
				$if_has_ref=1;	#表示当前的query找到了ref
			}
		}
	
	}
	$last_query=$cols[0];
}
close(IN2);
close(OUT);



