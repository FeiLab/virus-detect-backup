#!/usr/bin/perl -w 
#findRedundancy(),提取一个query对应的所有hit的所有hsp
use strict; 
use IO::File; 
use Getopt::Long; #这个模块可以接受完整参数
use Cwd;
my $usage = <<_EOUSAGE_;

#########################################################################################
# removeRedundancy.pl --input <FILE> --strand_specific --min_overlap [INT] --max_end_clip [INT] 
#                      --cpu_num [INT] --mis_penalty [INT] --gap_cost [INT] --gap_extension [INT]
# Required[1]:
#  --input              A name of an input file containing sequences in fasta format
#
# Megablast-related options[7]:
#  --strand_specific    Only for sequences assembled from strand specific RNA-seq [Not selected]
#  --min_overlap        The minimum overlap length between two contigs to be combined [30]
#  --max_end_clip       The maximum length of end clips [4]
#  --cpu_num         Number of processors to use [8] 
#  --mis_penalty     Penalty for a nucleotide mismatch [-1]
#  --gap_cost        Cost to open a gap [2] 
#  --gap_extension   Cost to extend a gap [1]     
##########################################################################################

_EOUSAGE_
	;
#################
##   全局变量  ##
#################	
our $input;            #需要处理的文件
our $strand_specific;  #专门用于strand specific转录组
our $min_overlap = 30; #hsp合并时，最短的overlap
our $max_end_clip = 4; #hsp合并时，两端允许的最小clip
our $min_identity;     #hsp合并需要满足的最小identity，随hsp长度而调整

our $filter_query = "F";     #默认不需要去除简单序列
our $word_size = int($min_overlap/3);
our $cpu_num = 8;          #megablast使用的cpu数目
our $hits_return = 10;     #megablast返回的hit数目，这个不需要用户设定
our $mis_penalty = -1;     #megablast中，对错配的罚分，必须是负整数
our $gap_cost = 2;         #megablast中，对gap open的罚分，必须是正整数
our $gap_extension = 1;    #megablast中，对gap open的罚分，必须是正整数

################################
##   设置所有目录和文件的路径 ##
################################
our $WORKING_DIR=cwd();#工作目录就是当前目录
our $BIN_DIR=$WORKING_DIR."/bin";#所有可执行文件所在的目录

##################
## 程序参数处理 ##
##################
&GetOptions( 'input=s' => \$input, 
	'strand_specific!' => \$strand_specific,
	'min_overlap=i' => \$min_overlap,
	'max_end_clip=i' => \$max_end_clip,
	'cpu_num=i' => \$cpu_num ,
	'mis_penalty=i' => \$mis_penalty,
	'gap_cost=i' => \$gap_cost,
	'gap_extension=i' => \$gap_extension											 
			 );
unless ($input) {#这个参数通过输入得到
die $usage;}
			 
#################
##  主程序开始 ##
#################
			 
my (@all_data, $name, $sequence); 
$name = ''; $sequence = ''; 

open(IN, "$input");
while (<IN>) {
	chomp; 
	s/\s+$//; 
	if (/^>(\S+)/) {
		my $tk = $1; 
		if ($name ne '') {#如果已经读入了序列名称（即前一次的），这样可以排除第一次（前一次为空）情况
			$sequence =~ s/\s//g; #去掉sequence中的空字符
			push(@all_data, [$name, length($sequence), $sequence]);
		}
		$name = $tk; #开始读入一条新序列，首先读入序列名称
		$sequence = ''; #并清空序列
	}else{
		$sequence .= $_; 
	}
}
close(IN);

if ($name ne '') { #最后一次读入的fasta数据，存入数组
	$sequence =~ s/\s//g;
	push(@all_data, [$name, length($sequence), $sequence]); 
}
my %inset;  #用于存储非冗余序列，这样会再次存储@all_data中的序列，应该改进为只保存index
my $restset = ''; #用于存储冗余序列，这些数据的保存仅为了校对用
@all_data = sort { -1*($a->[1] <=> $b->[1]) } @all_data; #@all_data中数据按照序列长度降序排序

my $contig_count=1; 
for my $tr (@all_data) {
	if (scalar(keys %inset)  == 0) {#第一次把最长的序列先放入inset集合
		$inset{$tr->[0]} = $tr->[2]; 
	}else{
		my $query = ">$tr->[0]\n$tr->[2]\n"; 
		my $return_string = &ifRedundant(\%inset, \$query);
		my @return_col = split(/\t/, $return_string);#第1列是name，第2列是r或n，或是合并后序列
		if ($return_col[1] eq "r") {#如果query是完全冗余，扔掉
			$restset .= $query;
		}elsif($return_col[1] eq "n"){#如果query是非冗余，包括进来
			$contig_count++;
			$inset{$tr->[0]} = $tr->[2]; 
		}
		else{#如果query是部分冗余，合并后替换
			my @names = split(/\:/, $return_col[0]);#第一列是(hit_name:query_name)		
			$inset{$names[0]} = $return_col[1];; #原来hit_name对应的记录被新序列覆盖
			$restset .= $query;
		}
	}
}

open(OUT1, ">inset");#将所有非冗余序列输出
while (my ($k,$v) = each %inset) {
	print OUT1 ">$k\n$v\n"; 
}
close(OUT1); 

open(OUT2, ">restset");#将所有冗余序列输出
print OUT2 $restset; 
close(OUT2);

print "@@\t".$input."\t".$contig_count."\n";


#######################
##     子程序开始    ##
#######################
sub ifRedundant {
	my ($inset, $query) = @_; 
	
	#保存所有被查询序列到文件tem
	open(OUT, ">tem");
	while (my ($k,$v) = each %$inset) {
	print OUT ">$k\n$v\n"; 
	}
	close(OUT);
	
	#保存查询序列到文件query
	open(OUT, ">query");
	print OUT $$query; 
	close(OUT);
	
    #执行blast，如果需要调试输出就用process_cmd()
	system("$BIN_DIR/formatdb -i tem -p F");#tem是inset文件里面的所有序列	
	my $blast_program = $BIN_DIR."/megablast";
	#my $blast_param = "-i query -d tem -o tem.paired -F F -a 8 -W $word_size -q -1 -G 2 -E 1 -b 5";
	my $blast_param = "-i query -d tem -o tem.paired -F $filter_query -a $cpu_num -W $word_size -q $mis_penalty -G $gap_cost -E $gap_extension -b $hits_return";
	if ($strand_specific)
	{
		$blast_param .= " -S 1";
	}
	system($blast_program." ".$blast_param);
	
	#从blast结果中找到redundant信息
	my $result = findRedundancy($inset, $query);
#   if($$query =~ />NOVEL1\n/){#调试用
#	    print STDERR $result."good\n";
#		die "this is >NOVEL1";
#	}
	return $result;
}

sub findRedundancy
{
	my ($inset, $query) = @_;	
	my ($query_name, $query_length, $hit_name, $hit_length, $hsp_length, $identity, $evalue, $score, $strand, $query_start, $query_end, $hit_start, $hit_end, $query_to_end, $hit_to_end);
	my %hsp = ();  #存储提取一个query的所有hsp，最后一起处理
	my $hsp_count=0; 
	my $query_sequenc; 
	my $subject_sequenc; 	
	my $is_hsp = 1;

	my $bfh = IO::File->new("tem.paired") || "Can not open blast result file: tem.paired $!\n";
	while(<$bfh>)
	{
		chomp;
		if (/Query=\s(\S+)/ || eof)#一旦遇到(新)Query Name或者文件最末（防止只有一个Query），就输出前面一个Query所有的hsp
		{
			#########################################
			#           记录前一个hsp记录           #
			#########################################
			if ($hsp_count > 0)#如果不是第一次出现（即已存有hsp），则保存（也可以输出）前面一个hsp结果
			{
				my $hsp_info =  $query_name."\t".$query_length."\t".$hit_name."\t".$hit_length."\t".
						$hsp_length."\t".$identity."\t".$evalue."\t".$score."\t".$strand."\t".
						$query_start."\t".$query_end."\t".$hit_start."\t".$hit_end;
				$hsp{$hsp_count} = $hsp_info;
			}

			#########################################
			#  分析属于上一对query-hit的所有hsp	#
			#########################################
			if (scalar(keys(%hsp)) > 0)
			{
				foreach my $hsp_num (sort {$a<=>$b} keys %hsp)
				{
					my @one_hit = split(/\t/, $hsp{$hsp_num});
					unless(scalar(@one_hit) == 13) { die "Error! the parsed blast may lost some info:\n $hsp{$hsp_num} \n"; }
					
					my ($query_name, $query_length, $hit_name, $hit_length, $hsp_length, $identity, $evalue, $score, $strand,
					    $query_start, $query_end, $hit_start, $hit_end) =
					   ($one_hit[0], $one_hit[1], $one_hit[2], $one_hit[3], $one_hit[4], $one_hit[5], $one_hit[6], $one_hit[7],
					    $one_hit[8], $one_hit[9], $one_hit[10],$one_hit[11],$one_hit[12]);

					my $query_to_end = $query_length - $query_end;
					my $hit_to_end = $hit_length - $hit_end;
					my $hit_to_start = $hit_length - $hit_start;

					$identity =~ s/%//;
					if($hsp_length <= 50) {$min_identity = 95;}
					elsif($hsp_length > 50 && $hsp_length <= 100) {$min_identity = 96;}
					else{$min_identity = 97;}
                    
					#下面判断的顺序非常重要
					if ($identity < $min_identity)#hsp的identity不够，不能合并，不能判定为非冗余，需要看下一个
					{
						next;#这样保证了，必须满足最小identity，才去判断下面条件，否则就跳过去了
					}
					if ($query_start -1 <= $max_end_clip  && $query_to_end <= $max_end_clip)#query被subject包括，判定冗余
					{
					    #my $hit_seq = $inset->{$hit_name}; 
					    #print "type1\t".$hit_name."\t".$query_name."\t".$hit_seq."\n";#输出合并信息供人工校对，调试用
						return $hit_name."\tr";
					}
					if($hsp_length >= $min_overlap)#第三种情况(有overlap)需要合并
					{

					    my $combined_seq;
					    (my $query_seq = $$query) =~ s/^>[^\n]+\n//;#替换掉前面的query name
					        $query_seq =~ s/\s//g;
						my $hit_seq = $inset->{$hit_name}; 
						if ($strand==1)
						{
						    #下列是query在前的情况
							if($query_start -1 > $max_end_clip  && $query_to_end <= $max_end_clip && $hit_start <= $max_end_clip)
							{ 
						       my $query_string = substr($query_seq, 0, $query_start); 							   
							   my $hit_string = substr($hit_seq, $hit_start, $hit_to_start);
							   #print "type2\t".$hit_name."\t".$query_name."\t".$query_string."\t".$hit_string."\n";#输出合并信息供人工校对，调试用
							   $combined_seq = $hit_name.":".$query_name."\t".$query_string.$hit_string;
							   return $combined_seq;
							} 
							#下列是hit在前的情况
							if($query_start -1 <= $max_end_clip  && $query_to_end > $max_end_clip && $hit_to_end <= $max_end_clip)
							{ 
							   my $hit_string = substr($hit_seq, 0, $hit_end); 
							   my $query_string = substr($query_seq, $query_end, $query_to_end); 
							   #print "type3\t".$hit_name."\t".$query_name."\t".$hit_string."\t".$query_string."\n";#输出合并信息供人工校对，调试用
							   $combined_seq = $hit_name.":".$query_name."\t".$hit_string.$query_string;
							   return $combined_seq;
							} 
						}
						if ($strand==-1)
						{
							#下列是query在前的情况
							if($query_start -1 > $max_end_clip  && $query_to_end <= $max_end_clip && $hit_to_end <= $max_end_clip)
							{ 
							   my $query_string = substr($query_seq, 0, $query_start); 							   
							   my $hit_string = substr($hit_seq, 0, $hit_end-1);
							   rcSeq(\$hit_string, 'rc'); #求序列的反向互补
							   #print "type4\t".$hit_name."\t".$query_name."\t".$query_string."\t".$hit_string."\n";#输出合并信息供人工校对，调试用
							   $combined_seq = $hit_name.":".$query_name."\t".$query_string.$hit_string;
							   return $combined_seq;
							} 
							#下列是hit在前的情况
							if($query_start -1 <= $max_end_clip  && $query_to_end > $max_end_clip && $hit_start-1 <= $max_end_clip)
							{ 
							   my $hit_string = substr($hit_seq, $hit_start, $hit_to_start); 
							   rcSeq(\$hit_string, 'rc'); #求序列的反向互补
							   my $query_string = substr($query_seq, $query_end, $query_to_end); 
							   #print "type5\t".$hit_name."\t".$query_name."\t".$hit_string."\t".$query_string."\n";#输出合并信息供人工校对，调试用
							   $combined_seq = $hit_name.":".$query_name."\t".$hit_string.$query_string;
							   return $combined_seq;
							} 
						}
					}#第三种情况需要合并
				}#循环结束，没有遇到符合冗余的判断，判定为非冗余
				return $hit_name."\tn";
			}
			
			#####################################
			#  开始记录一个新的query	    #
			#####################################
			%hsp = ();$hsp_count = 0;
			$query_name = $1; $query_length = ""; $hit_name = ""; $hit_length = "";
		}
		elsif (/\s+\((\S+)\sletters\)/)#取Query Length
		{
			$query_length = $1;
			$query_length =~ s/,//g;
		}
		
		elsif (/>(\S+)/)#一旦遇到Hit Name
		{
			#########################################
			#           记录前一个hsp记录           #
			#########################################
			if ($hsp_count > 0 || eof)#如果不是在Query后第一次出现（即已存有hsp），则保存前面一个hsp结果（也可以输出）
			{
				my $hsp_info =  $query_name."\t".$query_length."\t".$hit_name."\t".$hit_length."\t".
						$hsp_length."\t".$identity."\t".$evalue."\t".$score."\t".$strand."\t".
						$query_start."\t".$query_end."\t".$hit_start."\t".$hit_end;
				$hsp{$hsp_count} = $hsp_info;	
				$is_hsp = 0;
			}
			#################################
			#  开始记录一个新的hit	        #
			#################################
		    $hit_name = $1; $hit_length = "";
		}
		elsif (/\s+Length = (\d+)/)
		{
				$hit_length = $1;
				$hit_length =~ s/,//g;
		}

		elsif (/Score =\s+(\S+) bits.+Expect(\(\d+\))? = (\S+)/)#一旦遇到hsp
		{
			if ($hsp_count > 0 && $is_hsp == 1)
			{	
				my $hsp_info =  $query_name."\t".$query_length."\t".$hit_name."\t".$hit_length."\t".
								$hsp_length."\t".$identity."\t".$evalue."\t".$score."\t".$strand."\t".
								$query_start."\t".$query_end."\t".$hit_start."\t".$hit_end;
				$hsp{$hsp_count} = $hsp_info;
			}

			#################################
			#  开始记录一个新的hsp		#
			#################################
			$is_hsp = 1;
			$hsp_count++; 
			$score = $1; $evalue = $3;
			$evalue = "1$evalue" if ($evalue =~ m/^e/);
			$query_start = 0; $query_end = 0; $hit_start = 0; $hit_end = 0;

		}
		elsif (/\s+Identities = (\d+)\/(\d+)\s+\((\S+)\)/ && $hsp_count >= 1)
		{
			$identity = $1/$2*100;
			$identity = sprintf("%.".(2)."f", $identity);
			if ( $1 > $2 ) { $hsp_length = $1; } else { $hsp_length = $2; }
		}

		elsif (/\s+Strand = (\S+) \/ (\S+)/ && $hsp_count >= 1) 
		{
			if ( $2 eq "Plus" ) { $strand = 1;} else { $strand = -1;}
		}

		elsif (/Query\:\s(\d+)\s+(\S+)\s(\d+)/ && $hsp_count >= 1)
		{
			if ($query_start == 0) { $query_start = $1; $query_start =~ s/,//g;}
			$query_end = $3;
			$query_end =~ s/,//g;
		}

		elsif (/Sbjct\:\s(\d+)\s+(\S+)\s(\d+)/ && $hsp_count >= 1) 
		{
			if ( $strand == -1 )#永远保证$hit_start>=$hit_end
			{
				if ($hit_end == 0) { $hit_end = $1; $hit_end =~ s/,//g;  };
				$hit_start = $3;
				$hit_start =~ s/,//g;
			}
			else
			{
				if ($hit_start == 0) { $hit_start = $1; $hit_start =~ s/,//g; };
				$hit_end = $3;
				$hit_end =~ s/,//g;
			}
		}
		else
		{
			next;
		}
	}
	$bfh->close;
	return "null\tn";
}

sub rcSeq {
        my $seq_r = shift;
        my $tag = shift; defined $tag or $tag = 'rc'; # $tag = lc($tag);
        my ($Is_r, $Is_c) = (0)x2;
        $tag =~ /r/i and $Is_r = 1;
        $tag =~ /c/i and $Is_c = 1;
        #$tag eq 'rc' and ( ($Is_r,$Is_c) = (1)x2 );
        #$tag eq 'r' and $Is_r = 1;
        #$tag eq 'c' and $Is_c = 1;
        !$Is_r and !$Is_c and die "Wrong Input for function rcSeq! $!\n";
        $Is_r and $$seq_r = reverse ($$seq_r);
        # $Is_c and $$seq_r =~ tr/acgturyksbdhvnACGTURYKSBDHVN/tgcaayrmwvhdbnTGCAAYRMWVHDBN/;  # 2007-07-18 refer to NCBI;
        $Is_c and $$seq_r =~ tr/acgturykmbvdhACGTURYKMBVDH/tgcaayrmkvbhdTGCAAYRMKVBHD/; # edit on 2010-11-14;
        return 0;
}

sub process_cmd {
	my ($cmd) = @_;	
	print "CMD: $cmd\n";
	my $ret = system($cmd);	
	if ($ret) {
		print "Error, cmd: $cmd died with ret $ret";
	}
	return($ret);
}