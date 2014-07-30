#!/usr/bin/perl
#use strict;   #变量必须定义才可以使用
#use warnings; #对未初始化变量是给出警告
use Getopt::Long; #这个模块可以接受完整参数

#参数处理
my $usage = <<_EOUSAGE_;

#########################################################################################
# fastq_clipper.pl --filelist <FILE> --suffix <String> --adapter <String> --adapterLen INT --distance <INT> --levenshtein
# Required:
#  --filelist       a file containing all fastq file name for processing without suffix
#  --adapter        the reverse complimentary of 3' end adapter(5'->3')
#  --adapterLen     the adapter offset from the 5' end
#  --distance       the largest edit distance can be allowed between a adapter and reads
#  --levenshtein       using levenshtein distance, default is hamming distance
# Options(1):
#  --suffix	    the suffix of the input files
##########################################################################################

_EOUSAGE_
	;
#定义变量接受命令行参数	
our $filelist;        #包括需要处理的所有fastq文件名称列表
our $suffix="fastq";          #输入文件（fastq格式）的后缀名
our $adapter;         #adapter序列,所有的数据都是同一个adapter
our $adapter_offset;  #原始adapter序列的前adapter_offset个bp，否则就用adapter全长用来匹配
our $distance;        #adapter与read匹配时可允许的最大编辑距离，目前最多是1
our $levenshtein=0;#如果设定，此值自动转为1，表示使用levenshtein算法求编辑距离


#从命令行参数向变量传值
&GetOptions( 'filelist=s' => \$filelist,
             'adapter=s' => \$adapter,
	         'suffix=s' => \$suffix,
	         'adapterLen=i' => \$adapter_offset,
             'distance=i' => \$distance,
	         'levenshtein!' => \$levenshtein
			 );
unless ($filelist&&$adapter) {
	die $usage;
}
print "distance is $distance\n";

#定义全局变量
our $adapter_substr = substr($adapter, 0, $adapter_offset);#取出实际用于匹配的adapter部分，从第1个字符开始，共得到$adapter_offset个字符

#主程序开始

open(IN1,$filelist) || die "Can't open the sample list file\n";
while(<IN1>){
	chomp;	
	#print $levenshtein."\n";#调试专用
	trim_fastq($_);	
}
close(IN1);

print "all the files have been processed\n";
#主程序结束

#子程序结束
###########################################################
#	Process fastq file                                #
###########################################################

sub trim_fastq{
	my $basename =shift;        #得到第1个参数，样本文件名称
	my $input_file =$basename.".$suffix";                  #输入文件的名称
	my $output_file1 =$basename.".trimmed".$distance;      #输出文件1的名称
    my $output_file2 =$basename.".unmatched".$distance;    #输出文件2的名称
	my $output_file3 =$basename.".null".$distance;      #输出文件3的名称
	
	my $match_position;     #@match_result的第2个元素，barcode匹配到read上的的位置，
	my $id1;                #fastq文件中每个read的第1行
	my $DNAseq;             #fastq文件中每个read的第2行
	my $id2;                #fastq文件中每个read的第3行
	my $i = 0;              #用于读fastq文件循环内部，i表示fastq文件的第i行
	
    	open(IN2, $input_file) || die "Can't open the fastq file\n";
	open(OUT1, ">$output_file1") || die "can't create $output_file1 $!\n";#文件句柄与文件名映射
	open(OUT2, ">$output_file2") || die "can't create $output_file2 $!\n";#文件句柄与文件名映射
	open(OUT3, ">$output_file3") || die "can't create $output_file3 $!\n";#文件句柄与文件名映射
	while(<IN2>){
		chomp;
		$i ++;
		if($i%4 == 1){    #如果当前行是ID行
			$id1 = $_;
		}
		if($i%4 == 2){    #如果当前行是DNA序列,注意只有这里接受返回值
			$DNAseq = $_; #提取该DNA序列


			$match_position = adapter_match($DNAseq,$adapter_substr);#DNA序列与barcode序列匹配

			#print "$DNAseq,$adapter_substr,$match_position\n";

		}
		if($i%4 == 3){    #如果当前行是另一个ID行
			$id2 = $_;
		}
		if($i%4 == 0){    #如果当前行是质量序列（4的整数倍），根据匹配结果，开始输出		  

			if ($match_position ==0){#表示整条DNA都是污染，保存到".null"文件
				print OUT3 $id1."\n".$DNAseq."\n".$id2."\n".$_."\n";#输出整条序列
			}
			
			elsif($match_position ==length($DNAseq)){#没有找到匹配，保存到".unmatched"文件
				print OUT2 $id1."\n".$DNAseq."\n".$id2."\n".$_."\n";#输出整条序列
			}	
			
			else{#trim掉一部分，保存到".trimmed"文件
				my $trimmed_DNAseq = substr($DNAseq, 0, $match_position); #保留DNA序列上，barcode匹配点之前的部分
		        my $trimmed_qual = substr($_, 0, $match_position);        #保留质量序列上，barcode匹配点之前的部分
				print OUT1 $id1."\n".$trimmed_DNAseq."\n".$id2."\n".$trimmed_qual."\n";#输出trimmed后剩下的序列
			}		
		}#输出结束	
	}
	close(OUT1);
	close(OUT2);
	close(OUT3);
	close(IN2);
	my $trimmed_reads= count_reads($output_file1);            #统计并显示fastq文件的reads个数
	my $untrimmed_reads=count_reads($output_file2);           #统计并显示trim后文件的reads个数，不含空行
	my $null_reads=count_reads($output_file3);
	print $basename."\t".$null_reads."\t".$trimmed_reads."\t".$untrimmed_reads."\n";#显示文件名等信息
	
}

sub adapter_match{
	my $read_seq = shift;    #得到第1个参数，read序列
	my $adapter_seq = shift; #得到第2个参数，adapter序列
	my $read_len = length($read_seq);	
	my $position = $read_len;	#返回值初始值为全长，表示没有找到匹配
	my $min_offset= length($adapter_seq);
 
	if($levenshtein==1){#如果采用levenshtein距离
		for(my $i = 0; $i< $read_len-$adapter_offset+1; $i++){ 
			my $read_substr = substr($read_seq, $i, $adapter_offset);#检查过，没问题
			my $editDistance = levenshtein($read_substr,$adapter_seq);

			if ($editDistance <= $distance){
				$position=$i;  #记下这个位置，返回为了下一步trim
				last;          #找到了，就不继续找了，跳出循环	
			}				
		}			
	}
	else{#如果采用hamming距离
		for(my $i = 0; $i< $read_len-$min_offset+1; $i++){            
			my $read_substr = substr($read_seq, $i, $adapter_offset);#注意如果剩下的长度不够$adapter_offset，有多少截取多少而不报错
			my $editDistance = hamming($read_substr,$adapter_seq);#$read_substr有可能短于$adapter_offset，短的必须在前面

			if ($editDistance <= $distance){
				$position=$i;  #记下这个位置，返回为了下一步trim
				last;          #找到了，就不继续找了，跳出循环	
			}			
		}
	}		
	return $position;
}
##############################################
#    统计并显示fastq文件中的read数量         #
##############################################
sub count_reads{
        my $file = shift;
        my $count = `wc -l < $file`; #得到文件$file的行数
        die "wc failed: $?" if $?;
        chomp($count); 
        my $read = (int($count))/4;  #行数除以4，得到read数量
        return $read;
}
###############################################
#    求2个字符串之间的levenshtein距离         #
###############################################
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
            $mat{$i}{$j} = min($mat{$i-1}{$j} + 1,
                                $mat{$i}{$j-1} + 1,
                                $mat{$i-1}{$j-1} + $cost);
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
###########################################
#    求2个字符串之间的hamming距离         #
###########################################
sub hamming($$) { length( $_[ 0 ] ) - ( ( $_[ 0 ] ^ $_[ 1 ] ) =~ tr[\0][\0] ) }
