#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;
use FindBin;
use File::Basename;

my $usage = <<_EOUSAGE_;

#########################################################################################
# bwa-alignAndCorrect.pl --file_list <FILE> --reference <FILE> --coverage <Float> 
#                 		 --max_dist[INT] --max_open [INT] --max_extension [INT] --len_seed [INT] --dist_seed [INT] --thread_num [INT]
# 
# Required(3):
#  --file_list A txt file containing a list of input file names without any suffix
#  --reference A fasta file containing all the reference sequences
#  --coverage  A reference sequence covered by reads will be output as results
# 
# BWA-related options(6):
#  --max_dist      Maximum edit distance [1]  
#  --max_open      Maximum number of gap opens [1]  
#  --max_extension Maximum number of gap extensions [1]  
#  --len_seed      Take the first INT subsequence as seed [15] 
#  --dist_seed     Maximum edit distance in the seed [1]  
#  --thread_num    Number of threads (multi-threading mode) [8]
###########################################################################################

_EOUSAGE_

;
#################
#  global vars  #
#################
our $file_list;			# list of sample name
our $reference;			# reference sequence, fasta format
our $index_name;		# index name for reference
our $coverage;  		# cutoff of mapping reads / reference
our $max_dist = 1; 		# BWA allowed max distance 
our $max_open = 1; 		# BWA allowed max gap opens
our $max_extension = 1; # BWA allowed max gap extension
our $len_seed = 15; 	# length of seed
our $dist_seed = 1; 	# BWA allowed seed max distance
our $thread_num = 8; 	# thread number

################################
# set folder and file path     #
################################
our $WORKING_DIR	= cwd();								# current folder
our $DATABASE_DIR	= ${FindBin::RealBin}."/../databases";	# database folder
our $BIN_DIR		= ${FindBin::RealBin};					# programs folder
our $seq_info		= $DATABASE_DIR."/vrl_genbank.info";	# genbank info

###################
# get input paras #
###################
GetOptions( 
	'file_list=s'		=> \$file_list,
	'reference=s' 		=> \$reference,
	'coverage=f' 		=> \$coverage,
	'max_dist=i' 		=> \$max_dist,
	'max_open=i' 		=> \$max_open,
	'max_extension=i' 	=> \$max_extension,
	'len_seed=i' 		=> \$len_seed,
	'dist_seed=i' 		=> \$dist_seed,			 
	'thread_num=i' 		=> \$thread_num 
);

die $usage unless ($file_list && $reference && $coverage);	# required parameters
$index_name = basename($reference);							# remove folder name
#$index_name =~ s/\.\S*$//;									# remove suffix

#################
# main          #
#################
main: {
    # create bwa and fasta index file for reference file (plant virus)   
    process_cmd("$BIN_DIR/bwa index -p $DATABASE_DIR/$index_name -a bwtsw $reference 2> bwa.log") unless (-e "$DATABASE_DIR/$index_name.amb"); 
    process_cmd("$BIN_DIR/samtools faidx $reference") unless (-e "$reference.fai");
    
	# parse samples in file list 
	my $sample;
    my $i=0;
    open(IN, "$file_list") || die $!;
    while (<IN>) {
		$i=$i+1;
		chomp;
		$sample=$_;
		die "Error, file $sample does not exist\n" unless -s $sample;		
		print "#processing sample $i by $0: $sample\n";
		
		#aligment -> sam -> bam -> sorted bam -> pileup
		process_cmd("$BIN_DIR/bwa aln -n $max_dist -o $max_open -e $max_extension -i 0 -l $len_seed -k $dist_seed -t $thread_num $DATABASE_DIR/$index_name $sample 1> $sample.sai 2>> bwa.log") unless (-s "$sample.sai");
		process_cmd("$BIN_DIR/bwa samse -n 10000 -s $DATABASE_DIR/$index_name $sample.sai $sample 1> $sample.pre.sam 2>> bwa.log") unless (-s "$sample.pre.sam");			
		process_cmd("$BIN_DIR/SAM_filter_out_unmapped_reads.pl $sample.pre.sam > $sample.sam1") unless (-s "$sample.sam1");
		process_cmd("$BIN_DIR/samFilter.pl $sample.sam1 > $sample.sam") unless (-s "$sample.sam");	# only keep the best hits of reads alignment to reference
		process_cmd("$BIN_DIR/samtools view -bt $reference.fai $sample.sam > $sample.bam") unless (-s "$sample.bam");
		process_cmd("$BIN_DIR/samtools sort $sample.bam $sample.sorted") unless (-s "$sample.sorted.bam");
		process_cmd("$BIN_DIR/samtools mpileup -f $reference $sample.sorted.bam > $sample.pre.pileup") unless (-s "$sample.pre.pileup");

		#如果一个reference覆盖的区域不到其长度的某个比例，则去除这个reference所对应的所有行
		my $depth_cutoff;
		$depth_cutoff = pileup_filter("$sample.pre.pileup", "$seq_info", "$coverage", "$sample.pileup") unless (-s "$sample.pileup");
        
		my $file_size = -s "$sample.pileup";	# get the size of pileup file
		if($file_size != 0){
			&process_cmd("java -cp $BIN_DIR extractConsensus $sample 0 40 1");	#提取连续片段（depth>=1,length>=40），文件名含有1			
			my $result_contigs=$sample.".contigs1.fa";							#这个是corrected的序列
			my $count = `grep \'>\' $result_contigs | wc -l`;
			print "@".$sample."\t".$count;										#统计得到多少序列，java程序运行后，得到的序列有可能为0
			if($count!=0){ 														#如果$sample.contigs1.fa的序列数量不为0，改名称后转移到aligned文件夹
				system("$BIN_DIR/renameFasta.pl --inputfile $result_contigs --outputfile $sample.contigs2.fa --prefix ALIGNED");#每条序列的名字要统一命名
				system("mv $sample.contigs2.fa ./aligned/$sample.contigs1.fa");	#把结果文件移动到aligned文件夹
				system("rm $result_contigs");
			}
			else{																#如果$sample.contigs1.fa的序列数量为0，直接转移到aligned文件夹
				system("mv $result_contigs ./aligned/$sample.contigs1.fa");		#把结果文件移动到aligned文件夹	
			}
		}	
		else{																	#如果文件大小是0，或者$sample.pileup不存在
			print "@".$sample."\t0\n";											#得到0条序列
			system("touch ./aligned/$sample.contigs1.fa");						#建立一个空的结果文件			
		}			
		system("rm $sample.sai");
		system("rm $sample.pre.sam");
		system("rm $sample.sam1");
		system("rm $sample.sam");
		system("rm $sample.bam");
		system("rm $sample.sorted.bam");
		system("rm $sample.pre.pileup");
		system("rm $sample.pileup");	
		system("rm bwa.log");
	}
	close(IN);
	print "###############################\n";
	print "All the input files have been processed by $0\n";
	#system("touch bwa-alignAndCorrect.run.finished");							#建立这个文件，表示结束标志
}

#################
# subroutine    #
#################
sub process_cmd {
	my ($cmd) = @_;	
	print "CMD: $cmd\n";
	my $ret = system($cmd);	#成功就返回0，否则就失败
	if ($ret) {
		print "Error, cmd: $cmd died with ret $ret";
	}
	return($ret);
}
sub pileup_filter {
	my ($inputfile,$seq_info,$ratio,$outputfile) = @_;
	open(IN1, "$seq_info" );#读入此文件中序列的长度信息
	my %hash;
	my @a;
	my $id;
	while(my $line=<IN1>){#读入length文件
		chomp($line);
		@a= split(/\t/,$line);
		$id=$a[0];#第1列是序列id
		$hash{$id}=$a[1];#第2列是序列长度，存入一个hash表中
	}
	close IN1;

	open(IN2,  "$inputfile" );
	open(OUT, ">$outputfile" );

	my $pre_reference="a";
	my $current_reference;
	my $i=1;#从1开始计数，得到当前reference sequence的长度
	my @data;
	my $len_threshold=0;
	my $line;
	my $totalDepth=0;
	my $totalPositions=0;
	while (<IN2>) {
	    $line = $_;    
	    chomp($line);
	    my @a = split( /\t/, $line);
	    $current_reference= $a[0];#取当前的参考序列id
	    $totalDepth=$totalDepth+$a[3];#累加这个文件中所有位点的depth
	    $totalPositions=$totalPositions+1;#累加这个文件中所有位点的数量
	    if($current_reference eq $pre_reference){#如果当前的参考序列id与前面的相同，表示处于同一个参考序列内部
			$i++;#@a数组的指针加1，后面会加入对应的元素	
	    }
	    else{#id改变，开始判断并处理
		if($pre_reference ne "a"){#如果不是第1个query，第1个实际是空
			$len_threshold=$hash{$pre_reference}*$ratio;
			#print $pre_reference."\t".$len_threshold."\t".$i."\n";#输出上一个reference上覆盖碱基的实际长度，和长度阈值（总长*要求的覆盖度）
			if($i>=$len_threshold){#如果此序列的实际长度超过阈值
				foreach my $data (@data) {
					print OUT "$data\n";#输出到结果文件中
				}
			}
			$i=1;#重新开始计数
			@data=();  
		} 
	    }
	    push @data, $line;#每次循环都要存pileup文件的一行
	    $pre_reference=$current_reference;
	}
	print "the average depth of $inputfile is \t", 1.0*$totalDepth/$totalPositions, "\n";

	#别忘了处理剩下的部分内容
	$len_threshold=$hash{$current_reference}*$ratio;
	#print $current_reference."\t".$len_threshold."\t".$i."\n";
	if($i>=$len_threshold){
		foreach my $data (@data) {
		    print OUT "$data\n";
		}
	}
	close(IN2);
	close(OUT);
	my $aveDepth=1.0*$totalDepth/$totalPositions;
	return	$aveDepth;
}