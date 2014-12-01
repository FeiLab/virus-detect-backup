#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;
use File::Basename;
my $usage = <<_EOUSAGE_;

#########################################################################################
# paired-bowtie2pileup.pl --file_list <FILE> --reference [FILE]  --out_type <String> --mismatchs [INT]
#                  --len_seed [INT] --dist_seed [INT] --quality_sum [INT] --thread_num [INT]                
# Required:
#  --file_list The name of a txt file containing a list of input file names without any suffix
#  --reference The name of a fasta file containing the reference sequences
#  --out_type  output mapped, sam, bam or pileup files
#
# Bowtie-related options(6): 
#  --mismatchs     alignments may have no more than v mismatches(v<=3), ignoring qualities
#  --len_seed      Take the first INT subsequence as seed [36] 
#  --dist_seed     Maximum edit distance in the seed [1] 
#  --quality_sum   Max sum of mismatch quals across alignment in the seed region [70]  
#  --thread_num    Number of threads (multi-threading mode) [8]  
##########################################################################################

_EOUSAGE_
	;
###############################
##   ȫ�ֱ�������Ĭ�����ã�  ##
###############################
our $file_list;#�������д�����������ļ������ݣ����б��ļ����޺�׺��
our $reference;#����ȫ���ο����е��ļ����ƣ�FASTA��ʽ��
our $index_name;#�ο����е���������
our $out_type;  #����ļ�����
our $mismatchs; #ȫ�ִ���
our $len_seed = 36; #����������
our $dist_seed = 1; #��������������༭����
our $quality_sum = 70; #���������д���λ�������Phred qualityֵ֮������
our $thread_num = 8; #������õ��߳����� 

################################
##   ��������Ŀ¼���ļ���·�� ##
################################
our $WORKING_DIR=cwd();#����Ŀ¼���ǵ�ǰĿ¼
our $BIN_DIR=$WORKING_DIR."/bin";#���п�ִ���ļ����ڵ�Ŀ¼

#################
## �����������##
#################
&GetOptions( 'file_list=s' => \$file_list,#�������д�������������ı��ļ����ƣ��޺�׺��
		'reference=s' => \$reference,#�ο�ת¼����ļ����ƣ�FASTA��ʽ��
		'out_type=s' => \$out_type,
		'mismatchs=i' => \$mismatchs,
		'len_seed=i' => \$len_seed,
		'dist_seed=i' => \$dist_seed,	
		'quality_sum=i' => \$quality_sum,		
		'thread_num=i' => \$thread_num
			 );

unless ($file_list&&$reference&&$out_type) {#������Ҫ2������
	die $usage;
}
$index_name = basename($reference);#ȥ��Ŀ¼���ƣ�ֻ�����ļ�����
$index_name =~ s/\.\S*$//;#ȥ���ļ���׺��

#################
##  ������ʼ ##
#################
main: {
    #����bowtieΪ�ο����н�������,������һ����alignment
    &process_cmd("$BIN_DIR/bowtie-build -f $reference $index_name") unless (-e "$index_name.ebwt");#���������ļ���aligment�ã����������Ҳ��ɾ��   
    &process_cmd("$BIN_DIR/samtools faidx $reference") unless (-e "$reference.fai");#���������ļ�������bam��ʽת���ã����������Ҳ��ɾ��      
    my $sample;
    my $i=0;
    open(IN, "$file_list");

    while (<IN>) {
		$i=$i+1;
		chomp;
		my @paired_seq = split(/\t/, $_);
		$sample="a"; #ÿ��ѭ������һ�У��������붼�Ǵ���������ļ��������޺�׺����
		print "#processing sample $i by $0: $sample\n";
		
		#aligment -> sam -> bam -> sorted bam -> pileup
		if($mismatchs){
			&process_cmd("$BIN_DIR/bowtie $index_name -v $mismatchs -p $thread_num -I 50 -X 600 -1 paired_seq[1].clean -2 paired_seq[2].clean -S -a $sample.pre.sam") unless (-s "$sample.pre.sam");		

		}
		else{
			&process_cmd("$BIN_DIR/bowtie $index_name -n $dist_seed -l $len_seed -e $quality_sum -p $thread_num -I 50 -X 600 -1 paired_seq[1].clean -2 paired_seq[2].clean -S -a $sample.pre.sam") unless (-s "$sample.pre.sam");	
		}
		&process_cmd("$BIN_DIR/SAM_filter_out_unmapped_reads.pl $sample.pre.sam $sample.unmapped $sample.mapped > $sample.sam") unless (-s "$sample.sam");	
		if($out_type eq "mapped"){next;}
		&process_cmd("samtools view -bt $reference.fai $sample.sam > $sample.bam") unless (-s "$sample.bam");
		if($out_type eq "sam"){next;}
		&process_cmd("samtools sort $sample.bam $sample.sorted") unless (-s "$sample.sorted.bam");
		if($out_type eq "bam"){next;}
		&process_cmd("$BIN_DIR/samtools sort $sample.bam $sample.sorted") unless (-s "$sample.sorted.bam");
		&process_cmd("$BIN_DIR/samtools mpileup -f $reference $sample.sorted.bam > $sample.pileup") unless (-s "$sample.pileup");
		if($out_type eq "pileup"){next;}
		my $file_size= -s "$sample.pileup";#����pileup�ļ���С�ǲ���0���������洦������
		if($file_size!=0){#����ļ���С����0��
			&process_cmd("java -cp $BIN_DIR extractConsensus $sample 0 40 1");#��ȡ����Ƭ�Σ�depth>=1,length>=40�����ļ�������1			
		}
	}
	close(IN);
	print "###############################\n";
	print "All the input files have been processed by the $0\n";
}
system("rm *.ebwt");
####
sub process_cmd {
	my ($cmd) = @_;	
	print "CMD: $cmd\n";
	my $ret = system($cmd);	
	if ($ret) {
		die "Error, cmd: $cmd died with ret $ret";
	}
	return($ret);
}
