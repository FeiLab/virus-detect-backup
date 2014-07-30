#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;
my $usage = <<_EOUSAGE_;

#########################################################################################
# virus_itentify.pl --file_list <FILE> --contig_type <String> --file_type <String> --reference [FILE]
#                   --diff_ratio --word_size [INT] --exp_value <Float> --identity_percen <Float>
#                   --cpu_num  [INT] --mis_penalty [INT] --gap_cost[INT] --gap_extension [INT]
#                                  
# Required(2):
#  --file_list The name of a txt file containing a list of input file names without any suffix
#  --contig_type The type of contig files(aligned��assembled or combined) 
#
# Options(3):
#  --file_type The format of files containing reads are fasta or fastq  [fastq]
#  --reference The name of a fasta file containing all of the virus reference sequences  [vrl_genbank.fasta] 
#  --diff_ratio The hits with distance less than 0.2 will be combined into one  [0.2] 
# 
# blast-related options(7):
#  --word_size [11] 
#  --exp_value [1e-5]
#  --identity_percen [80] 
#  --cpu_num         Number of processors to use [8] 
#  --mis_penalty     Penalty for a nucleotide mismatch [-1]
#  --gap_cost        Cost to open a gap [2] 
#  --gap_extension   Cost to extend a gap [1] 
#
###########################################################################################

_EOUSAGE_

	;

################################
##   ��������Ŀ¼���ļ���·�� ##
################################
our $WORKING_DIR=cwd();#����Ŀ¼���ǵ�ǰĿ¼
our $DATABASE_DIR=$WORKING_DIR."/databases";#�������ݿ��ļ����ڵ�Ŀ¼
our $BIN_DIR=$WORKING_DIR."/bin";#���п�ִ���ļ����ڵ�Ŀ¼
our $seq_info= $DATABASE_DIR."/vrl_genbank.info";
our $result_dir= $WORKING_DIR."/result";

###############################
##   ȫ�ֱ�������Ĭ�����ã�  ##
###############################
our $file_list;#�������д�����������ļ������ݣ����б��ļ����޺�׺��
our $contig_type;#����contig������
our $reference= "vrl_genbank.fasta";#����ȫ���ο����е��ļ����ƣ�FASTA��ʽ��
our $file_type= "fastq";

our $diff_ratio= 0.25;
our $word_size = 11;
our $cpu_num = 8;          #megablastʹ�õ�cpu��Ŀ
our $mis_penalty = -1;     #megablast�У��Դ���ķ��֣������Ǹ�����
our $gap_cost = 2;         #megablast�У���gap open�ķ��֣�������������
our $gap_extension = 1;    #megablast�У���gap open�ķ��֣�������������
our $exp_value = 1e-5;    #
our $identity_percen = 25;    #tblastx�Ե������������ȶ�ʱhsp����Сͬһ��

our $filter_query = "F";     #Ĭ�ϲ���Ҫȥ�������У��������Ҫ�û��趨
our $hits_return = 500;   #megablast���ص�hit��Ŀ���������Ҫ�û��趨
our $input_suffix='clean';           #�����ļ���׺��clean��unmapped��
#################
## �����������##
#################
&GetOptions( 'file_list=s' => \$file_list,
	'contig_type=s' => \$contig_type,
	'file_type=s' => \$file_type,
	'reference=s' => \$reference,
	'diff_ratio=f' => \$diff_ratio,		
	'word_size=i' => \$word_size,
	'exp_value=f' => \$exp_value,
	'identity_percen=f' => \$identity_percen,
	'cpu_num=i' => \$cpu_num,
	'mis_penalty=i' => \$mis_penalty,
	'gap_cost=i' => \$gap_cost,
	'gap_extension=i' => \$gap_extension
	 );

unless ($file_list) {#������Ҫ1������
die $usage;
}

#################
##  ������ʼ ##
#################
main: {
system("rm *.finished");
system("rm *.result");
my $sample;
my $i=0;
open(IN,$file_list) || die "Can't open the file $file_list\n";
my $format="-q";# Ĭ�������ļ�����fastq��ʽ
while(<IN>){
	#ÿ��ѭ������һ�У��������붼�Ǵ���������ļ��������޺�׺����
	$i=$i+1;
	chomp;
	$sample=$_; 
	print "#processing sample $i by $0: $sample\n";
	#ÿ����������һ���Լ��Ľ���ļ���
	my $sample_dir=$result_dir."_$sample";
	system("mkdir $sample_dir");
	&process_cmd("cp ./$contig_type/$sample.$contig_type.fa $sample_dir/contig_sequences.fa");#������copy������ļ���
		
	system("$BIN_DIR/formatdb -i $DATABASE_DIR/$reference -p F") unless (-e "$DATABASE_DIR/$reference.nhr");
	my $blast_program = $BIN_DIR."/megablast";
	my $blast_param = "-i ./$contig_type/$sample.$contig_type.fa -d $DATABASE_DIR/$reference -o $sample.blastn.paired -F $filter_query -a $cpu_num -W $word_size -q $mis_penalty -G $gap_cost -E $gap_extension -b $hits_return -e $exp_value";				
	&process_cmd($blast_program." ".$blast_param);
	&process_cmd("$BIN_DIR/blast_parse_table22.pl $sample.blastn.paired $sample.blastn.table1");#����blast���
	&process_cmd("$BIN_DIR/blast_filter3.pl $sample.blastn.table1 $sample.blastn.table 0.75 5");#��ÿ��query��hit�ԣ�ֻѡȡeֵ��ߵĽ��


	die;
	#system("rm $sample.blastn.paired");#�����꣬����ɾ��
	#system("rm $sample.blastn.table1");#�����꣬����ɾ��
	#s
	#known��novel contig����
	#&process_cmd("$BIN_DIR/query_filter1.pl $sample.blastn.table ./$contig_type/$sample.$contig_type.fa $sample.novel.contigs 60 50 > $sample.known.contigs");#known��novel contigs�ֿ�


	# split known and novel contigs
	
	&process_cmd("$BIN_DIR/query_filter1.pl $sample.blastn.table ./$contig_type/$sample.$contig_type.fa $sample.novel.contigs 80 25 > $sample.known.contigs");



	my $file_size= -s "$sample.known.contigs";#���ݴ��ļ���С�ǲ���0���������洦������
	if($file_size==0){#����ļ���С��0���ͽ�������ѭ��
		system("rm $sample.known.contigs");#�����������ܹ�align����֪�������contigs��table��ʽ��
		system("rm $sample.novel.contigs");#���������в���align����֪�������contigs��fasta��ʽ��
		system("rm $sample.blastn.table");
		system("touch $sample_dir/no_virus_detected");#��������ļ�����ʾ��������û�м�⵽����
		next;
	}
	#�ȿ���û��known��virus
	&process_cmd("$BIN_DIR/uniqComb.pl $sample.blastn.table -index $sample.known.contigs -col 0 -newCol 0 -exist > $sample.known.table");#known contigs��Ӧ��blast���
	&process_cmd("$BIN_DIR/hit_cov1.pl $sample.known.table $sample.known.cov 60 0.5 > $sample.known.block");
	&process_cmd("cut -f1 $sample.known.cov > hit_virus.list");
	&process_cmd("$BIN_DIR/extractFromFasta.pl -i $DATABASE_DIR/$reference --type list --query hit_virus.list --output1 hit_virus.fa --output2 remainding.fa");
	
	#�ٵõ�average depth��Ϣ
	my $sample_reads= $sample.".$input_suffix";#read�ļ�����Ҫalign����contigs��
	&process_cmd("bowtie-build --quiet hit_virus.fa virus");
	if($file_type eq "fasta"){$format="-f"};	
	&process_cmd("bowtie --quiet virus -v 1 -p 8 -a --best --strata $format $sample_reads -S --sam-nohead $sample.sam");
	&process_cmd("$BIN_DIR/samtools faidx hit_virus.fa");	
	&process_cmd("$BIN_DIR/samtools view -bt hit_virus.fa.fai $sample.sam > $sample.bam");
	&process_cmd("$BIN_DIR/samtools sort $sample.bam $sample.sorted");
	&process_cmd("$BIN_DIR/samtools mpileup -f hit_virus.fa $sample.sorted.bam > $sample.pileup");	
	&process_cmd("$BIN_DIR/pileup_depth.pl $sample.pileup $sample.known.depth");
	&process_cmd("$BIN_DIR/ColLink.pl $sample.known.cov -keyC1 0 -Col1 1,2 -keyC2 0 -add -f1 $sample.known.depth > 1.tem");#��depth�ļ���ȡ��2��3�мӵ��ļ�$sample.known.cov�ĺ���	
	&process_cmd("$BIN_DIR/ColLink.pl 1.tem -keyC1 0 -Col1 2,3 -keyC2 0 -add -f1 $seq_info > $sample.known.identified1");#��$seq_info��ȡ��3��4�мӵ��ļ�1.tem�ĺ���		
	&process_cmd("$BIN_DIR/hit_filter2.pl --input $sample.known.identified1 --diff_ratio $diff_ratio --output $sample.known.identified");#������������contig���ظ�������ϲ�hit	
	&process_cmd("$BIN_DIR/query_filter2.pl $sample.known.identified $sample.known.table $sample.contigs.table");#��contig��Ӧ��hit���ˣ����߶�Ӧ��known.identified��hit�ϣ����߶�Ӧ��e value��ߵ�hit����		
	&process_cmd("$BIN_DIR/ColLink.pl $sample.contigs.table -keyC1 0 -Col1 2,3 -keyC2 2 -add -f1 $seq_info > $sample.contigs.info");#���$seq_info�ļ��еĵ�3��4�У�genus��spicies��
	&process_cmd("$BIN_DIR/ColLink.pl $sample.contigs.info -keyC1 0 -Col1 1 -keyC2 0 -add -f1 $sample.known.contigs > $sample.contigs.table1");#���$sample.known.contigs�ļ��еĵ�2�У�contig���У�
	&process_cmd("$BIN_DIR/arrange_col1.pl $sample.contigs.table1 > $sample_dir/$sample.known.xls");#���¶Ա�񲼾֣��õ�known contigs�ıȶ��ļ�
	
	&process_cmd("cut -f1 $sample.known.identified > known.references.list");#ȡ��hit name����3�У�,����������sam�ļ�
	&process_cmd("$BIN_DIR/uniqComb.pl $sample.contigs.table1 -index known.references.list -col 0 -newCol 2 -exist > $sample.known.table1");#ֻҪ.known.identified�е�hit��Ӧ��
	&process_cmd("$BIN_DIR/extractFromFasta.pl -i $DATABASE_DIR/$reference --type list --query known.references.list --output1 $sample_dir/known.references.fa --output2 remainding.fa");#��ʾsam��	
	&process_cmd("$BIN_DIR/blastTable2sam.pl $sample.known.table1 > $sample_dir/$sample.known.sam");#������Ӧhit��sam�ļ�	
	
	system("mkdir $sample_dir/known_references");#����һ����Ŀ¼�����ڷ�����html�ļ�
	&process_cmd("$BIN_DIR/hit_cov1.pl $sample.known.table1 $sample.known.cov 60 0.5 > $sample.known.block");#��������$sample.known.cov
	&process_cmd("$BIN_DIR/ColLink.pl $sample.known.cov -keyC1 0 -Col1 6,7,8,9 -keyC2 0 -add -f1 $sample.known.identified > $sample.known.identified1");#��������$sample.known.identified1
	&process_cmd("$BIN_DIR/arrange_col2.pl $sample.known.identified1 > $sample.known.identified");#��������$sample.known.identified����ʽ�б仯��
	&process_cmd("$BIN_DIR/plot_results.pl $sample.known.identified $sample.known.table1 result_$sample known");

	#system("rm $sample.blastn.table");
	system("rm $sample.known.table");
	system("rm $sample.known.cov");
	system("rm $sample.known.block");#����ļ�ֻ����У��
	system("rm $sample.sam");
	system("rm $sample.bam");
	system("rm $sample.sorted.bam");
	system("rm $sample.pileup");
	system("rm $sample.known.depth");
	system("rm $sample.known.identified1");
	system("rm $sample.known.identified");
	system("rm $sample.contigs.table");
	system("rm $sample.contigs.info");
	system("rm $sample.contigs.table1");
	system("rm known.references.list");
	system("rm $sample.known.table1");
	system("rm *.ebwt");

	if($contig_type ne "aligned")#�����Ҫ��һ�����novel�Ĳ���
	{
		#���Ȱ�novel��contigs��tblastx�Աȵ�������
		$blast_program = $BIN_DIR."/blastall -p tblastx";
		$blast_param = "-i $sample.novel.contigs -d $DATABASE_DIR/$reference -o $sample.novel.paired -F $filter_query -a $cpu_num -e $exp_value";				
		&process_cmd($blast_program." ".$blast_param);
		&process_cmd("$BIN_DIR/blast_parse_table44.pl $sample.novel.paired $sample.novel1.table");
		&process_cmd("$BIN_DIR/blast_filter5.pl $sample.novel1.table $sample.novel.table 0 5");#��ÿ��query��hit�ԣ�ֻѡȡeֵ��ߵĽ��		
		system("rm $sample.novel.paired");
		system("rm $sample.novel1.table");
		my $file_size= -s "$sample.novel.table";#���ݴ��ļ���С�ǲ���0���������洦������
		if($file_size==0){#����ļ���С��0���ͽ�������ѭ��
			system("rm $sample.novel.table");
			system("rm $sample.novel.contigs");			
			system("touch $sample_dir/no_novel_virus_detected");#��������ļ�����ʾ��������û�м���²���
			system("cut -f1 $sample.known.contigs > known.contigs.list");
			&process_cmd("$BIN_DIR/extractFromFasta.pl -i $sample_dir/contig_sequences.fa --type list --query known.contigs.list --output1 known.contigs.fa --output2 $sample_dir/unknown.contigs.fa");
			system("rm known.contigs.list");
			system("rm $sample.known.contigs");
			system("rm known.contigs.fa");
			next;
		}
		
		#����ע���hsp��Ҫ��Ҫ��
		&process_cmd("$BIN_DIR/hit_cov1.pl $sample.novel.table $sample.novel.cov $identity_percen 0 > $sample.novel.block");
		&process_cmd("cut -f1 $sample.novel.cov > hit_virus.list");
		&process_cmd("$BIN_DIR/extractFromFasta.pl -i $DATABASE_DIR/$reference --type list --query hit_virus.list --output1 hit_virus.fa --output2 remainding.fa");
		
		#�ٵõ�average depth��Ϣ
		my $sample_reads= $sample.".$input_suffix";#read�ļ�����Ҫalign����contigs��
		&process_cmd("bowtie-build --quiet hit_virus.fa virus");				
		&process_cmd("bowtie --quiet virus -v 1 -p 8 -a --best --strata $format $sample_reads -S --sam-nohead $sample.sam");
		&process_cmd("$BIN_DIR/samtools faidx hit_virus.fa");	
		&process_cmd("$BIN_DIR/samtools view -bt hit_virus.fa.fai $sample.sam > $sample.bam");
		&process_cmd("$BIN_DIR/samtools sort $sample.bam $sample.sorted");
		&process_cmd("$BIN_DIR/samtools mpileup -f hit_virus.fa $sample.sorted.bam > $sample.pileup");	
		&process_cmd("$BIN_DIR/pileup_depth.pl $sample.pileup $sample.novel.depth");
		&process_cmd("$BIN_DIR/ColLink.pl $sample.novel.cov -keyC1 0 -Col1 1,2 -keyC2 0 -add -f1 $sample.novel.depth > 1.tem");#��depth�ļ���ȡ��2��3�мӵ��ļ�$sample.novel.cov�ĺ���
		&process_cmd("$BIN_DIR/ColLink.pl 1.tem -keyC1 0 -Col1 2,3 -keyC2 0 -add -f1 $seq_info > $sample.novel.identified1");#��$seq_info��ȡ��3��4�мӵ��ļ�1.tem�ĺ���			
		&process_cmd("$BIN_DIR/hit_filter2.pl --input $sample.novel.identified1 --diff_ratio $diff_ratio --output $sample.novel.identified");#���ǵõ��������ļ�
		&process_cmd("$BIN_DIR/query_filter2.pl $sample.novel.identified $sample.novel.table $sample.contigs.table");#�ѵ�һ�ļ��д��ڵ�hit��Ӧ��contig��������	
		&process_cmd("$BIN_DIR/ColLink.pl $sample.contigs.table -keyC1 0 -Col1 2,3 -keyC2 2 -add -f1 $seq_info > $sample.contigs.info");#���hit��ע����Ϣ
		&process_cmd("$BIN_DIR/fasta2tab.pl $sample.novel.contigs $sample.novel.contigs1");#��ʽת��
		&process_cmd("$BIN_DIR/ColLink.pl $sample.contigs.info -keyC1 0 -Col1 1 -keyC2 0 -add -f1 $sample.novel.contigs1 > $sample.contigs.table1");#���contig��������Ϣ		
		&process_cmd("$BIN_DIR/arrange_col1.pl $sample.contigs.table1 > $sample_dir/$sample.novel.xls");#���¶Ա�񲼾�
		
		&process_cmd("cut -f1 $sample.novel.identified > novel.references.list");#ȡ��hit name����3�У�,����������sam�ļ�
		&process_cmd("$BIN_DIR/uniqComb.pl $sample.contigs.table1 -index novel.references.list -col 0 -newCol 2 -exist > $sample.novel.table1");#ֻҪ.known.identified�е�hit��Ӧ��
	
		system("mkdir $sample_dir/novel_references");#����һ����Ŀ¼�����ڷ�����html�ļ�		
		&process_cmd("$BIN_DIR/hit_cov2.pl $sample.novel.table1 $sample.novel.cov $identity_percen 0 > $sample.novel.block");#��������$sample.novel.cov
		&process_cmd("$BIN_DIR/ColLink.pl $sample.novel.cov -keyC1 0 -Col1 6,7,8,9 -keyC2 0 -add -f1 $sample.novel.identified > $sample.novel.identified1");#��������$sample.novel.identified1
		&process_cmd("$BIN_DIR/arrange_col2.pl $sample.novel.identified1 > $sample.novel.identified");#��������$sample.novel.identified����ʽ�б仯��
		&process_cmd("$BIN_DIR/plot_results.pl $sample.novel.identified $sample.novel.table1 result_$sample novel");
		
		#������ȡ�Ȳ���knownҲ����novel��contigs
		system("cut -f1 $sample.known.contigs > known.contigs.list");
		system("cut -f1 $sample.novel.table | sort | uniq > novel.contigs.list");
		system("cat known.contigs.list novel.contigs.list > final.contigs.list");#������hit������
		&process_cmd("$BIN_DIR/extractFromFasta.pl -i $sample_dir/contig_sequences.fa --type list --query final.contigs.list --output1 known.contigs.fa --output2 $sample_dir/unknown.contigs.fa");
		
		system("rm $sample.novel.table");
		system("rm $sample.novel.cov");
		system("rm $sample.novel.block");#����ļ�ֻ����У��
		system("rm $sample.sam");
		system("rm $sample.bam");
		system("rm $sample.sorted.bam");
		system("rm $sample.pileup");
		system("rm $sample.novel.depth");
		system("rm $sample.novel.identified1");
		system("rm $sample.novel.identified");
		system("rm $sample.contigs.table");#���ǹ�known������ͬ���ļ�
		system("rm $sample.contigs.info");#���ǹ�known������ͬ���ļ�
		system("rm $sample.contigs.table1");#���ǹ�known������ͬ���ļ�
		system("rm novel.references.list");
		system("rm $sample.novel.table1");
		system("rm *.ebwt");
		system("rm $sample.known.contigs");#�������ɾ��
		system("rm $sample.novel.contigs");
		system("rm $sample.novel.contigs1");	
		system("rm known.contigs.list");
		system("rm novel.contigs.list");
		system("rm final.contigs.list");
		system("rm known.contigs.fa");
	}
}
close(IN);
system("rm hit_virus.list");
system("rm hit_virus.fa");
system("rm hit_virus.fa.fai");
system("rm remainding.fa");
system("rm *.tem");
system("rm *.log");
print "###############################\n";
print "All the samples have been processed by $0\n";
system("touch virus_identify.$contig_type.finished");#��������ļ�����ʾ������־

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
