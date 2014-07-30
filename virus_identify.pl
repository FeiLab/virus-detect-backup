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
#  --contig_type The type of contig files(aligned、assembled or combined) 
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
##   设置所有目录和文件的路径 ##
################################
our $WORKING_DIR=cwd();#工作目录就是当前目录
our $DATABASE_DIR=$WORKING_DIR."/databases";#所有数据库文件所在的目录
our $BIN_DIR=$WORKING_DIR."/bin";#所有可执行文件所在的目录
our $seq_info= $DATABASE_DIR."/vrl_genbank.info";
our $result_dir= $WORKING_DIR."/result";

###############################
##   全局变量（含默认设置）  ##
###############################
our $file_list;#包括所有待处理的输入文件（数据）的列表文件（无后缀）
our $contig_type;#输入contig的类型
our $reference= "vrl_genbank.fasta";#包括全部参考序列的文件名称（FASTA格式）
our $file_type= "fastq";

our $diff_ratio= 0.25;
our $word_size = 11;
our $cpu_num = 8;          #megablast使用的cpu数目
our $mis_penalty = -1;     #megablast中，对错配的罚分，必须是负整数
our $gap_cost = 2;         #megablast中，对gap open的罚分，必须是正整数
our $gap_extension = 1;    #megablast中，对gap open的罚分，必须是正整数
our $exp_value = 1e-5;    #
our $identity_percen = 25;    #tblastx以蛋白质序列来比对时hsp的最小同一性

our $filter_query = "F";     #默认不需要去除简单序列，这个不需要用户设定
our $hits_return = 500;   #megablast返回的hit数目，这个不需要用户设定
our $input_suffix='clean';           #数据文件后缀（clean或unmapped）
#################
## 输入参数处理##
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

unless ($file_list) {#至少需要1个参数
die $usage;
}

#################
##  主程序开始 ##
#################
main: {
system("rm *.finished");
system("rm *.result");
my $sample;
my $i=0;
open(IN,$file_list) || die "Can't open the file $file_list\n";
my $format="-q";# 默认数据文件都是fastq格式
while(<IN>){
	#每次循环读入一行，后续代码都是处理该样本文件（名称无后缀）。
	$i=$i+1;
	chomp;
	$sample=$_; 
	print "#processing sample $i by $0: $sample\n";
	#每个样本都有一个自己的结果文件夹
	my $sample_dir=$result_dir."_$sample";
	system("mkdir $sample_dir");
	&process_cmd("cp ./$contig_type/$sample.$contig_type.fa $sample_dir/contig_sequences.fa");#把数据copy到结果文件夹
		
	system("$BIN_DIR/formatdb -i $DATABASE_DIR/$reference -p F") unless (-e "$DATABASE_DIR/$reference.nhr");
	my $blast_program = $BIN_DIR."/megablast";
	my $blast_param = "-i ./$contig_type/$sample.$contig_type.fa -d $DATABASE_DIR/$reference -o $sample.blastn.paired -F $filter_query -a $cpu_num -W $word_size -q $mis_penalty -G $gap_cost -E $gap_extension -b $hits_return -e $exp_value";				
	&process_cmd($blast_program." ".$blast_param);
	&process_cmd("$BIN_DIR/blast_parse_table22.pl $sample.blastn.paired $sample.blastn.table1");#解析blast结果
	&process_cmd("$BIN_DIR/blast_filter3.pl $sample.blastn.table1 $sample.blastn.table 0.75 5");#对每个query和hit对，只选取e值最高的结果


	die;
	#system("rm $sample.blastn.paired");#解析完，立刻删除
	#system("rm $sample.blastn.table1");#解析完，立刻删除
	#s
	#known与novel contig分离
	#&process_cmd("$BIN_DIR/query_filter1.pl $sample.blastn.table ./$contig_type/$sample.$contig_type.fa $sample.novel.contigs 60 50 > $sample.known.contigs");#known和novel contigs分开


	# split known and novel contigs
	
	&process_cmd("$BIN_DIR/query_filter1.pl $sample.blastn.table ./$contig_type/$sample.$contig_type.fa $sample.novel.contigs 80 25 > $sample.known.contigs");



	my $file_size= -s "$sample.known.contigs";#根据此文件大小是不是0，进入下面处理流程
	if($file_size==0){#如果文件大小是0，就结束本次循环
		system("rm $sample.known.contigs");#保存了所有能够align到已知病毒库的contigs（table格式）
		system("rm $sample.novel.contigs");#保存了所有不能align到已知病毒库的contigs（fasta格式）
		system("rm $sample.blastn.table");
		system("touch $sample_dir/no_virus_detected");#建立这个文件，表示此样本中没有检测到病毒
		next;
	}
	#先看有没有known的virus
	&process_cmd("$BIN_DIR/uniqComb.pl $sample.blastn.table -index $sample.known.contigs -col 0 -newCol 0 -exist > $sample.known.table");#known contigs对应的blast结果
	&process_cmd("$BIN_DIR/hit_cov1.pl $sample.known.table $sample.known.cov 60 0.5 > $sample.known.block");
	&process_cmd("cut -f1 $sample.known.cov > hit_virus.list");
	&process_cmd("$BIN_DIR/extractFromFasta.pl -i $DATABASE_DIR/$reference --type list --query hit_virus.list --output1 hit_virus.fa --output2 remainding.fa");
	
	#再得到average depth信息
	my $sample_reads= $sample.".$input_suffix";#read文件，需要align到新contigs上
	&process_cmd("bowtie-build --quiet hit_virus.fa virus");
	if($file_type eq "fasta"){$format="-f"};	
	&process_cmd("bowtie --quiet virus -v 1 -p 8 -a --best --strata $format $sample_reads -S --sam-nohead $sample.sam");
	&process_cmd("$BIN_DIR/samtools faidx hit_virus.fa");	
	&process_cmd("$BIN_DIR/samtools view -bt hit_virus.fa.fai $sample.sam > $sample.bam");
	&process_cmd("$BIN_DIR/samtools sort $sample.bam $sample.sorted");
	&process_cmd("$BIN_DIR/samtools mpileup -f hit_virus.fa $sample.sorted.bam > $sample.pileup");	
	&process_cmd("$BIN_DIR/pileup_depth.pl $sample.pileup $sample.known.depth");
	&process_cmd("$BIN_DIR/ColLink.pl $sample.known.cov -keyC1 0 -Col1 1,2 -keyC2 0 -add -f1 $sample.known.depth > 1.tem");#从depth文件提取第2、3列加到文件$sample.known.cov的后面	
	&process_cmd("$BIN_DIR/ColLink.pl 1.tem -keyC1 0 -Col1 2,3 -keyC2 0 -add -f1 $seq_info > $sample.known.identified1");#从$seq_info提取第3、4列加到文件1.tem的后面		
	&process_cmd("$BIN_DIR/hit_filter2.pl --input $sample.known.identified1 --diff_ratio $diff_ratio --output $sample.known.identified");#根据所包括的contig的重复情况，合并hit	
	&process_cmd("$BIN_DIR/query_filter2.pl $sample.known.identified $sample.known.table $sample.contigs.table");#把contig对应的hit过滤，或者对应到known.identified的hit上，或者对应到e value最高的hit上面		
	&process_cmd("$BIN_DIR/ColLink.pl $sample.contigs.table -keyC1 0 -Col1 2,3 -keyC2 2 -add -f1 $seq_info > $sample.contigs.info");#填加$seq_info文件中的第3、4列（genus和spicies）
	&process_cmd("$BIN_DIR/ColLink.pl $sample.contigs.info -keyC1 0 -Col1 1 -keyC2 0 -add -f1 $sample.known.contigs > $sample.contigs.table1");#填加$sample.known.contigs文件中的第2列（contig序列）
	&process_cmd("$BIN_DIR/arrange_col1.pl $sample.contigs.table1 > $sample_dir/$sample.known.xls");#重新对表格布局，得到known contigs的比对文件
	
	&process_cmd("cut -f1 $sample.known.identified > known.references.list");#取出hit name（第3列）,将用于生成sam文件
	&process_cmd("$BIN_DIR/uniqComb.pl $sample.contigs.table1 -index known.references.list -col 0 -newCol 2 -exist > $sample.known.table1");#只要.known.identified中的hit对应的
	&process_cmd("$BIN_DIR/extractFromFasta.pl -i $DATABASE_DIR/$reference --type list --query known.references.list --output1 $sample_dir/known.references.fa --output2 remainding.fa");#显示sam用	
	&process_cmd("$BIN_DIR/blastTable2sam.pl $sample.known.table1 > $sample_dir/$sample.known.sam");#产生对应hit的sam文件	
	
	system("mkdir $sample_dir/known_references");#产生一个子目录，用于放最后的html文件
	&process_cmd("$BIN_DIR/hit_cov1.pl $sample.known.table1 $sample.known.cov 60 0.5 > $sample.known.block");#重新生成$sample.known.cov
	&process_cmd("$BIN_DIR/ColLink.pl $sample.known.cov -keyC1 0 -Col1 6,7,8,9 -keyC2 0 -add -f1 $sample.known.identified > $sample.known.identified1");#重新生成$sample.known.identified1
	&process_cmd("$BIN_DIR/arrange_col2.pl $sample.known.identified1 > $sample.known.identified");#重新生成$sample.known.identified（格式有变化）
	&process_cmd("$BIN_DIR/plot_results.pl $sample.known.identified $sample.known.table1 result_$sample known");

	#system("rm $sample.blastn.table");
	system("rm $sample.known.table");
	system("rm $sample.known.cov");
	system("rm $sample.known.block");#这个文件只用于校对
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

	if($contig_type ne "aligned")#如果需要进一步检查novel的病毒
	{
		#首先把novel的contigs用tblastx对比到病毒库
		$blast_program = $BIN_DIR."/blastall -p tblastx";
		$blast_param = "-i $sample.novel.contigs -d $DATABASE_DIR/$reference -o $sample.novel.paired -F $filter_query -a $cpu_num -e $exp_value";				
		&process_cmd($blast_program." ".$blast_param);
		&process_cmd("$BIN_DIR/blast_parse_table44.pl $sample.novel.paired $sample.novel1.table");
		&process_cmd("$BIN_DIR/blast_filter5.pl $sample.novel1.table $sample.novel.table 0 5");#对每个query和hit对，只选取e值最高的结果		
		system("rm $sample.novel.paired");
		system("rm $sample.novel1.table");
		my $file_size= -s "$sample.novel.table";#根据此文件大小是不是0，进入下面处理流程
		if($file_size==0){#如果文件大小是0，就结束本次循环
			system("rm $sample.novel.table");
			system("rm $sample.novel.contigs");			
			system("touch $sample_dir/no_novel_virus_detected");#建立这个文件，表示此样本中没有检测新病毒
			system("cut -f1 $sample.known.contigs > known.contigs.list");
			&process_cmd("$BIN_DIR/extractFromFasta.pl -i $sample_dir/contig_sequences.fa --type list --query known.contigs.list --output1 known.contigs.fa --output2 $sample_dir/unknown.contigs.fa");
			system("rm known.contigs.list");
			system("rm $sample.known.contigs");
			system("rm known.contigs.fa");
			next;
		}
		
		#这里注意对hsp的要求都要低
		&process_cmd("$BIN_DIR/hit_cov1.pl $sample.novel.table $sample.novel.cov $identity_percen 0 > $sample.novel.block");
		&process_cmd("cut -f1 $sample.novel.cov > hit_virus.list");
		&process_cmd("$BIN_DIR/extractFromFasta.pl -i $DATABASE_DIR/$reference --type list --query hit_virus.list --output1 hit_virus.fa --output2 remainding.fa");
		
		#再得到average depth信息
		my $sample_reads= $sample.".$input_suffix";#read文件，需要align到新contigs上
		&process_cmd("bowtie-build --quiet hit_virus.fa virus");				
		&process_cmd("bowtie --quiet virus -v 1 -p 8 -a --best --strata $format $sample_reads -S --sam-nohead $sample.sam");
		&process_cmd("$BIN_DIR/samtools faidx hit_virus.fa");	
		&process_cmd("$BIN_DIR/samtools view -bt hit_virus.fa.fai $sample.sam > $sample.bam");
		&process_cmd("$BIN_DIR/samtools sort $sample.bam $sample.sorted");
		&process_cmd("$BIN_DIR/samtools mpileup -f hit_virus.fa $sample.sorted.bam > $sample.pileup");	
		&process_cmd("$BIN_DIR/pileup_depth.pl $sample.pileup $sample.novel.depth");
		&process_cmd("$BIN_DIR/ColLink.pl $sample.novel.cov -keyC1 0 -Col1 1,2 -keyC2 0 -add -f1 $sample.novel.depth > 1.tem");#从depth文件提取第2、3列加到文件$sample.novel.cov的后面
		&process_cmd("$BIN_DIR/ColLink.pl 1.tem -keyC1 0 -Col1 2,3 -keyC2 0 -add -f1 $seq_info > $sample.novel.identified1");#从$seq_info提取第3、4列加到文件1.tem的后面			
		&process_cmd("$BIN_DIR/hit_filter2.pl --input $sample.novel.identified1 --diff_ratio $diff_ratio --output $sample.novel.identified");#这是得到的最终文件
		&process_cmd("$BIN_DIR/query_filter2.pl $sample.novel.identified $sample.novel.table $sample.contigs.table");#把第一文件中存在的hit对应的contig保留下来	
		&process_cmd("$BIN_DIR/ColLink.pl $sample.contigs.table -keyC1 0 -Col1 2,3 -keyC2 2 -add -f1 $seq_info > $sample.contigs.info");#添加hit的注释信息
		&process_cmd("$BIN_DIR/fasta2tab.pl $sample.novel.contigs $sample.novel.contigs1");#格式转换
		&process_cmd("$BIN_DIR/ColLink.pl $sample.contigs.info -keyC1 0 -Col1 1 -keyC2 0 -add -f1 $sample.novel.contigs1 > $sample.contigs.table1");#添加contig的序列信息		
		&process_cmd("$BIN_DIR/arrange_col1.pl $sample.contigs.table1 > $sample_dir/$sample.novel.xls");#重新对表格布局
		
		&process_cmd("cut -f1 $sample.novel.identified > novel.references.list");#取出hit name（第3列）,将用于生成sam文件
		&process_cmd("$BIN_DIR/uniqComb.pl $sample.contigs.table1 -index novel.references.list -col 0 -newCol 2 -exist > $sample.novel.table1");#只要.known.identified中的hit对应的
	
		system("mkdir $sample_dir/novel_references");#产生一个子目录，用于放最后的html文件		
		&process_cmd("$BIN_DIR/hit_cov2.pl $sample.novel.table1 $sample.novel.cov $identity_percen 0 > $sample.novel.block");#重新生成$sample.novel.cov
		&process_cmd("$BIN_DIR/ColLink.pl $sample.novel.cov -keyC1 0 -Col1 6,7,8,9 -keyC2 0 -add -f1 $sample.novel.identified > $sample.novel.identified1");#重新生成$sample.novel.identified1
		&process_cmd("$BIN_DIR/arrange_col2.pl $sample.novel.identified1 > $sample.novel.identified");#重新生成$sample.novel.identified（格式有变化）
		&process_cmd("$BIN_DIR/plot_results.pl $sample.novel.identified $sample.novel.table1 result_$sample novel");
		
		#下面提取既不是known也不是novel的contigs
		system("cut -f1 $sample.known.contigs > known.contigs.list");
		system("cut -f1 $sample.novel.table | sort | uniq > novel.contigs.list");
		system("cat known.contigs.list novel.contigs.list > final.contigs.list");#所有有hit的序列
		&process_cmd("$BIN_DIR/extractFromFasta.pl -i $sample_dir/contig_sequences.fa --type list --query final.contigs.list --output1 known.contigs.fa --output2 $sample_dir/unknown.contigs.fa");
		
		system("rm $sample.novel.table");
		system("rm $sample.novel.cov");
		system("rm $sample.novel.block");#这个文件只用于校对
		system("rm $sample.sam");
		system("rm $sample.bam");
		system("rm $sample.sorted.bam");
		system("rm $sample.pileup");
		system("rm $sample.novel.depth");
		system("rm $sample.novel.identified1");
		system("rm $sample.novel.identified");
		system("rm $sample.contigs.table");#覆盖过known产生的同名文件
		system("rm $sample.contigs.info");#覆盖过known产生的同名文件
		system("rm $sample.contigs.table1");#覆盖过known产生的同名文件
		system("rm novel.references.list");
		system("rm $sample.novel.table1");
		system("rm *.ebwt");
		system("rm $sample.known.contigs");#这里才能删除
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
system("touch virus_identify.$contig_type.finished");#建立这个文件，表示结束标志

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
