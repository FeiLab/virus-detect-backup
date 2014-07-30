#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;
use IO::File;
use Bio::SeqIO;
use FindBin;
use File::Basename;

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
#  --diff_ratio The hits with distance less than 0.25 will be combined into one  [0.25] 
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
# New options(3):
#  --hsp_cover		The coverage of hsp should be more than this cutoff for query or hit [0.75]
#  --diff_contig_cover	The coverage for different contigs [0.5]
#  --diff_contig_length	The length of different contigs [100]
#
###########################################################################################

_EOUSAGE_

;

################################
## set folder and file path   ##
################################
my $WORKING_DIR = cwd();				# current folder : working folder
my $BIN_DIR = ${FindBin::RealBin};			# bin folder
my $result_dir = $WORKING_DIR."/result";		# result folder
my $tf = $WORKING_DIR."/temp";				# temp folder

my $DATABASE_DIR = ${FindBin::RealBin}."/../databases";	# database folder
my $seq_info = $DATABASE_DIR."/vrl_genbank.info";	# virus sequence info
###############################
##  global vars		     ##
###############################
my $file_list;			# temp list of input file
my $file_type= "fastq";
my $contig_type;		# aligned assembled conbined
my $reference= "";		# virus reference, fasta format

my $hsp_cover = 0.75;		# for blast filter
my $diff_ratio= 0.25;		# ratio for number of diff contig
my $diff_contig_cover = 0.5;	# for hit filter
my $diff_contig_length = 100;	# for hit filter
my $word_size = 11;
my $cpu_num = 8;		# megablast: thread number
my $mis_penalty = -1;		# megablast: penalty for mismatch
my $gap_cost = 2;		# megablast: penalty for gap open
my $gap_extension = 1;		# megablast: penalty for gap extension
my $exp_value = 1e-5;		#
my $identity_percen = 25;	# tblastx以蛋白质序列来比对时hsp的最小同一性

my $filter_query = "F";		# megablast: F - disable remove simple sequence
my $hits_return = 500;		# megablast: hit number
my $input_suffix='';

########################
# get input parameters #
########################
GetOptions( 
	'file_list=s'		=> \$file_list,
	'contig_type=s' 	=> \$contig_type,
	'file_type=s' 		=> \$file_type,
	'reference=s' 		=> \$reference,
	'diff_ratio=f' 		=> \$diff_ratio,
	'hsp_cover=f'		=> \$hsp_cover,
	'diff_contig_cover=f'	=> \$diff_contig_cover,
	'diff_contig_length=i'	=> \$diff_contig_length,
	'word_size=i' 		=> \$word_size,
	'exp_value=f' 		=> \$exp_value,
	'identity_percen=f' 	=> \$identity_percen,
	'cpu_num=i' 		=> \$cpu_num,
	'mis_penalty=i' 	=> \$mis_penalty,
	'gap_cost=i' 		=> \$gap_cost,
	'gap_extension=i' 	=> \$gap_extension
);

die $usage unless $file_list;	# required parameter

#################
#    main       #
#################
main: {

my $sample;
my $i=0;
my $format = "-q"; # default input file is fastq format

my $in =IO::File->new($file_list) || die "Can't open the list file $file_list $!\n";
while(<$in>)
{
	chomp;
	$i=$i+1; $sample=$_; 
	my $sample_base = basename($sample);
	print "#processing sample $i : $sample_base using $0\n";

	# create folder for each sample
	my $sample_dir = $result_dir."_".$sample_base;
	system("mkdir $sample_dir") unless -e $sample_dir;

	# copy aligned, assemblied, or combined contigs to each folder, named as contig_sequences
	my $contigs = $sample.".".$contig_type;
	process_cmd("cp $contigs $sample_dir/contig_sequences.fa");

	# format the reference plant virus sequence, blast it against contig sequences	
	system("$BIN_DIR/formatdb -i $DATABASE_DIR/$reference -p F") unless (-e "$DATABASE_DIR/$reference.nhr");
	my $blast_program = $BIN_DIR."/megablast";
	my $blast_param = "-i $contigs -d $DATABASE_DIR/$reference -o $sample.blastn.paired -F $filter_query -a $cpu_num -W $word_size -q $mis_penalty -G $gap_cost -E $gap_extension -b $hits_return -e $exp_value";				
	process_cmd($blast_program." ".$blast_param);

	# parse blast results to table format (NOT USE SEARCH::IO)
	# output table format:
	# query_name \t query_length \t hit_name \t hit_length \t hsp_length \t identity \t evalue \t score \t strand \t 
	# query_start \t query_end \t hit_start \t hit_end \t identity2 \t aligned_query \t aligned_hit \t aligned_string\n
	process_cmd("$BIN_DIR/blast_parse_table22.pl $sample.blastn.paired $sample.blastn.table1");

	# get paired query and hit according to e value	
	# the coverage of hsp should be more than 0.75 for query or hit
	# keep the best pair of <query,hit> according to evalue
	# 一个query对应的hits中 只保留比最低evalue那个hsp的identity少于drop_off(5%)以内的
	#process_cmd("$BIN_DIR/blast_filter3.pl $sample.blastn.table1 $sample.blastn.table 0.75 5");
	# *** do not understand why using the last parameter:1 for this function	
	blast_filter("$sample.blastn.table1", "$sample.blastn.table", $hsp_cover, 5, 1);
	
	# separate known and novel contig
	# 1. get the all the hsps with > 60 identity for one query
	# 2. combine all the hsps for one query
	# 3. get ratio of combined hsp and query
	# 4. separate know if ratio > 50 and novel contigs  
	process_cmd("$BIN_DIR/query_filter1.pl $sample.blastn.table $contigs $sample.novel.contigs 60 50 > $sample.known.contigs");

	my $file_size= -s "$sample.known.contigs";	#根据此文件大小是不是0，进入下面处理流程
	if($file_size==0){				#如果文件大小是0，就结束本次循环
		#system("rm $sample.known.contigs");	#保存了所有能够align到已知病毒库的contigs（table格式）
		#system("rm $sample.novel.contigs");	#保存了所有不能align到已知病毒库的contigs（fasta格式）
		#system("rm $sample.blastn.table");
		system("touch $sample_dir/no_virus_detected");#建立这个文件，表示此样本中没有检测到病毒
		next;
	}

	# get the blast result (tab delimit format ) for known contigs
	process_cmd("$BIN_DIR/uniqComb.pl $sample.blastn.table -index $sample.known.contigs -col 0 -newCol 0 -exist > $sample.known.table");

	# get coverage for each hit  
	#process_cmd("$BIN_DIR/hit_cov1.pl $sample.known.table $sample.known.cov 60 0.5 > $sample.known.block");
	hit_cov("$sample.known.table", "$sample.known.cov", "$sample.known.block", 60, 0.5, 1);

	# get hit virus sequence ID
	process_cmd("cut -f1 $sample.known.cov > $tf/hit_virus.list");

	# generate hit virus sequence and remainding sequence
	process_cmd("$BIN_DIR/extractFromFasta.pl -i $DATABASE_DIR/$reference --type list --query $tf/hit_virus.list --output1 $tf/hit_virus.fa --output2 $tf/remainding.fa");

	# align reads to hit virus to get average depth
	my $sample_reads= $sample;
	process_cmd("bowtie-build --quiet $tf/hit_virus.fa $tf/virus");

	if($file_type eq "fasta"){$format="-f"};	
	process_cmd("bowtie --quiet $tf/virus -v 1 -p 8 -a --best --strata $format $sample_reads -S --sam-nohead $sample.sam");
	process_cmd("$BIN_DIR/samtools faidx $tf/hit_virus.fa");	
	process_cmd("$BIN_DIR/samtools view -bt $tf/hit_virus.fa.fai $sample.sam > $sample.bam");
	process_cmd("$BIN_DIR/samtools sort $sample.bam $sample.sorted");
	process_cmd("$BIN_DIR/samtools mpileup -f $tf/hit_virus.fa $sample.sorted.bam > $sample.pileup");	
	process_cmd("$BIN_DIR/pileup_depth.pl $sample.pileup $sample.known.depth");

	# extract the column 2, and 3 form sample.know.depth, 
	# add them to file sample.known.cov, generate new file 1.tem	
	process_cmd("$BIN_DIR/ColLink.pl $sample.known.cov -keyC1 0 -Col1 1,2 -keyC2 0 -add -f1 $sample.known.depth > $tf/1.tem");

	# extract the column 3, and 4 from seq_info
	# add them to file 1.tem, generate new file sample.known.identified1
	process_cmd("$BIN_DIR/ColLink.pl $tf/1.tem -keyC1 0 -Col1 2,3 -keyC2 0 -add -f1 $seq_info > $sample.known.identified1");	

	# combine hit accoridng to contig ? how it works
	process_cmd("$BIN_DIR/hit_filter2.pl --input1 $sample.known.identified1 --input2 $sample.known.table --diff_ratio $diff_ratio --diff_contig_cover $diff_contig_cover --diff_contig_length $diff_contig_length --output $sample.known.identified");

	##把contig对应的hit过滤，或者对应到known.identified的hit上，或者对应到e value最高的hit上面	
	process_cmd("$BIN_DIR/query_filter2.pl $sample.known.identified $sample.known.table $sample.contigs.table");	

	#填加$seq_info文件中的第3、4列（genus和spicies）
	process_cmd("$BIN_DIR/ColLink.pl $sample.contigs.table -keyC1 0 -Col1 2,3 -keyC2 2 -add -f1 $seq_info > $sample.contigs.info");

	#填加$sample.known.contigs文件中的第2列（contig序列） 
	process_cmd("$BIN_DIR/ColLink.pl $sample.contigs.info -keyC1 0 -Col1 1 -keyC2 0 -add -f1 $sample.known.contigs > $sample.contigs.table1");

	#重新对表格布局，得到known contigs的比对文件
	arrange_col1("$sample.contigs.table1", "$sample_dir/$sample_base.known.xls");

	#取出hit name（第3列）,将用于生成sam文件
	process_cmd("cut -f1 $sample.known.identified > $tf/known.references.list");

	#只要.known.identified中的hit对应的
	process_cmd("$BIN_DIR/uniqComb.pl $sample.contigs.table1 -index $tf/known.references.list -col 0 -newCol 2 -exist > $sample.known.table1");

	# #显示sam用	
	process_cmd("$BIN_DIR/extractFromFasta.pl -i $DATABASE_DIR/$reference --type list --query $tf/known.references.list --output1 $sample_dir/known.references.fa --output2 $tf/remainding.fa");

	#产生对应hit的sam文件
	process_cmd("$BIN_DIR/blastTable2sam.pl $sample.known.table1 > $sample_dir/$sample_base.known.sam");

	# generate folder for contig alignment with html format	
	system("mkdir $sample_dir/known_references");

	#重新生成$sample.known.cov
	#process_cmd("$BIN_DIR/hit_cov1.pl $sample.known.table1 $sample.known.cov 60 0.5 > $sample.known.block");
	hit_cov("$sample.known.table1", "$sample.known.cov", "$sample.known.block", 60, 0.5, 1);

	#重新生成$sample.known.identified1
	process_cmd("$BIN_DIR/ColLink.pl $sample.known.cov -keyC1 0 -Col1 6,7,8,9 -keyC2 0 -add -f1 $sample.known.identified > $sample.known.identified1");

	#重新生成$sample.known.identified（格式有变化）
	arrange_col2("$sample.known.identified1",  "$sample.known.identified");
	process_cmd("$BIN_DIR/plot_results.pl $sample.known.identified $sample.known.table1 $sample_dir known");

	#system("rm $sample.blastn.table");
	#system("rm $sample.known.table");
	#system("rm $sample.known.cov");
	#system("rm $sample.known.block");#这个文件只用于校对
	#system("rm $sample.sam");
	#system("rm $sample.bam");
	#system("rm $sample.sorted.bam");
	#system("rm $sample.pileup");
	#system("rm $sample.known.depth");
	#system("rm $sample.known.identified1");
	#system("rm $sample.known.identified");
	#system("rm $sample.contigs.table");
	#system("rm $sample.contigs.info");
	#system("rm $sample.contigs.table1");
	#system("rm known.references.list");
	#system("rm $sample.known.table1");
	#system("rm *.ebwt");

	# check novel virus of contig type is not aligned 
	if($contig_type ne "aligned")
	{
		# compare noval contigs against virus database using tblastx
		$blast_program = $BIN_DIR."/blastall -p tblastx";
		$blast_param = "-i $sample.novel.contigs -d $DATABASE_DIR/$reference -o $sample.novel.paired -F $filter_query -a $cpu_num -e $exp_value";				
		process_cmd($blast_program." ".$blast_param);
		process_cmd("$BIN_DIR/blast_parse_table44.pl $sample.novel.paired $sample.novel1.table");


		# get best paired query and hit according to e value	
		# the coverage of hsp should be more than 0 for query or hit
		# keep the best pair of <query,hit> according to evalue
		# 一个query对应的hits中 只保留比最低evalue那个hsp的identity少于drop_off(5%)以内的
		#process_cmd("$BIN_DIR/blast_filter5.pl $sample.novel1.table $sample.novel.table 0 5");
		blast_filter("$sample.novel1.table", "$sample.novel.table", 0, 5, 3);
		#system("rm $sample.novel.paired");
		#system("rm $sample.novel1.table");

		# check the blast result
		# of the blast do not have any hit (file_size is 0), exit it
		my $file_size= -s "$sample.novel.table";
		if($file_size==0){
			#system("rm $sample.novel.table");
			#system("rm $sample.novel.contigs");			
			system("touch $sample_dir/no_novel_virus_detected"); # create this sign for no virus identified
			system("cut -f1 $sample.known.contigs > $tf/known.contigs.list");
			process_cmd("$BIN_DIR/extractFromFasta.pl -i $sample_dir/contig_sequences.fa --type list --query $tf/known.contigs.list --output1 $tf/known.contigs.fa --output2 $sample_dir/unknown.contigs.fa");
			#system("rm known.contigs.list");
			#system("rm $sample.known.contigs");
			#system("rm known.contigs.fa");
			next;
		}
		
		#这里注意对hsp的要求都要低
		# process_cmd("$BIN_DIR/hit_cov1.pl $sample.novel.table $sample.novel.cov $identity_percen 0 > $sample.novel.block");
		hit_cov("$sample.novel.table", "$sample.novel.cov", "$sample.novel.block", $identity_percen, 0, 1);
		
		process_cmd("cut -f1 $sample.novel.cov > $tf/hit_virus.list");
		process_cmd("$BIN_DIR/extractFromFasta.pl -i $DATABASE_DIR/$reference --type list --query $tf/hit_virus.list --output1 $tf/hit_virus.fa --output2 $tf/remainding.fa");
		
		#再得到average depth信息
		my $sample_reads= $sample.$input_suffix;#read文件，需要align到新contigs上
		process_cmd("bowtie-build --quiet $tf/hit_virus.fa $tf/virus");				
		process_cmd("bowtie --quiet $tf/virus -v 1 -p 8 -a --best --strata $format $sample_reads -S --sam-nohead $sample.sam");
		process_cmd("$BIN_DIR/samtools faidx $tf/hit_virus.fa");	
		process_cmd("$BIN_DIR/samtools view -bt $tf/hit_virus.fa.fai $sample.sam > $sample.bam");
		process_cmd("$BIN_DIR/samtools sort $sample.bam $sample.sorted");
		process_cmd("$BIN_DIR/samtools mpileup -f $tf/hit_virus.fa $sample.sorted.bam > $sample.pileup");	
		process_cmd("$BIN_DIR/pileup_depth.pl $sample.pileup $sample.novel.depth");


		#从depth文件提取第2、3列加到文件$sample.novel.cov的后面
		process_cmd("$BIN_DIR/ColLink.pl $sample.novel.cov -keyC1 0 -Col1 1,2 -keyC2 0 -add -f1 $sample.novel.depth > $tf/1.tem");

		#从$seq_info提取第3、4列加到文件1.tem的后面
		process_cmd("$BIN_DIR/ColLink.pl $tf/1.tem -keyC1 0 -Col1 2,3 -keyC2 0 -add -f1 $seq_info > $sample.novel.identified1");			
		#这是得到的最终文件
		process_cmd("$BIN_DIR/hit_filter2.pl --input1 $sample.novel.identified1 --input2 $sample.novel.table --diff_ratio $diff_ratio --diff_contig_cover $diff_contig_cover --diff_contig_length $diff_contig_length --output $sample.novel.identified");

		#把第一文件中存在的hit对应的contig保留下来
		process_cmd("$BIN_DIR/query_filter2.pl $sample.novel.identified $sample.novel.table $sample.contigs.table");

		#添加hit的注释信息	
		process_cmd("$BIN_DIR/ColLink.pl $sample.contigs.table -keyC1 0 -Col1 2,3 -keyC2 2 -add -f1 $seq_info > $sample.contigs.info");

		# convert fasta sequence to tab delimit format
		# column1: seqID
		# column2: sequence
		fasta2tab("$sample.novel.contigs", "$sample.novel.contigs1");

		#添加contig的序列信息
		process_cmd("$BIN_DIR/ColLink.pl $sample.contigs.info -keyC1 0 -Col1 1 -keyC2 0 -add -f1 $sample.novel.contigs1 > $sample.contigs.table1");		
		# re-generate table 
		#process_cmd("$BIN_DIR/arrange_col1.pl $sample.contigs.table1 > $sample_dir/$sample.novel.xls");
		arrange_col1("$sample.contigs.table1", "$sample_dir/$sample_base.novel.xls");	

		# get hit name (3rd col) for generate sam file	
		process_cmd("cut -f1 $sample.novel.identified > $tf/novel.references.list");

		#只要.known.identified中的hit对应的
		process_cmd("$BIN_DIR/uniqComb.pl $sample.contigs.table1 -index $tf/novel.references.list -col 0 -newCol 2 -exist > $sample.novel.table1");

		#产生一个子目录，用于放最后的html文件	
		system("mkdir $sample_dir/novel_references");

		#重新生成$sample.novel.cov
		# process_cmd("$BIN_DIR/hit_cov2.pl $sample.novel.table1 $sample.novel.cov $identity_percen 0 > $sample.novel.block");
		hit_cov("$sample.novel.table1", "$sample.novel.cov", "$sample.novel.block", $identity_percen, 0, 3);


		#重新生成$sample.novel.identified1
		process_cmd("$BIN_DIR/ColLink.pl $sample.novel.cov -keyC1 0 -Col1 6,7,8,9 -keyC2 0 -add -f1 $sample.novel.identified > $sample.novel.identified1");

		# re-generate sample.novel.identified (format changed)
		#process_cmd("$BIN_DIR/arrange_col2.pl $sample.novel.identified1 > $sample.novel.identified");#重新生成$sample.novel.identified（格式有变化）
		arrange_col2("$sample.novel.identified1", "$sample.novel.identified");
		process_cmd("$BIN_DIR/plot_results.pl $sample.novel.identified $sample.novel.table1 $sample_dir novel");
		
		#下面提取既不是known也不是novel的contigs
		system("cut -f1 $sample.known.contigs > $tf/known.contigs.list");
		system("cut -f1 $sample.novel.table | sort | uniq > $tf/novel.contigs.list");
		system("cat $tf/known.contigs.list $tf/novel.contigs.list > $tf/final.contigs.list");#所有有hit的序列
		process_cmd("$BIN_DIR/extractFromFasta.pl -i $sample_dir/contig_sequences.fa --type list --query $tf/final.contigs.list --output1 $tf/known.contigs.fa --output2 $sample_dir/unknown.contigs.fa");
		
		#system("rm $sample.novel.table");
		#system("rm $sample.novel.cov");
		#system("rm $sample.novel.block");#这个文件只用于校对
		#system("rm $sample.sam");
		#system("rm $sample.bam");
		#system("rm $sample.sorted.bam");
		#system("rm $sample.pileup");
		#system("rm $sample.novel.depth");
		#system("rm $sample.novel.identified1");
		#system("rm $sample.novel.identified");
		#system("rm $sample.contigs.table");#覆盖过known产生的同名文件
		#system("rm $sample.contigs.info");#覆盖过known产生的同名文件
		#system("rm $sample.contigs.table1");#覆盖过known产生的同名文件
		#system("rm novel.references.list");
		#system("rm $sample.novel.table1");
		#system("rm *.ebwt");
		#system("rm $sample.known.contigs");#这里才能删除
		#system("rm $sample.novel.contigs");
		#system("rm $sample.novel.contigs1");	
		#system("rm known.contigs.list");
		#system("rm novel.contigs.list");
		#system("rm final.contigs.list");
		#system("rm known.contigs.fa");
	}
}
$in->close;
print "###############################\n";
print "All the samples have been processed by $0\n";
}

# put folder new folder

#################################################################
# kentnf: subroutine						#
#################################################################
sub process_cmd 
{
	my ($cmd) = @_;	
	print "CMD: $cmd\n";
	my $ret = system($cmd);	
	if ($ret) {
		die "Error, cmd: $cmd died with ret $ret";
	}
	return($ret);
}


=head2
 blast_filter : 
 blast_filter(input, output, coverage, identity_dropoff, ???param??? );
=cut
sub blast_filter
{
	my ($input, $output, $coverage, $identity_dropoff, $param) = @_;

	my $in = IO::File->new($input) || die $!;
	my $out = IO::File->new(">".$output) || die $!;

	my $last_query="";
	my $last_hit="";
	my $current_query;
	my $current_hit;
	my $high_identity;
	my $current_identity;

	while (<$in>) 
	{
		my @cols = split(/\t/, $_);

		if( $param*$cols[4]/$cols[1] >= $coverage || $param*$cols[4]/$cols[3] >= $coverage)	#query或hit至少要被hsp覆盖一定比例
		{ 

			$current_query=$cols[0];
			$current_hit=$cols[2];	
			$current_identity=$cols[5];

			if($current_query ne $last_query){	#如果遇到一个新的query
				$high_identity=$cols[5];	#记录最高的identity
				print $out $_;			#query的第一个hit肯定要，因为是evalue最低的
			}else{
				if ($current_hit ne $last_hit && ($current_identity >= $high_identity-$identity_dropoff))#相同的<query,hit>中，只保留evalue最低的(即第一个)
				{
					print $out $_;	
				}	
			}
			$last_query=$cols[0];
			$last_hit=$cols[2];
		}
	}
	
	$in->close;
	$out->close;
}

=head2
 hit_cov : 
 usage: hit_cov(input, output1, output2, identify, query_cov, output2, ???param???);

 # input file is the output of blast_parse_table2.pl
 # format of input file:
 # query_name \t query_length \t hit_name \t hit_length \t hsp_length \t identity \t evalue \t score \t	
 # strand \t query_start \t query_end \t hit_start \t hit_end ... \n
 # index 2、3、4、5、11、12 are required
 # for each pair of [query, hit] should meet the requirement of cutoff_identity and query_cov
 # ratio = all query coverage / hit_length
 # hit: query1, query2, ... queryN

=cut
sub hit_cov
{
	my ($input, $output1, $output2, $cutoff_identity, $query_cov, $param) = @_;

	# put hit info to hash hash
	# %bkl
	# key: hit_name, 
	# value: arrays of block and query
	# 	 for each element in array is another array include three element [hit_start, hit_end, hash_of_query]
	#        
	#        % hash_of_query
	#        key: queryID
	#        value: 1
	#
 	# %hit_len
	# key: hit_name
	# value: hit_length
	my %blk;
	my %hit_len;

	my $in = IO::File->new($input) || die $!;
	my $out1 = IO::File->new(">".$output1) || die $!;
	my $out2 = IO::File->new(">".$output2) || die $!;
	while (<$in>) 
	{
		chomp; 
		my @ta = split(/\t/, $_);
		die "Error in line $_" if scalar @ta < 13;
		my ($query_name, $query_length,	$hit_name, $hit_length,	$hsp_length, $identity,	$evalue,
		$score, $strand, $query_start, $query_end, $hit_start, $hit_end) = (@ta[0..12]);

		# query covered = param * hsp length / query length
		my $query_covered= $param * $hsp_length / $query_length; 
		if ( $identity >= $cutoff_identity && $query_covered >= $query_cov)
		{
			#建立hit和(hit_start,hit_end)之间的映射
			push(@{$blk{$hit_name}}, [$hit_start, $hit_end]); 
			$blk{$ta[2]}[-1][2]{$ta[0]}=1; 
			defined $hit_len{$hit_name} or $hit_len{$hit_name} = $hit_length; 
		}
	}
	$in->close;

	#print OUT join("\t", qw/hit_name hit_len total_cov cov_rate% query_names/)."\n";
	print $out2 join("\t", qw/hit_name block_len block_start block_end covered_len query_names/)."\n"; 

	foreach my $tk (sort keys %blk) # sort by hit name 
	{
		# structure of array '@o'
		# [array] 		, [array]		  ... [array]
		# [start, end, %contigs], [start, end, %contigs], ... [start, end, %contigs]

		my @o; #存储没有overlap的hit上的block

		# sort by hit start, hit end, and query name
		foreach my $ar (sort {$a->[0]<=>$b->[0] || $a->[1]<=>$b->[1] || $a->[2] cmp $b->[2];} @{$blk{$tk}}) 
		{
			#if ($tk eq 'D10663')
			#{
			#	print "$$ar[0]\t$$ar[1]\n";
			#	my %sss = %{$$ar[2]};
			#	foreach my $k (sort keys %sss) {
			#		print $k."\n";
			#	}
			#	die;
			#}

			if (scalar(@o) == 0) {
				push(@o, [@$ar]); 
			}elsif ($o[-1][1] >= $ar->[0]-1) {	#把一个hit上重叠的区域合并
				for my $qname (keys %{$ar->[2]}) {
					$o[-1][2]{$qname}=1; 
				}
				$o[-1][1] < $ar->[1] and $o[-1][1] = $ar->[1]; 
			}else{
				push(@o, [@$ar]); 
			}
		}

		my $total_cov = 0; # hit totoal coverage 
		my %aa; 
		for my $ar (@o) {
			my @query_names = sort keys %{$ar->[2]};
			# output to file
			# hit_name \t hit_len \t hit_block_start \t hit_block_end \t hit_block_length \t query names
			print $out2 join("\t", $tk, $hit_len{$tk}, $ar->[0], $ar->[1], $ar->[1]-$ar->[0]+1, join(",", @query_names))."\n";
			
			$total_cov += ($ar->[1]-$ar->[0]+1); 		# 一个hit上所有非重叠block的长度之和
			@aa{@query_names} = (1) x scalar(@query_names); # 将所有contigs的名称存入hash表
		}

		# output file: known.cov
		# format:
		# hit_id \t hit_length \t hit_converage_len \t hit_coverage_% \t contigs_name \t contigs_num
		print $out1 join("\t", $tk, $hit_len{$tk}, $total_cov, $total_cov/$hit_len{$tk}*1.0, join(",", sort keys %aa), scalar(keys %aa))."\n"; 
		
	}
	$out1->close;
	$out2->close;
}

=head2
 query_filter1 : 
 query_filter1(inputfile1, inputfile2, output1, output2, identity, min_ratio);
 # inputFile1 is the output of blast_parse_table2.pl
 # format of inputFile1
 # query_name \t query_length \t hit_name \t hit_length \t hsp_length \t identity \t evalue \t score \t strand \t 
 # query_start \t query_end \t hit_start \t hit_end \n
 # the index 0, 1, 9, 10 are required
 # all hit must meet the cutoff_identity
 # ratio = all hsp of hit / query_length
 # 把query被覆盖超过一定ratio的记录下来，然后从一个fasta文件中提取记录中不包括的相应序列，输出到output1
 # 剩下的序列的名称输出到ouput2
 # 如果hit都是病毒，可以用这个ratio来判断所得到contig是不是病毒
=cut
=head
sub query_filter1
{
	my ($input1, $input2, $output1, $output2, $identity, $min_ratio) = @_;

	my %blk; 
	my %query_len; 
	my %query_filtered;
	
	my $in1 = IO::File->new($input1) || die $!;
	while(<$in1>)
	{	
		chomp; 
		my @ta = split(/\t/, $_); 
		if ($ta[5]>=$cutoff_identity){			# only use record meet the cutoff_identity
			push(@{$blk{$ta[0]}}, [@ta[9,10]]); 	# create map between queryName and (query_start,query_end)
			defined $query_len{$ta[0]} or $query_len{$ta[0]} = $ta[1];
		}
	}
	$in1->close;

	# print join("\t", qw/query_name block_len block_start block_end covered_len/)."\n"; #向标准输出输出的各列名称
	# print OUT join("\t", qw/query_name query_len total_cov cov_rate%/)."\n";   #向结果文件输出的各列名称

	foreach my $tk (sort keys %blk) # sort by query name，然后提取所有需要filtered掉的query名称
	{
		my @o; #存储没有overlap的query上的block
		for my $ar (sort { $a->[0]<=>$b->[0] || $a->[1]<=>$b->[1];} @{$blk{$tk}}) 
		{
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
		} else {
			if($flag eq "on") {#表示为"on"
				print $_;#后面的序列需要继续向标准输出输出
			} else {
				print OUT $_;#否则，后面的序列需要继续向OUT输出
			}
		}
	}
	close(IN2);
	close(OUT);
}
=cut

=head2
 fasta2tab: convert fasta sequence to tab delimit file
=cut
sub fasta2tab
{
	my ($input, $output) = @_;

	my $out = IO::File->new(">".$output) || die $!;
	my $in = Bio::SeqIO->new(-format=>'fasta', -file=>$input);
	while(my $inseq = $in->next_seq)
	{
		print $out $inseq->id."\t".$inseq->seq."\n";
	}
	$out->close;
}

=head2
 arrange_col1 : 

# 首先调换列次数，然后排序，最后加表头输出
# 下列column需要互相调换,先根据decription，后根据Hit_ID排序
# Contig_ID	Contig_length	Hit_ID	Hit_length	strand	Contig_start	Contig_end	Hit_start	Hit_end	identity	genus	 description	contig_sequence
=cut
sub arrange_col1 
{
	my ($input, $output) = @_;

	my @all_data;
	my $in = IO::File->new($input) || die "Can not open input file $input $!\n";
	while(<$in>) {
		chomp;
		my @a = split(/\t/, $_);
		push(@all_data, [@a[0,19,1,2,3,17,18,9,10,11,12,13,6,8]]);	# re-order the data, then put them to array 	
	}
	$in->close;
	
	@all_data = sort { ($a->[5] cmp $b->[5]) || ($a->[3] cmp $b->[3])} @all_data; # sort according to Genus and hit

	my $out = IO::File->new(">".$output) || die "Can not open output file $output $!\n";	
	print $out "Contig_ID\tContig_Seq\tContig_Len\tHit_ID\tHit_Len\tGenus\tDescription\tContig_start\tContig_end\tHit_start\tHit_end\tHsp_identity\tE_value\tHsp_strand\n";
	foreach my $data ( @all_data ) {
		print $out join("\t", @$data)."\n";
	}
	$out->close;
}
=head2
 arrange_col2 : 
=cut
sub arrange_col2
{
	my ($input, $output) = @_;

	my @all_data;
	my $in = IO::File->new($input) || die "Can not open input file $input $!\n";
	while(<$in>) {
		chomp;
		my @a = split(/\t/, $_);
		my $coverage= 1.0 * $a[7] / $a[1];				# get coverage
		push(@all_data, [@a[0,1,7],$coverage,@a[4,5,6,8,9]]);		# re-order the data, then put them to array
	}
	$in->close;
	
	@all_data = sort { ($a->[7] cmp $b->[7]) || ($b->[2] cmp $a->[2])} @all_data; # sort according to Genus and hit_covered(bp)

	my $out = IO::File->new(">".$output) || die "Can not open output file $output $!\n";	
	foreach my $data ( @all_data ) {
		print $out join("\t", @$data)."\n";
	}
	$out->close;	
}



