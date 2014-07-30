#!/usr/bin/env perl
use strict;
use Getopt::Long qw(:config no_ignore_case bundling);
use Bio::SeqIO;
use Cwd;
my $usage = <<_EOUSAGE_;

#########################################################################################
# runVelvet.pl --velvet_dir [String] --parameters <File> --input_suffix <String>  --file_type <String> --output_suffix <String> 
# Required:
#  --parameters            A txt file containing a list of input file names with related parameters
#  --input_suffix          The suffix name of all the input files

# Options:
#  --velvet_dir           The directory to install the program velvet [workingdir/bin]
#  --file_type             The format of input files(fastq or fasta)  [fastq]
###########################################################################################


_EOUSAGE_
;

#################
# global vars   #
#################
our $velvet_dir;	# velvet的安装目录
our $parameters;	# 包括所有待处理样本和相应参数的文本文件名称（无后缀）
our $input_suffix = '';	# 输入文件的后缀名称
our $file_type="fastq";	# 输入文件的类型，当前只支持fastq和fasta格式
our $output_suffix; 	# 输出文件的后缀名称

################################
# set folder and file path     #
################################
our $WORKING_DIR=cwd();		# current folder
$velvet_dir=$WORKING_DIR."/bin";# velvet path 

##################
# get input para #
##################
GetOptions( 
	'velvet_dir=s'		=> \$velvet_dir,
	'parameters=s' 		=> \$parameters,
	'input_suffix=s' 	=> \$input_suffix,
	'file_type=s' 		=> \$file_type,
	'output_suffix=s' 	=> \$output_suffix
);
	     
die $usage unless $parameters;

# main
main: {
        my $i=0;
	my $resultDir;	
	open(IN, "$parameters");
	open(OUT1, ">velvet_assembly.result");

	while (<IN>) {
	
	$i=$i+1;
	chomp;	
	my @a = split "\t";		
	runVelvet($a[0],$a[1],$a[2]);				#每次循环读入一行，样本，Hash_len，Cov_cutoff。
	print "#processing sample $i: $a[0]\n";
	$resultDir=$a[0]."_".$a[1]."_".$a[2];
	my $contigName="./assembled/".$a[0].".".$output_suffix;	#输出文件名称
	system("mv ./$resultDir/contigs.fa $contigName");	#将结果文件转移到assembled文件夹
	system("rm ./$resultDir -r");
	my $count = `grep \'>\' $contigName `;#统计输出文件中的contig个数
	print OUT1 $a[0]."\t".$count;
	}
	close(IN);
	close(OUT1);
	system("rm velvet.log");
	print "###############################\n";
	print "All the samples have been processed by $0\n";
	#system("touch velvet_assembly.run.finished");
}

# subroutine
sub runVelvet {
	my $sample=shift;
	my $hash_length=shift;
	my $cov_cutoff=shift;
	my $outputDir=$sample."_".$hash_length."_".$cov_cutoff;
	my $file = $sample.$input_suffix;
	process_cmd($velvet_dir."/velveth $outputDir $hash_length -$file_type $file > velvet.log");
	process_cmd($velvet_dir."/velvetg $outputDir -cov_cutoff $cov_cutoff -min_contig_lgth 30 > velvet.log");	
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
