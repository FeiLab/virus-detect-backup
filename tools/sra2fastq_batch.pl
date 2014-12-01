#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;
#��ͨ�����ѭ��ȥ���࣬Ȼ������base correction
my $usage = <<_EOUSAGE_;

#########################################################################################################
# sra2fastq_batch.pl --file_list <FILE>
#
# Required(1):
#  --file_list       A txt file containing a list of input file names without any suffix
#########################################################################################################

_EOUSAGE_
	;
	
#################
##   ȫ�ֱ���  ##
#################
our $file_list; #�������д�������������ı��ļ����ƣ��޺�׺��

#################
## �����������##
#################
&GetOptions( 'file_list=s' => \$file_list		 			 
			 );
			 
unless ($file_list) {#�������ͨ������õ�
die $usage;}

################################
##   ��������Ŀ¼���ļ���·�� ##
################################
our $WORKING_DIR=cwd();#����Ŀ¼���ǵ�ǰĿ¼
our $BIN_DIR=$WORKING_DIR."/tools";#���п�ִ���ļ����ڵ�Ŀ¼

#################
##  ������ʼ ##
#################
open(IN1,$file_list) || die "Can't open the file $file_list\n";
my $sample;
my $j=0;
while(<IN1>){
	$j=$j+1;
	chomp;
	$sample=$_; #ÿ��ѭ������һ�У��������붼�Ǵ���������ļ��������޺�׺����
	process_cmd ("$BIN_DIR/fastq-dump $_.sra");		
}
close(IN1);
print "###############################\n";
print "All the samples have been processed by $0\n";

#################
##    �ӳ���   ##
#################
sub process_cmd {
	my ($cmd) = @_;	
	print "CMD: $cmd\n";
	my $ret = system($cmd);	#�ɹ��ͷ���0�������ʧ��
	if ($ret) {
		die "Error, cmd: $cmd died with ret $ret";
	}
	return($ret);
}