#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
my $usage = <<_EOUSAGE_;

#########################################################################################
# extractFromFasta.pl --inputfile <FILE> --type <String> --query <String> --output1 <FILE> --output2 <FILE>
# Required:
#  --inputfile        a fasta file containing all sequences
#  --type             "single" or "list"
#  --query            an sequence ID or a list of sequence IDs
#  --output1          a file to contain all the sequences in query
#  --output2          a file to contain all the remaiding sequences
###########################################################################################

_EOUSAGE_

	;


my $fasta_file;#序列所在的文件
my $type;   #"single"表示用序列ID（字符串）查找,"list"表示用一个序列ID列表查找
my $query;  #对应"single"的序列ID，或者对应"list"的列表文件
my $output1;#query包括的输出到这个文件
my $output2;#剩下的的输出到这个文件

&GetOptions( 'inputfile=s' => \$fasta_file,#
             'type=s' => \$type,
             'query=s' => \$query,
	     'output1=s' => \$output1,
	     'output2=s' => \$output2,
			 );

unless ($fasta_file&&$type&&$query) {
	die $usage;
}


#提取需要保留的序列的name到%list变量中
my %list=();
if ($type eq "list") {#如果需要提取的fasta序列的name是一个列表
	open(IN, $query);#就要读这个列表所在的文件
	while(<IN>) {
		chomp($_);
		#$_=~s/>//; #序列名称如果以>开头，替换为空
		$list{$_} = 1;#每次存入一个name进入变量$list
	}
	close(IN);
}
elsif($type eq "single") {#如果需要提取的fasta序列的只有一个name
	#$query=~s/>//; #序列名称如果以>开头，替换为空
	$list{$query}=1;#存入一个name进入变量$list
}
else {
	die &usage();
}

#读fasta file文件，
open(IN, $fasta_file);
open(OUT1, ">$output1");
open(OUT2, ">$output2");
my $flag = "off";
while(<IN>) {
	if($_ =~ /^>(\S+)/) {
		my $head = $1;
		if(defined $list{$head}) {#如果包括这个name
			print OUT1 $_;#就输出到OUT1
			$flag = "on";#同时改变标志，表示后面的序列需要继续向OUT1输出
		}
		else {#如果不包括这个name
			if($type eq "single" && $flag eq "on") {#表示name找到过一次，就不再继续找了
				exit;
			}
			print OUT2 $_;#就输出到OUT2
			$flag = "off";#同时改变标志，表示后面的序列需要继续向OUT2输出
		}
	}
	else {
		if($flag eq "on") {#表示为"on"
			print OUT1 $_;#后面的序列需要继续向OUT1输出
		}
		else {
			print OUT2 $_;#否则，后面的序列需要继续向OUT2输出
		}
	}
}
close(IN);
close(OUT1);
close(OUT2);
