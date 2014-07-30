#!/usr/bin/perl
#自己编模块解析blastn结果（速度快），从paired输出格式到table格式
#bioperl的定义：$hsp->start('hit')总是不大于$hsp->end('hit')，这与blast输出的table格式不同
#这个程序保持了bioperl的格式，必须加一列strand信息

if (@ARGV < 2)
{
  print "usage: blast_parse_table22.pl input output\n";
  exit(0);
}

our $input = $ARGV[0];
our $output = $ARGV[1];

my ($query_name, $query_length, $hit_name, $hit_length, $hsp_length, $identity, $evalue, $score, $strand, $query_start, $query_end, $hit_start, $hit_end, $query_to_end, $hit_to_end);
my ($identity2, $aligned_query,$aligned_hit,$aligned_string); #增加4列信息
my (%hsp, $hsp_count); 
%hsp = ();$hsp_count=0; 
my $query_sequenc; 
my $subject_sequenc; 	
my $is_hsp = 1;

open(IN, "$input");
open(OUT, ">$output");
#print OUT "query_name\tquery_length\thit_name\thit_length\thsp_length\tidentity\tevalue\tscore\tstrand\tquery_start\tquery_end\thit_start\thit_end\tidentity2\taligned_query\taligned_hit\taligned_string;
while (<IN>) 
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
					$query_start."\t".$query_end."\t".$hit_start."\t".$hit_end."\t".$identity2."\t".$aligned_query."\t".$aligned_hit."\t".$aligned_string;
			$hsp{$hsp_count} = $hsp_info;
		}

		#########################################
		#  输出属于上一对query-hit的所有hsp		#
		#########################################
		if (scalar(keys(%hsp)) > 0)
		{
		    my $aa = scalar(keys(%hsp));
			foreach my $hsp_num (sort {$a<=>$b} keys %hsp)
			{
				print OUT $hsp{$hsp_num}."\n";
			}#循环结束
		}
		#####################################
		#  开始记录一个新的query			#
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
					$query_start."\t".$query_end."\t".$hit_start."\t".$hit_end."\t".$identity2."\t".$aligned_query."\t".$aligned_hit."\t".$aligned_string;
			$hsp{$hsp_count} = $hsp_info;
			$is_hsp = 0;			
		}
		#################################
		#  开始记录一个新的hit			#
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
			                $query_start."\t".$query_end."\t".$hit_start."\t".$hit_end."\t".$identity2."\t".$aligned_query."\t".$aligned_hit."\t".$aligned_string;
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
		$aligned_query = ""; $aligned_hit= ""; $aligned_string= "";

	}
	
	elsif (/\s+Identities = (\d+)\/(\d+)\s+\((\S+)\)/ && $hsp_count >= 1)
	{
		$identity = $1/$2*100;
		$identity = sprintf("%.".(2)."f", $identity);
		if ( $1 > $2 ) { $hsp_length = $1; $hsp_length =~ s/,//g; } else { $hsp_length = $2; $hsp_length =~ s/,//g; }
		$identity2 = $1."/".$2."(".$3.")";
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
		chomp($_);
		$aligned_query .= $_.",";
		my $next_line=<IN>;
		chomp($next_line); 
		$aligned_string .= $next_line.",";
	}

	elsif (/Sbjct\:\s(\d+)\s+(\S+)\s(\d+)/ && $hsp_count >= 1)
	{
		if ( $strand == -1 )
		{
			if ($hit_end == 0) { $hit_end = $1; $hit_end =~ s/,//g; };
			$hit_start = $3;
			$hit_start =~ s/,//g;
		}
		else
		{
			if ($hit_start == 0) { $hit_start = $1; $hit_start =~ s/,//g;};
			$hit_end = $3;
			$hit_end =~ s/,//g;
		}
		chomp($_);
		$aligned_hit .= $_.",";
	}
	else
	{
		next;
	}
}
close(IN);
close(OUT);
