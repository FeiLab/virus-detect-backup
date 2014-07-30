package Util;

use strict;
use warnings;

sub process_cmd {
	my ($cmd, $debug) = @_;	
	if (defined $debug && $debug > 0) { print "[CMD FOR DEBUG]: $cmd\n"; }
	my $ret = system($cmd);	
	if ($ret) {
		die "[ERR IN CMD]: \n$cmd\n Died with ret: $ret\n";
	}
	return($ret);
}

sub detect_FileType{
	my $file = shift;
	open(FH, $file) || die $!;
	my $line = <FH>;
	my $file_type;
	if	($line =~ m/^>/) { $file_type = 'fasta'; }
	elsif	($line =~ m/^@/) { $file_type = 'fastq'; }
	else	{ die "Error, can not detect the file type for file: $file\n"; }
	return $file_type;
}

sub print_user_message {
	my @message = @_;
	#print "\n";
	my $time = get_time();
	foreach my $line (@message) { 
		print $time." ".$line."\n";
	}
	#print "\n";
}

sub print_user_submessage{
	my @message = @_;
	foreach my $line (@message) { 
		print "   ".$line."\n";
	}
}

sub get_time {
	my $time = `date +"%D %T"`;
	chomp $time;
	$time = "[$time]";
	return $time
}

1;
