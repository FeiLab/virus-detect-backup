package Util;

use strict;
use warnings;

sub process_cmd {
	my ($cmd) = @_;	
	print "CMD: $cmd\n";
	my $ret = system($cmd);	
	if ($ret) {
		die "Error, cmd: $cmd died with ret $ret";
	}
	return($ret);
}

sub detect_FileType{
	my $file = shift;
	open(FH, $file) || die $!;
	my $line = <FH>;
	my $file_type;
	if	($line =~ m/^>/) { $file_type = 'fa'; }
	elsif	($line =~ m/^@/) { $file_type = 'fq'; }
	else	{ die "Error, can not detect the file type for file: $file\n"; }
	return $file_type;
}

1;
