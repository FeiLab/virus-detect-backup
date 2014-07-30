#!/usr/bin/perl -w 
use strict; 

if (@ARGV < 1)
{
  print "usage: errorcheck.pl inputfile\n";
  exit(0);
}

our $input = $ARGV[0]; 
our @errors=("Error","Aborted","cannot","No","core dumped"); 

for my $err (@errors) {
	print "# This is information of ".$err."\n";
	my $return = `grep -n \'$err\' $input`;
	print $return."\n";
}



