#! /usr/bin/perl

use warnings;
use strict;

	my ($in,$out) = @ARGV;
	my ($line,$id);

	#~ Passo. Abra os arquivo de entrada e saída
	open IN, "<$in" or die "Could not open input file $in: $!\n";
	open OUT, ">$out" or die "Could not open input file $in: $!\n";
	
	#~ Passo. Enquanto houver entrada
	while($line = <IN>){
		chomp($line);
		#~ Passo. Se a linha for do id e possuir espaço	
		if ($line =~ /^@(.+)\:(\d+)\:(.+)\:(\d+)\:(\d+)\:(\d+)\:(\d+)\s(.+)/){
			$line =~ s/$8//g;
		}
		print OUT "$line\n";
	}
	
	close IN;
	close OUT;	
		
exit 0;


