#! usr/bin/perl

use warnings;
use strict;

#########################################################################
# Programa: numberab.pl													#	
# Data: 24/09/2014						     							#
# Autor: Heidi Muniz Silva					     						#
# Este programa recebe um arquivo contendo sequências proteicas,		# 
# executa o comando wget que envia o arquivo de entrada para o programa #
# Abnum. O programa Abnum numera e identifica as sequências de 			#	
# acordo com esquemas de numeração e identificação de imunoglobulinas. 	#
# O conteúdo de cada arquivo gerado pelo programa Abnum é redirecionado #
# para um arquivo de saída, de modo que contenha todas as sequências	#
# reconhecidas e numeradas, com seus respectivos identificadores.		#
#########################################################################




if (@ARGV < 2)
{
	die "numberab.pl: perl numberab.pl <inputfile> <outputfile>\n";
}
my ($input,$output1) = @ARGV;

# Passo. Abra o arquivo contendo as sequências

open IN, "<$input" or die "numberab.pl: Could not open input file $input: $! \n";

# Passo. Leia as sequências

my ($line,$header,$n,%seqs);
$n = 0;

while ($line = <IN>)
{
	chomp($line);
	if ($line =~ /^>(.+)\|P1(.+)/)
	{
		$header = $1;
		$seqs{$header}{number} = $n;
		$n++;
	}
	else
	{
		$seqs{$header}{seq} = $line;
	}
}

close IN;

# Passo. Numere as sequências e identifique as sequências

foreach $header (sort {$seqs{$a}{number} <=> $seqs{$b}{number}} (keys %seqs))
{
	# Passo. Abra o arquivo de saída, em que serão armazenadas as sequência numeradas

	open OUT, ">>$output1" or die "numberab.pl: Could not open output file $output1: $! \n";
	
	# Passo. Envie a sequência para o programa Abnum, que executa numeração
	
	system("wget \'http://www.bioinf.org.uk/cgi-bin/abnum/abnum.pl?plain=1&aaseq=$seqs{$header}{seq}&scheme=-k\'");
	 
	my $out;
	
	# Passo. Obtenha o tamanho do arquivo em bytes
	
	$out =`du -b abnum* \| cut -f1`;
	chomp($out);
	# Passo. Caso o arquivo possua mais de 1 byte de tamanho,
	# redirecione seu conteúdo para o arquivo de saída $output1
	
	if ($out > 1)
	{
		# Passo. Imprima o identificador da sequência atual
		
		print OUT ">$header\n";
		
		# Passo. Redirecione o conteúdo do arquivo gerado pelo programa Abnum para o arquivo de saída $output1
		
		system("cat abnum* >> $output1");
		
		close OUT;
		
	}
	
	
	# Passo. Remova os arquivos gerados pelos programa Abnum
	system("rm abnum*");
}
