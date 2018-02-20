#! usr/bin/perl

use warnings;
use strict;

#~ ----------------------------------------------------------------------
# Programa: get_nsequences.pl
#~ Este programa copia as primeiras n sequências do arquivo de entrada, e
#~ as escreve num arquivo de saída.
#~ ----------------------------------------------------------------------

	#~ Passo. Declaração de variáveis
	my ($number,$line,$header,%seq,$n,$seq);
	my ($input,$output,$total);
	$number = 0;
	$header = "";
	$seq = "";

	# Passo. Leitura

	if (@ARGV < 3){	
		die "get_nsequences: perl get_nsequences.pl <inputfile.fasta> <outputfile> <integer>\n";
	}

	($input,$output,$n) = @ARGV;
		
	# Passo. Abra o arquivo de entrada
	open IN, "<$input" or die "get_nsequences: Could not open input file $input: $! \n";

	#~ Passo. Conte o número de sequências do arquivo
	$total = `grep -c -P ^'>' $input`;
	
	#~ Passo. Se o número de sequências solicitado for maior que o total de sequências
	if($n > $total){
		#~ Passo. Atribua o conteúdo de $total a $n
		$n = $total;
	}
	#~ Passo. Enquanto não ler n sequências
	while ($n > 0){
		#~ Passo. Leia a linha
		$line = <IN>;
		#~ Passo. Se tiver terminado de ler os dados da sequência anterior
		if(length($seq) > 1){
			# Passo. Se linha for do id ou tiver terminado de ler o arquivo
			if((!defined($line)) || ($line =~ /^>(.+)/)){
				#~ Passo. Armazene os dados da sequência anterior no hash
				$seq{$header}{read} = $seq;
				$seq{$header}{number} = $number;
				$number++;
				$n--;
			}
		}
		#~ Passo. Se a linha tiver algum conteúdo definido
		if(defined($line)){
			#~ Passo. Remova o '\n' da linha
			chomp($line);
			#~ Passo. Se a linha for do id
			if($line =~ /^>(.+)/){
			#~ Passo. Armazene a linha em $header
			$header = $1;
			#~ Passo. Inicialize $seq
			$seq = "";			
			}
			else{
				# Passo. Concatene a linha em $seq
				$seq .= $line;
			}
		}
	}	
	close IN;

	open OUT, ">$output" or die "get_nsequences: Could not open output file $output: $!\n";

	
	foreach $line (sort{$seq{$a}{number} <=> $seq{$b}{number}} keys(%seq)){
		#~ Passo. Imprima o id e a sequência no arquivo de saída
		print OUT ">$line\n$seq{$line}{read}\n";
	}
	
	close OUT;
	
	exit(0);
