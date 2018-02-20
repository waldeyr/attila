#! /usr/bin/perl

#~ Programa: rscript_creator.pl
#~ Data: 03/02/2016
#~ Este programa recebe como entrada os seguintes argumentos:
#~ $i = path do arquivo da biblioteca original, antes da filtragem, em formato fasta
#~ $f = path do arquivo da biblioteca final, antes da filtragem, em formato fasta
#~ $csv = path do arquivo csv contendo o número de reads por etapa
#~ $scriptlen = path do script r para criar gráfico de proporção de reads com tamanho adequado
#~ $scriptask = path do script r para criar gráfico de número de reads por etapa
#~ $dir = diretório onde serão armazenados os pngs dos gráficos
#~ $readlen = tamanho mínimo de read
#~ O programa gera dois scripts R. Um deles cria um gráfico de proporção de reads com tamanho adequado
#~ e o outro, um gráfico de número de reads por etapa. 

	use warnings;
	use strict;
#~ --------------------------------------------------------------------------------------------------------------------------------------------------------------	


	my ($i, $f, $csv, $scriptlen, $scriptask, $dir, $readlen, $libtype) = @ARGV;
	my ($line, $len, $seq, $n, $input, $total, $good, $pgood, $pbad, $lib, $aux);
	
	if (@ARGV < 8){
		die "rscript_creator.pl: perl script_creator.pl <input_initial.fasta> <input_final.fasta> <input.csv> <script1.r> <script2.r> <path> <integer> <string>"
	}
#~ --------------------------------------------------------------------------------------------------------------------------------------------------------------
	#~ Crie o script R 1
#~ ---------------------------------------------------------------------------------------------------------------------------------------------------------------
	#~ Passo. Escreva o header do script R 1
	open(OUT, '>:encoding(UTF-8)', $scriptlen) or die "Could not create script1.r: $!\n";
	print OUT "#! /usr/bin/env  Rscript\n";
	print OUT "library(\"ggplot2\")\n";
	print OUT "library(\"scales\")\n";
	print OUT "input <- read.table(header= TRUE, stringsAsFactors = FALSE, text = '\n";
	print OUT "Library\t\t\t\t\tLength\t\t\t\t\t\treads\n";
	
	$n = 0;
	while($n <= 1){
		#~ Passo. Inicialize os contadores
		$good = 0;
		$pgood = 0;
		$pbad = 0;
		#~ Passo. Selecione o arquivo da biblioteca original
		if ($n == 0){
			$input = $i;
			$lib = "Initial";
		}
		#~ Passo. Selecione o arquivo da biblioteca final
		else{	
			$input = $f;
			$lib = "Final";
		}
		#~ Passo. Calcule o número de reads da biblioteca 
		$aux = `grep -cP "^>" $input`;
		$total = $aux;
		#~ Passo. Abra o arquivo 
		open IN, "<$input" or die "Could not open input file $input: $!\n";
		$seq = "";
		#~ Passo. Enquanto houver entrada
		while ($aux > 0){
			$line = <IN>;
			if (defined($line)){
				chomp($line);
			}
			#~ Passo. Se já tiver terminado de o arquivo ou de ler uma sequência
			if (!defined($line) || $line =~ /^>(.+)/){
				#~ Passo. Calcule o tamanho da sequência atual
				$len = length($seq);
				if ($len > 1){
					#~ Passo. Se a sequência tiver tamanho adequado
					if ($len >= $readlen){
						#~ Passo. Incremente o contador de sequências de tamanho adequado
						$good++;
					}
					$aux--;
				}					
				$seq = "";
			}
			else{
				$seq .= $line;
			}
		}
		close IN;
		#~ Passo. Calcule porcentagens de good e bad
		$pgood = $good / $total;
		$pgood = sprintf("%.4f", $pgood);
		$pbad = 1 - $pgood;
		$pbad = sprintf("%.4f", $pbad);
		
		#~ Passo. Escreva os dados no script R 1
		print OUT "$lib\t\t\t\t\tInadequate\t\t\t\t\t$pbad\n";
		print OUT "$lib\t\t\t\t\tAdequate\t\t\t\t\t$pgood\n";
		$n++;
	}
	
	$aux = "";
	$aux = $dir . "proportion_read_length_" . $libtype . ".png";
	#~ Passo. Escreva as configurações do gráfico no script R 1
	print OUT "')\n";
	print OUT "png(file=\"$aux\", res=300, width=7, height=7, units = 'in')\n";
	print OUT "m <- ggplot(data = input,aes(x= factor(Library, levels = unique(Library)),y=reads,fill=Length)) + geom_bar(width=0.3,stat=\"identity\") + scale_y_continuous(labels=percent)\n";
	print OUT "m <- m +  scale_fill_discrete(name=\"Read Length\")\n";
	print OUT "m <- m + xlab(\"Library\") + ylab(\"Reads\")\n";
	print OUT "m <- m + ggtitle(\"Proportion of Reads with Adequate Length\")\n";
	print OUT "print(m)\n";
	print OUT "garbage <- dev.off()\n";
	
	close OUT;
	
#~ --------------------------------------------------------------------------------------------------------------------------------------------------------------
	#~ Crie o script R 2
#~ ---------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	#~ Passo. Abra o arquivo csv
	open IN, "<$csv" or die "Could not open input file $input: $!\n";
	
	#~ Passo. Abra o script R 2 e escreva o header
	open(OUT, '>:encoding(UTF-8)', $scriptask) or die "Could not create script2.r: $!\n";
	
	print OUT "#! /usr/bin/env Rscript\n";
	print OUT "library(\"ggplot2\")\n";
	print OUT "library(\"scales\")\n";
	print OUT "input <- read.table(header= TRUE, stringsAsFactors = FALSE, text = '\n";
	print OUT "Library\t\t\t\t\tProcess\t\t\t\t\t\treads\n";
	
	while ($line = <IN>){
			
		chomp($line);
		if($line =~ /^R0,(\d+),raw$/){
			print OUT "Initial\t\t\t\t\tNone\t\t\t\t\t$1\n";
		}
		if($line =~ /^RN,(\d+),raw$/){
			print OUT "Final\t\t\t\t\tNone\t\t\t\t\t$1\n";
		}
		if($line =~ /^R0,(\d+),joining$/){
			print OUT "Initial\t\t\t\t\tJoining\t\t\t\t\t$1\n";
		}	
		if($line =~ /^RN,(\d+),joining$/){
			print OUT "Final\t\t\t\t\tJoining\t\t\t\t\t$1\n";
		}
		if($line =~ /^R0,(\d+),filtering$/){
			print OUT "Initial\t\t\t\t\tFiltering\t\t\t\t\t$1\n";
		}
		if($line =~ /^RN,(\d+),filtering$/){
			print OUT "Final\t\t\t\t\tFiltering\t\t\t\t\t$1\n";
		}
		if($line =~ /^R0,(\d+),translation$/){
			print OUT "Initial\t\t\t\t\tTranslation\t\t\t\t\t$1\n";
		}
		if($line =~ /^RN,(\d+),translation$/){
			print OUT "Final\t\t\t\t\tTranslation\t\t\t\t\t$1\n";
		}
		if($line =~ /^R0,(\d+),frequency$/){
			print OUT "Initial\t\t\t\t\tFrequency_Calculation\t\t\t\t\t$1\n";
		}
		if($line =~ /^RN,(\d+),frequency$/){
			print OUT "Final\t\t\t\t\tFrequency_Calculation\t\t\t\t\t$1\n";
		}
		if($line =~ /^Selected,(\d+),enrichment$/){
			print OUT "Final\t\t\t\t\tEnrichment\t\t\t\t\t$1\n";
		}						
	}
	
	close IN;
	
	$aux = "";
	$aux = $dir . "reads_by_task_" . $libtype . ".png";
	print OUT "')\n";
	print OUT "png(file=\"$aux\", res=300, width=9, height=7, units = 'in')\n";
	print OUT "m <- ggplot(data = input,aes(x= factor(Process, levels = unique(Process)),y=reads,fill=Library)) + geom_bar(stat=\"identity\",position=position_dodge())\n";
	print OUT "m <- m + xlab(\"Task\") + ylab(\"Reads\")\n";
	print OUT "m <- m + ggtitle(\"Total Number of Reads by Task\")\n";
	print OUT "options(scipen=999)\n";
	print OUT "print(m)\n";
	print OUT "garbage <- dev.off()\n";
	
	close OUT;
	
	exit(0);
