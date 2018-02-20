#!/usr/bin/perl

	use warnings;
	use warnings;
#~ ---------------------------------------------------------------------	
	#~ Programa: statscript_creator.pl
	#~ Data: 23/06/2016
	#~ Este programa recebe como entrada um arquivo contendo as n primeiras
	#~ sequências enriquecidas, um arquivo contendo as sequências numeradas 
	#~ pelo Abnum, e o tamanho das bibliotecas original e final. O programa
	#~ obtém o id e as frequências inicial e final de cada clone, e então
	#~ cria um script R que executa o teste de diferença de proporção e o cál
	#~ culo de intervalo de confiança para cada sequência. O script R cria 
	#~ um arquivo de saída onde são armazenados os resultados dos testes
	#~ estatísticos.
#~ ---------------------------------------------------------------------

	#~ Declaração de variáveis
	my ($enriched,$abnum,$path,$scriptr,$sizer0,$sizern,$libtype,$ncandidates) = @ARGV;
	my ($line,%seqs,$count);
	
	#~ Passo. Abra o arquivo de sequências numeradas pelo Abnum
	open IN, "<$abnum" or die "Could not open input file $abnum: $!\n";
	
	$count = 0;
	#~ Passo. Leia o arquivo e armazene os id das sequências
	while($line=<IN>){
		
		chomp($line);
		
		if ($line =~ /^>(.+)\|FRAME(.+)/){
				$seqs{$1}{freq1} = 0;
				$seqs{$1}{freq2} = 0;
				$seqs{$1}{count} = $count;
				$count++;
		}
	}

	close IN;
	
	
	#~ Passo. Abra o arquivo de sequências enriquecidas
	open IN, "<$enriched" or die "Could not open input file $enriched: $!\n";
	
	#~ Passo. Leia o arquivo e armazene os id das sequências
	while($line=<IN>){
		
		chomp($line);
		
		if ($line =~ /^>(.+)\|FRAME:\d\|FOLD-CHANGE:(\d+\.\d*)\|P1:(\d+\.\d*)\|P2:(\d+\.\d*)/){
			
				
				if (exists($seqs{$1})){
					$seqs{$1}{freq1} = $3 / 100000000000000 ;
					$seqs{$1}{freq2} = $4 / 100000000000000 ;
			}
		}
	}

	close IN;
	
	
	#~ Passo. Abra o arquivo de saída
	open OUT, ">$scriptr" or die "Could not create R script $scriptr: $!\n";
	
	#~ Passo. Imprima no script R
	print OUT "#! /usr/bin/env Rscript\n";
	if ($libtype == 0){
		print OUT "sink(\"$path/vhoutputRstats.txt\")\n";
	}
	else{
		print OUT "sink(\"$path/vloutputRstats.txt\")\n";
	}
	
	#~ Passo. Imprima o cálculo do alfa de Bonferroni
	if ($ncandidates > 0){
		print OUT "alfaBonf <- 0.05 / $ncandidates\n";
		print OUT "cat('alfaBonf:\\t', alfaBonf)\n";
		print OUT "cat('\\n')\n";
	
	
		#~ Passo. Para cada sequência, calcule o p-valor e o intervalo de confiança para 95%
		foreach $line (sort{$seqs{$a}{count} <=> $seqs{$b}{count}} keys %seqs){
				print OUT "cat('>$line\\n')\n";
				print OUT "pantes <- $seqs{$line}{freq1}\n";
				print OUT "pdepois <- $seqs{$line}{freq2}\n";
				print OUT "pantes_menos_pdepois <- pdepois - pantes\n";
				print OUT "n1 <- $sizer0\n";
				print OUT "n2 <- $sizern\n";
				print OUT "pantesdepois <- (n1 * pantes + n2 * pdepois) / (n1 + n2)\n";
				print OUT "erropadiff <- sqrt(pantesdepois * (1 - pantesdepois) * (1/n1 + 1/n2))\n";
				print OUT "Zobs <- (pantes_menos_pdepois - 0) / erropadiff\n";
				print OUT "pvalor <- 1 -pnorm(Zobs)\n";
				print OUT "cat('p-value:\\t', pvalor)\n";
				print OUT "cat('\\n')\n";
				
				print OUT "EPdiffIC <- sqrt(pantes*(1-pantes)/n1+pdepois*(1-pdepois)/n2)\n";
				print OUT "infIC <- pantes-pdepois-1.96*EPdiffIC\n";
				print OUT "cat('infIC:\\t', infIC)\n";
				print OUT "cat('\\n')\n";
				print OUT "supIC <- pantes-pdepois+1.96*EPdiffIC\n";
				print OUT "cat('supIC:\\t', supIC)\n";
				print OUT "cat('\\n')\n";
		}
	}
	else{
		print OUT "cat('No candidates\\n')\n";
	}
	close OUT;
	
	
	exit 0;
