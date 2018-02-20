#! /usr/perl/bin

use warnings;
use strict;
#~ ------------------------------------------------------------------------
#~ Programa: frequency_counter.pl
#~ Data: 24/04/2015
#~ Este programa recebe como entrada uma lista de sequências com 
#~ no seguinte formato: 
#~ >seq1
#~ seq
#~ #idcdr
#~ cdr*
#~ <
#~ Em que (<) marca o fim do arquivo. Atenção, este programa lê as 
#~ strings considerando que estas estão em apenas uma linha, não estão 
#~ subdivididas em linhas de n caracteres. O programa armazena os dados
#~ em um vetor hashes. A chave do hash é a substring contendo cdrs, 
#~ e o campo freq armazena o número de sequências únicas igual a uma dada
#~ chave. O programa usa um outro hash para armazenar todas as sequências e
#~ seus respectivos ids que possuem a mesma subtring de cdrs. 
#~ Finalmente, é impressa uma lista ordenada por frequência, contendo 
#~ o id da cdr e sua frequência, a substring contendo cdrs, e todas as 
#~ sequências e seus ids que contem a substring.
#~ ------------------------------------------------------------------------

	#~ Passo. Declaração de variáveis
	my ($id,$seq,$idcdr,$cdr,$line,$temp);
	my ($flag,$len,$count,$c,$i,$k,$libsize,$freq);
	my ($input,$inputfiltered,$output) = @ARGV;
	my (@seqs);
	my %lista;	
	$cdr = "";
	$i = 0;
	$count = 0;
	$k = 0;
#~ -----------------------------------------------------------------------
								#~ Leitura
#~ -----------------------------------------------------------------------

	#~ Passo. Se o usuário não informou todos os argumentos
	if (@ARGV < 3){
		#~ Passo. Imprima mensagem de erro e pare a execução
		die "frequency_counter3.pl: perl frequency_counter3.pl <inputfile_protein.fasta> <input_filtered.fasta> <outputfile>\n";
	}
	#~ Passo. Abra o arquivo de entrada
	open IN, "<$input" or die "Could not find input file $input\n";

	#~ Passo. Enquanto houver entrada
	while($line = <IN>){
		#~ Passo. Remova o '\n' da linha
		chomp($line);
		#~ Passo. Se linha for do id
		if($line =~ /^>(.+)/ || $line =~ /</){
			#~ Passo. Se tiver terminado de ler os dados da sequência anterior
			if ($cdr =~ /(.+)\*$/){
				#~ Passo. Obtenha o tamanho da substring cdr sem "*"
				$len = length($cdr) - 1;
				#~ Passo. Copie $cdr sem "*" para $temp
				$temp = substr($cdr,0,$len);
				#~ Passo. Busque a cdr nos hashes do elemento atual
				$c = busca_cdr(\$temp,\$k,\@seqs);
				#~ Passo. Se a cdr foi encontrada
				if ($c != -1){
					#~ Passo. Incremente a frequência da cdr
					$seqs[$c]{$temp}{freq} = $seqs[$c]{$temp}{freq} + 1;
					#~ Passo. Armazene id e sequência que contem a substring de cdrs
					$lista{$temp}{id}[$seqs[$c]{$temp}{freq}-1] = $id;
					$lista{$temp}{seq}[$seqs[$c]{$temp}{freq}-1] = $seq;
				}
				else{
					#~ Passo. Inicialize a frequência da cdr
					$seqs[$i]{$temp}{freq} = 1;;
					#~ Passo. Armazene id e sequência que contem a substring de cdrs
					$lista{$temp}{id}[0] = $id;
					$lista{$temp}{seq}[0] = $seq;
					#~ Passo. Incremente o contador de sequências
					$count++;
				}				
				#~ Passo. Armazene o id da cdr em lista
				$lista{$temp}{idcdr} = $idcdr;
				#~ Passo. Se 10000 já foram armazenadas no hash do elemento atual
				if ($count % 10000 == 0){
					#~ Passo. Vá para o próximo elemento do vetor de hashes
					$i++;
					#~ Passo. Volte um elemento
					$k = $i - 1;
				}
				else{
					#~ Passo. Copie $i para $k
					$k = $i;
				}				
			}
			#~ Passo. Armazene linha na temporária id
			$id = $line;
			#~ Passo. Inicialize seq
			$seq = "";
			#~ Passo. Inicialize flag
			$flag = 1;
		}
		#~ Passo. Senão se a linha for do idcdr
		elsif($line =~ /^#(.+)/){
			#~ Passo. Armazene linha na temporária idcdr
			$idcdr = $1;
			#~ Passo. Inicialize cdr
			$cdr = "";
			#~ Passo. Inicialize flag
			$flag = 2;
		}
		#~ Passo. Senão se a linha for de seq
		elsif($flag == 1){
			#~ Passo. Armazene linha em seq
			$seq = $line;
		}
		else{
			#~ Passo. Armazene a linha em cdr
			$cdr = $line;
		}
	}
	
	#~ Passo. Feche o arquivo de entrada
	close(IN);
#~ ------------------------------------------------------------------------
					#~ Processamento
#~ ------------------------------------------------------------------------	
	
	#~ Passo. Obtenha o tamanho da biblioteca de sequências filtradas
	$libsize = `grep -cP "^>" $inputfiltered`;
	chomp($libsize);
	$libsize = $libsize + 0;
	print "libsize : $libsize\n";
	#~ Passo. Se o último elemento do vetor totaliza 10000 hashes
	if ($count % 10000 == 0){
		#~ Passo. Volte um elemento
		$i = $i -1;
	}
	
	#~ Passo. Inicialize $c
	$c = 0;
	#~ Enquanto existirem elementos no vetor de hash
	while($c <= $i){
		#~ Passo. Para cada chave cdr
		foreach $line (keys %{$seqs[$c]}){
			#~ Passo. Copie a frequência da cdr para o hash lista
			$lista{$line}{freq} = $seqs[$c]{$line}{freq};
		}
		#~ Passo. Vá para o próximo elemento do vetor
		$c++;
	}
#~ ------------------------------------------------------------------------
						#~ Escrita
#~ ------------------------------------------------------------------------	
	#~ Passo. Abra o arquivo de saída
	open OUT, ">$output" or die "frequency_counter2.pl: Could not create output file $output\n";
	#~ Passo. Para cada cdr
	foreach $k (sort{$lista{$b}{freq} <=> $lista{$a}{freq}} keys %lista){
		#~ Passo. Inicialize o índice dos vetores de ids e seqs
		$c = 0;
		#~ Passo. Calcule a frequência relativa da cdr
		$freq = ($lista{$k}{freq} / $libsize) * 100000000000000;
		#~ Passo. Imprima o id da cdr e frequencia
		print OUT "#$lista{$k}{idcdr}|$libsize|$freq\n";
		#~ Passo. Imprima cdr 
		print OUT "*$k\n";
		#~ Passo. Enquanto existirem ids e sequências associadas a cdr atual
		while($c < $lista{$k}{freq}){
			#~ Passo. Imprima o id da sequência
			print OUT "$lista{$k}{id}[$c]\n";
			#~ Passo. Imprima a sequência
			print OUT "$lista{$k}{seq}[$c]\n";
			#~ Passo. Vá para o proximo id e próxima sequência
			$c++;
		}
	}
	#~ Passo. Feche o arquivo de saída
	close(OUT);
	
exit(0);

#~ ------------------------------------------------------------------------
				#~ Definição de subrotinas
#~ ------------------------------------------------------------------------

sub busca_cdr{
	
	#~ Passo. Declaração e atribuição de variáveis

	my $cdr = shift or die "frequency_counter2.pl: Could not receive string reference\n";
	$cdr = ${$cdr};
	my $i = shift or die "frequency_counter2.pl: Could not receive integer index reference\n";
	$i = ${$i};
	my $temp = shift or die "frequency_counter2.pl: Could not receive array\n";
	my @seqs = @{$temp};
	my $flag = -1;
	my $c = 0;
	#~ Passo. Enquanto não chegar ao fim do vetor e não achar a cdr
	while($c <= $i && $flag == -1)
	{
		#~ Passo. Se a cdr for chave de algum hash do elemento atual
		if(exists $seqs[$c]{$cdr})
		{
			#~ Passo. Modifique flag
			$flag = $c;
		}
		else
		{
			#~ Passo. Vá para o próximo elemento do vetor de hashes
			$c++;
		}
	}
	
	#~ Passo. Retorne o índice $c
	return $flag;
}
		
