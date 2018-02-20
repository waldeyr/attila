#! usr/bin/perl

	use warnings;
	use strict;
	

#~ --------------------------------------------------------------------------------------
#~ Programa: get_ntsequences.pl
#~ Este programa recebe como entrada 3 arquivos
#~ arquivos de entrada:
#~ arquivo 1: lista de sequências de aminoácidos numeradas pelo Abnum em formato fasta
#~ arquivo 2: lista de sequências de aminoácidos enriquecidas, enviada ao Abnum
#~ arquivo 3: arquivo de sequências de nucleotídeos do ciclo 3, cuja tradução originou
#~ as sequências do arquivo 2.
#~ O programa obtem as sequências de nucleotídeos correspondente à sequências
#~ de aminoácidos do arquivo 1. Para isso, os ids do arquivo 1 e a frame são 
#~ armazenados em um hash. Em seguida são armazenadas as sequências de aminoácidos do 
#~ arquivo 2 e de nucleotídeos do arquivo 3, cujo id é igual a um id do hash.
#~ Toda sequência do arquivo 1 é uma substring de uma sequência do arquivo 2.
#~ O programa encontra a posição inicial da subtring do arquivo 1, na sequência
#~ do arquivo 2; calcula o número de nucleotídeos necessários para codificar
#~ a substring, e extrai da sequência de nucleotídeos uma substring que começa
#~ na posição correspondente ao primeiro aminoácido da sequência do arquivo 1 e cujo
#~ tamanho é o triplo da sequência do arquivo 1. Finalmente, as substrings 
#~ de nucleotídeos são impressas em linhas de até 70 caracteres no arquivo de
#~ saída.
#~ ---------------------------------------------------------------------------------------
	
#~ ---------------------------------------------------------------------------------------
#~ 								Main
#~ ---------------------------------------------------------------------------------------


	#~ Passo. Declaração de variáveis
	my %seqs;
	my ($line,$line2,$count,$len,$i,$aux);
	my ($start,$file,$id,$lensub,$temp);
	my $input1 = $ARGV[0];
	my $output = $ARGV[3];
	$count = 0;
			
	#~ Passo. Se o usuário não informou o número correto de argumentos
	if(@ARGV < 4){
		#~ Passo. Imprima mensagem de erro e pare a execução
		die "get_ntsequences: perl get_ntsequences.pl <inputfile1> <inputfile1> <inputfile> <outputfile>\n";
	}

#~ ---------------------------------------------------------------------------------------
#~ 						Leitura do arquivo 1 
#~ ---------------------------------------------------------------------------------------
	#~ Passo. Abra o arquivo de entrada 1
	open IN, "<$input1" or die "get_ntsequences: Could not open input file $input1: $!\n";
	#~ Passo. Enquanto houver entrada
	while($line = <IN>){
		#~ Passo. Remova o '\n'
		chomp($line);
		#~ Passo. Se a linha for do id
		if($line =~ /^>(.+)\|FRAME:(\d)\|FOLD-CHANGE:((\d+\.\d+)?)$/){
			#~ Passo. Inicialize campo seq usando o id como chave e incremente o contador de leitura
			if(defined($1)){
				$seqs{$1}{id} = $1 . "|FOLD-CHANGE:" . $3;
				$seqs{$1}{subcdr} = "";
				$seqs{$1}{protein} = "";
				$seqs{$1}{dna} = "";
				$seqs{$1}{sense} = "";
				$seqs{$1}{frame} = $2;
				$seqs{$1}{cds} = "";
				$seqs{$1}{number} = $count;
				$count++;
			}
		}
		else{
			#~ Passo. Armazene a sequência numerada pelo Abnum
			$seqs{$1}{subcdr} = $line;
		}			
	}
	#~ Para cada sequência
	foreach $line(keys %seqs){
		#~ Passo. Remova o traço de gap
		$seqs{$line}{subcdr} =~ s/-//g;
	}
	#~ Passo. Feche a entrada
	close IN;
#~ ---------------------------------------------------------------------------------------
#~ 					Leitura dos arquivos 2 e 3
#~ ---------------------------------------------------------------------------------------
	$file = 1;
	#~ Enquanto existirem arquivos
	while($file <= 2){
		#~ Passo. Abra o arquivo de entrada 2
		open IN, "<$ARGV[$file]" or die "get_ntsequences: Could not open input file $ARGV[$file]: $!\n";
		$aux = $count;
		#~ Passo. Enquanto houverem sequências 
		while($aux > 0 && defined($line = <IN>)){
			#~ Passo. Remova o '\n'
			chomp($line);
			#~ Passo. Se a linha for do id
			if($line =~ /^>(.+)\|FRAME(.+)/){
				#~ Passo. Armazene $1 em id
				$id = $1;
				$temp = ">" . $1 . "|FRAME" . "$2";
			}
			else{
				#~ Passo. Se o id existir no hash
				if(exists $seqs{$id}){
					#~ Passo. Se a linha for de sequência de proteína
					if($file == 1){
						#~ Passo. Armazene a linha no campo proteina do hash
						$seqs{$id}{protein} = $line;
					}
					#~ Passo. Senão, se a linha for de sequência de dna
					else{
						#~ Passo. Identifique qual sentido da fita codificante
						if($temp =~ /^>(.+)\|FRAME:(\d)(.?)/){
							$i = $3;
							if($i =~ /\+/){
								#~ Passo. Fita plus
								$seqs{$id}{sense} = '+';
							}
							else{
								#~ Passo. Fita minus
								$seqs{$id}{sense} = '-';
						}
					}
						#~ Passo. Armazene a linha no campo dna do hash
						$seqs{$id}{dna}= $line;
					}
					#~ Passo. Decremente o número de sequências
					$aux--;
				}
			}
		}
		#~ Passo. Feche a entrada 
		close IN;
		#~ Passo. Vá para o próximo arquivo
		$file++;
	}
	
	#~ foreach $line(sort{$seqs{$a}{number} <=> $seqs{$b}{number}} keys %seqs){
		#~ print "$seqs{$line}{sense}\n";
	#~ }
#~ ---------------------------------------------------------------------------------------
#~ 						Processamento e Escrita
#~ ---------------------------------------------------------------------------------------
	#~ Passo. Para cada sequência do hash
	foreach $line(keys %seqs){
	#~ print "$seqs{$line}{sense}\n";
		$start = -1;
		#~ Passo. Obtenha a posição inicial da substring numerada na sequência da proteína
		if($seqs{$line}{protein} =~ /$seqs{$line}{subcdr}/){
			$start = $-[0];
		}
		#~ Passo. Obtenha a sequência de dna correspondente a substring numerada
		if($start != -1){
			$i = ($start * 3) + ($seqs{$line}{frame} - 1);
			$len = length($seqs{$line}{subcdr}) * 3;
			#~ Passo. Se a fita codificante for plus
			if($seqs{$line}{sense} =~ /\+/){
				$seqs{$line}{cds} = substr($seqs{$line}{dna},$i,$len);
			}
			else{
				$seqs{$line}{minus} = $seqs{$line}{dna};
				$seqs{$line}{minus} =~ tr/AGCT/TCGA/;
				$seqs{$line}{minus} = reverse($seqs{$line}{minus});
				$seqs{$line}{cds} = substr($seqs{$line}{minus},$i,$len);
				$seqs{$line}{cds} =~ tr/AGCT/TCGA/;
				$seqs{$line}{cds} = reverse($seqs{$line}{cds});
			}	
		}
	}
	
	#~ Passo. Abra o arquivo de saída
	open OUT, ">$output" or die "get_ntsequences: Could not create output file $output: $!\n";
	#~ Passo. Para cada sequência do hash
	foreach $line(sort{$seqs{$a}{number} <=> $seqs{$b}{number}} keys %seqs){
		#~ Passo. Imprima o id
		print OUT ">$seqs{$line}{id}\n";
		#~ Passo. Se foi possível recuperar a sequência de nucleotídeos
		if(length($seqs{$line}{cds}) > 1){

			$i = 0;
			$len = length($seqs{$line}{cds});
			$lensub = $len;
			#~ Passo. Enquanto houverem caracteres da string cds
			while($i < $len){
				#~ Passo. Se tiver 70 caracteres
				if($lensub / 70 > 0){
					#~ Passo. Obtenha substring de 70 caracteres
					$aux = substr($seqs{$line}{cds},$i,70);
					$i += 70;
				}
				#~ Passo. Senão, obtenha substring de $i até o último posição
				else{
					$aux = substr($seqs{$line}{cds},$i,$len-$i);
					$i = $len;
				}
				#~ Passo. Imprima a substring
				print OUT "$aux\n";
				#~ Passo. Obtenha  o número de caracteres restantes
				$lensub = $len - length($aux);
			}
		}
		else{
			print OUT "Could not retrieve nucleotide sequence.";
			if($seqs{$line}{protein} =~ /X/){
				print OUT " Unknown base (N).\n";
			}
			else{
				print OUT "\n";
			}
		}
	}
	#~ Passo. Feche a saída
	close OUT;
	
	exit(0);
	

