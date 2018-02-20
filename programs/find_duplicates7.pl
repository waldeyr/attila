#! usr/bin/perl

use warnings;
use strict;

#----------------------------------------------------------------------------------------------------------------------
# Programa: find_duplicates6.pl
# Este programa recebe dois arquivos de entrada, ambos contendo uma lista
# de sequências de cdr com suas respectivas frequências e ids. O programa busca
# sequências cuja frequência aumentou do arquivo 1 em relação ao arquivo 2, e as
# imprime num arquivo de saída. Este programa volta a fazer uma coisa, que eu pensei que
# era gambiarra, mas em perl não é. Então a estrutura de dados onde são armazenadas as
# sequências é um vetor de hash, porém, as chaves dos hashes, diferente das 
# outras versões, são as próprias sequências das cdrs.Este programa imprime uma lista ordenada de
# sequências de cdrs de acordo com o fold change (freq2/freq1).
##----------------------------------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------------------------------------
# 					Início da função principal
#------------------------------------------------------------------------------------------------------------------------

	# Passo. Se o programa receber menos argumentos que o necessário
	if (@ARGV < 3)
	{
	# 	Passo. Pare a execução do programa e imprima mensagem de erro
		die "find_duplicates7.pl: perl find_duplicates7.pl <inputfile1> <inputfile2> <outputfile>\n";
	}

	# Passo. Declaração de variáveis 
	my ($input1,$input2,$output) = @ARGV;
	my ($line,$header,$freq,$number,$x,$y);
	my ($i,$c,$f,$counter,$libsize,$temp);
	my ($idseq,$seq,$cdr,$k,$frame);
	my @cdrs;
	my %allcdrs;

	# Passo. Inicialize o contador de sequências lidas
	$number = 0;
	# Passo. Inicialize o índice do vetor de hashes
	$i = 0;

#---------------------------------------------------------------------------------------------------------------------------
#					Leitura do arquivo 1
#---------------------------------------------------------------------------------------------------------------------------

	# Passo. Abra o arquivo de entrada
	open IN, "<$input1" or die "find_duplicates7.pl: Could not open input file $input1: $!\n";

	# Passo. Enquanto houver entrada
	while ($line = <IN>){
	# 	Passo. Remova o caracter \n do fim da linha
		chomp($line);

	# 	Passo.Se a linha começar com '#'
		if ($line =~ /^#(.+)\|FRAME:(\d)\|(.+)\|(\d+)\|((\d+)\.?(\d+))$/){		
	# 		Passo. Armazene a primeira expressão regular em $header
			$header = $1;
			#~ Passo. Armazene a segunda expressão regular como tamanho da biblioteca
			$libsize = $4;
	# 		Passo. Armazene a terceira expressão regular em $freq
			$freq = $5;
		}
		else{
	# 		Passo. Caso a linha comece com '*'
			if($line =~ /^\*(.+)/){
	# 			Passo. Armazene $freq no campo freq, usando a sequência como chave 
				$cdrs[$i]{$1}{freq1} = $freq;
	# 			Passo. Inicialize o campo freq2
				$cdrs[$i]{$1}{freq2} = 0;
	# 			Passo. Incremente o contador de Leitura
				$number++;
		# 		Passo. Se 10000 sequências já foram lidas
				if ($number % 10000 == 0){
	#		 		Passo. Vá para o próximo hash do vetor
					$i++;
				}
			}

		}

	}

	# Passo. Feche o arquivo 1
	close IN;

#---------------------------------------------------------------------------------------------------------------------------
#					Leitura do arquivo 2 e processamento
#---------------------------------------------------------------------------------------------------------------------------
	# Passo. Abra o arquivo 2
	open IN, "<$input2" or die "find_duplicates7.pl: Could not open input file $input2: $!\n";

	$k = -1;
	# Passo. Enquanto houver entrada
	while ($line = <IN>){
	# 	Passo. Remova o caracter \n do fim da linha
		chomp($line);
		
		# 	Passo.Se a linha começar com '#'
		if ($line =~ /^#(.+)\|FRAME:(\d)\|(.+)\|(\d+)\|((\d+)\.?(\d+))$/){		
	# 		Passo. Armazene a quinta expressão regular em $freq
			$freq = $5;
		}
		#~ Passo. Senão se a linha for da susbtring de cdrs
		elsif ($line =~ /^\*(.+)/){
			#~ Passo. Armazene a linha em cdr
			$cdr = $1;
	# 		Passo. Se $number é múltiplo de 10000
			if ($number % 10000 == 0){
	# 			Passo. Atribua como último índice $i-1
				$f = $i - 1;
			}
			else{
	# 			Passo. Atribua como último índice $i
				$f = $i;
			}
	#		Passo. Busque a cdr2 no hash cdrs
			$c = busca_cdr(\@cdrs,$1,\$f) ;
	# 					
	# 		Passo. Se a sequência do ciclo 3 existir no ciclo 1
			if($c != -1){
				#~ Passo. k recebe o índice da cdr
 				$k = $c;
			}
			else{
				#~ Passo. k recebe o índice atual do vetor de hashes
				$k = $i;
				#~ Passo. Armazene 1/libsize no campo freq1 usando a sequência como chave
				$cdrs[$k]{$1}{freq1} = (1 / $libsize) * 100000000000000 ;
#				Passo. Incremente o contador de leitura
				$number++;
			}
			
			#~ Passo. Armazene freq no campo freq2
			$cdrs[$k]{$1}{freq2} = $freq;
# 			Passo. Se a freq2 for maior que a frequência 1
			if ($cdrs[$k]{$1}{freq2} > $cdrs[$k]{$1}{freq1}){
# 				Passo. Armazene o fold change da sequência no hash dif
				$temp = $cdrs[$k]{$1}{freq2} / $cdrs[$k]{$1}{freq1};
				$allcdrs{$1}{dif} = sprintf("%.4f", $temp);
# 				Passo. Armazene o índice do vetor de hashes onde se encontra esta sequência
				$allcdrs{$1}{index} = $k;
				#~ Passo. Inicialize os campo id e seq de allcdrs
				$allcdrs{$1}{id} = "";
				$allcdrs{$1}{seq} = "";
				
			}
# 			Passo. Se $number é múltiplo de 10000
			if ($number % 10000 == 0){
#	 				Passo. Vá para o proximo elemento do vetor
					$i++;
			}
		}
		#~ Passo. Se a linha for do id de uma sequência completa
		elsif($line =~ /^>(.+)/){
				#~ Passo. Armazene linha em idseq
				$idseq = $1;
		}
		else{
			#~ Passo. Armazene linha em seq
			$seq = $line;
			if(exists $allcdrs{$cdr}){
				if(length($seq) > length($allcdrs{$cdr}{seq})){
					$allcdrs{$cdr}{id} = $idseq;
					$allcdrs{$cdr}{seq} = $seq;
				}
			}
		}
	}
	# Passo. Feche o arquivo 2
	close(IN);

	# Passo. Se $number é múltiplo de 10000
	if ($number % 10000 == 0){
	# 	Passo. Atribua $i -1 como último índice do vetor
		$f = $i - 1;
	}
	else{
	# 	Passo. Atribua $i como último índice do vetor
		$f = $i;
	}

#---------------------------------------------------------------------------------------------------------------------------
#						Escrita
#---------------------------------------------------------------------------------------------------------------------------

	# Passo. Abra o arquivo de saída
	open OUT, ">$output" or die "find_duplicates7.pl: Could not open output file $output\n";

	# Passo. Obtenha a lista ordenada de chaves do hash %allcdrs de acordo com o valor do hash dif 
	foreach $header (sort{$allcdrs{$b}{dif} <=> $allcdrs{$a}{dif}} keys %allcdrs){
	# 	Passo. Atribua o índice do vetor de hash à variavel temporária $c
		$c = $allcdrs{$header}{index};
		#~ $cdrs[$c]{$header}{freq1} = $cdrs[$c]{$header}{freq1} / 100000000000000 ;
		#~ $cdrs[$c]{$header}{freq2} = $cdrs[$c]{$header}{freq2} / 100000000000000 ;
	# 	Passo. Imprima o id, a frequência do ciclo 1, a frequência do ciclo 3 e o fold change da sequência atual
		print OUT ">$allcdrs{$header}{id}|FOLD-CHANGE:$allcdrs{$header}{dif}|P1:$cdrs[$c]{$header}{freq1}|P2:$cdrs[$c]{$header}{freq2}\n";
	# 	Passo. Imprima a sequência da cdrs
		print OUT "$allcdrs{$header}{seq}\n";
	}

	# Passo. Feche o arquivo de saída
	close OUT;

	exit (0);

#----------------------------------------------------------------------------------------------------------------------------------------
# Subrotina busca_cdr: Esta subrotina recebe um vetor de hashes, uma string, e o último índice cujo elemento foi preenchido
# no vetor. A subrotina busca no vetor de hashes uma chave que seja igual a string recebida. Se achar uma chave, uma variável
# auxiliar armazena o valor do índice do vetor onde está o hash com tal chave. A subrotina devolve a variável auxiliar, que
# terá valor igual a -1 caso não seja encontrada nenhuma chave igual à string, ou o valor do índice do vetor onde foi encontrada a chave.
#----------------------------------------------------------------------------------------------------------------------------------------
sub busca_cdr{

	# Passo. Receba as variáveis da main
	my $temp = shift or die "find_duplicates7.pl: Could not receive array of hashes\n";
	my @cdrs = @$temp;
	my $cdr2 = shift or die "find_duplicates7.pl: Could not receive string\n";
	$temp = shift  or die "find_duplicates7.pl: Could not receive last array index";
	my $f = $$temp;
	my ($c,$flag);

	# Passo. Inicialize flag
	$flag = -1;
	# Passo. Inicialize $c
	$c = 0;
	# Passo. Enquanto não achar a chave em um hash e enquanto não chegar ao fim do vetor
	while($flag == -1 && $c <= $f){
	# 	Passo. Se a chave existe no hash atual
		if (exists $cdrs[$c]{$cdr2}){
	# 		Passo. Altere flag
			$flag = $c;
		}
		else{
	# 		Passo. Vá para o próximo elemento do vetor
			$c++;
		}
	}

	# Passo. Retorne o índice encontrado
	return $flag;

}


