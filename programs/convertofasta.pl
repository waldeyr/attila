#! usr/bin/perl

use warnings;
use strict;
#--------------------------------------------------------------------------------------------------------------
# Programa:convertofasta.pl
# Este programa recebe como entrada um arquivo que contem uma lista de sequências dotadas
# de id e sequência. A sequência se encontra em formato de colunas numeradas, em que cada 
# resíduo de aminoácido recebe um número de acordo com o esquema de numeração de imunoglobulinas
# de Kabat. O programa concatene cada um dos caracteres que representam os aminoácidos para formar 
# uma sequência de domínio variável completa e imprime a sequência em formato fasta.
#--------------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------------------
#				Início da função principal
#--------------------------------------------------------------------------------------------------------------

# Passo. Declaração de variáveis
my ($line,$number,$aux1);
my %cdrs;
# Passo. Se o usuário entrar com menos argumentos que o necessário
if (@ARGV < 2)
{
# 	Passo. Imprima mensagem de erro e pare a execução do programa	
	die "convertfasta.pl: perl convertofasta.pl <inputfile> <outputfile> \n";
}
my ($input,$output) = @ARGV;

#--------------------------------------------------------------------------------------------------------------
#					Leitura
#--------------------------------------------------------------------------------------------------------------

# Passo. Abra o arquivo de entrada
open IN, "<$input" or die "convertfasta.pl: Could not open input file $input:$!\n";

# Passo. Inicialize o contador da ordem de leitura
$number = 0;

# Passo. Enquanto houver entrada
while ($line = <IN>)
{
# 	Passo. Remova o "\n" da linha
	chomp($line);
# 	Passo. Se a linha for do id da sequência
	if ($line =~ /^>(.+)/)
	{
# 		Passo. Armazene o id em $aux1
		$aux1 = $1;
# 		Passo. Inicialize o hash sequência usando o id como chave
		$cdrs{$aux1}{seq} = "";
# 		Passo. Atribua o contador de leitura para o hash number
		$cdrs{$aux1}{number} = $number;
# 		Passo. Incremente o contador de leitura	
		$number++;
	}
# 	Passo. Senão
	else 
	{
# 		Passo. Se a linha for de uma coluna de resíduo numerado
		if ($line =~ /([A-Z]?)$/)
		{
# 			Passo. Concatene o último caracter da linha no campo seq
			if(defined($1)){
				$cdrs{$aux1}{seq} .= $1;
			}
		}
	}
}

# Passo. Feche o arquivo de saída
close (IN);

#--------------------------------------------------------------------------------------------------------------
#					Escrita
#--------------------------------------------------------------------------------------------------------------
# Passo. Abra o arquivo de saída
open OUT, ">$output" or die "convertfasta.pl: Could not open output file $output:$!\n";

# Para cada id do hash
foreach $line (sort{$cdrs{$a}{number} <=> $cdrs{$b}{number}} keys %cdrs)
{
# 	Passo. Imprima o id
	print OUT ">$line\n";
# 	Passo. Imprima a sequência
	print OUT "$cdrs{$line}{seq}\n";
}

# Passo. Feche o arquivo
close (OUT);

exit(0);
		
