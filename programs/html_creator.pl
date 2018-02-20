#! /usr/bin/perl

	use warnings;
	use strict;
#~ ----------------------------------------------------------------------------------------------------------------------------------	

	#~ Programa: html_creator.pl
	#~ Data: 14/09/2015
	#~ Este programa recebe como entrada 6 arquivos:
	#~ arquivo 1: saída do Abnum em formato fasta da biblioteca VH
	#~ arquivo 2: saída do Abnum em formato fasta da biblioteca VL
	#~ arquivo 3: saída do IgBlast da biblioteca VH
	#~ arquivo 4: saída do IgBlast da biblioteca VL
	#~ arquivo 5: arquivo contendo links ncbi das germlines VH ***
	# arquivo 6: arquivo contendo links ncbi das germlines VL ***
	#~ Este programa lê e armazena as seguintes informações 
	#~ arquivos 1 e 2: id, sequência, fold change
	#~ arquivos 3 e 4: germline, valor de identidade, posições inicial e final de cada região do domínio variável
	#~ arquivos 5 e 6: links das germlines do NCBI.
	#~ Estas informações são escritas	 num arquivo de saída em formato html, o qual conterá duas abas, uma para dados da biblioteca VH
	# e outra para dados da biblioteca VL, em que cada aba contem duas tabelas, uma tabela referente a
	# classificação de germlines e outra tabela mostrando as regiões do domínio variável identificadas 
	# de acordo com o alinhamento com as germlines.
	#~ *** Atenção: no momento, o NCBI removeu os bancos de germlines humanos, e portanto os links não funcionam. Diante
	#~ disso os trechos do código referentes aos links serão comentados.
#~ ----------------------------------------------------------------------------------------------------------------------------------	

	#~ Passo. Declaração de variáveis
	my ($in1,$in2,$in3,$in4,$out,$dir1,$dir2,$dir3,$dir4,$csvfilevh,$csvfilevl,$statfilevh,$statfilevl) = @ARGV;
	my ($line,$counter,$id,$start,$end,$cdr3end,$fr4,$flag,$lib,$header);
	my (%table1,%table2,$input,$countvh,$countvl,$path1,$path2,$aux,%perda);
	my ($perdavhr0,$perdavhrn,$perdavlr0,$perdavlrn,%stat,$alfavh,$alfavl,$temp);
	$flag = -1;
	
	if(@ARGV < 13){
		die "perl html_creator.pl <input1.fasta> <input2.fasta> <input3.txt> <input4.txt> <output.html> <pathtoplot1> <pathtoplot2> <pathtoplot3> <pathtoplot4> <csvfilevh> <csvfilevl> <statvh.txt> <statvl.txt>\n";
	}

#~ ----------------------------------------------------------------------------------------------------------------------------------	
				#~ Leitura do arquivo fasta 
#~ ----------------------------------------------------------------------------------------------------------------------------------	
	# Passo. Inicialize o tipo de biblioteca e o contador de leitura
	$lib = 0;
	$counter = 0;
	# Passo. Enquanto houver arquivos fasta para serem lidos
	while($lib <= 1){
		# Passo. Se a biblioteca for VH
		if($lib == 0){
			# Passo. Selecione o arquivo fasta de VH
			$input = $in1;
			
		}
		else{
			# Passo. Selecione o arquivo fasta de VL
			$input = $in2;	
			
		}
		#~ Passo. Abra o arquivo fasta 
		open IN, "<$input" or die "Could not open input file $input: $!\n";
		#~ Passo. Enquanto houver entrada
		while(defined($line = <IN>)){
			chomp($line);
			#~ Passo. Se a linha for do id
			if($line =~ /^>(.+)\|FRAME:(\d)\|FOLD-CHANGE:(\d+\.\d+)$/){
				$table1{$1}{fc} = sprintf("%.2f", $3);
				$table1{$1}{count} = $counter;
				$table2{$1}{count} = $counter;
				$counter++;
			}
			else{
				$table2{$1}{seq} = $line;
				$table2{$1}{seq} =~ s/-//g;
			}
		}
		if ($lib == 0){
			$countvh = $counter;
		}
		else{
			$countvl = $counter;
		}
		
		close IN;
		$lib++;
	}
	

#~ ----------------------------------------------------------------------------------------------------------------------------------	
		#~ Leitura da saída do IgBlast
#~ ----------------------------------------------------------------------------------------------------------------------------------	
	# Passo. Inicialize o tipo de biblioteca
	$lib = 0;
	# Passo. Enquanto houver arquivos do IgBlast para serem lidos
	while($lib <= 1){
		# Se a biblioteca for VH
		if($lib == 0){
			# Passo. Selecione a saída do IgBlast de VH
			$input = $in3;
		}
		else{
			# Passo. Selecione a saída do IgBlast de VL
			$input = $in4;
		}
		# Passo. Abra o arquivo de entrada
		open IN, "<$input" or die "Could not open input file $input: $!\n";
	
		#~ Passo. Enquanto houver entrada
		while(defined($line = <IN>)){
			chomp($line);
			#~ Passo. Se a linha for do id
			if($line =~ /^#\s+Query:\s+(.+)\|FRAME:(\d)\|FOLD-CHANGE:(\d+\.\d+)$/){
				#~ Passo. Armazene o id numa string temporária
				$id = $1;
				$flag = 0;
			}
			#~ Passo. Se a linha for da FR1 
			elsif($line =~ /^FR1\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+((\d+(\.)?\d+)?)/){
				#~ Passo. Armazene as posições inicial e final da região e obtenha a substring correspondente
				$start = $1 - 1;
				$end = $2;
				$table2{$id}{fr1} = substr($table2{$id}{seq},$start,($end - $start));
			}
			#~ Passo. Se a linha for da CDR1 
			elsif($line =~ /^CDR1\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+((\d+(\.)?\d+)?)/){
				#~ Passo. Armazene as posições inicial e final da região e obtenha a substring correspondente
				$start = $1 - 1;
				$end = $2;
				$table2{$id}{cdr1} = substr($table2{$id}{seq},$start,($end - $start));
				
			}
			#~ Passo. Se a linha for da FR2 
			elsif($line =~ /^FR2\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+((\d+(\.)?\d+)?)/){
				#~ Passo. Armazene as posições inicial e final da região e obtenha a substring correspondente
				$start = $1 - 1;
				$end = $2;
				$table2{$id}{fr2} = substr($table2{$id}{seq},$start,($end - $start));
				
			}
			#~ Passo. Se a linha for da CDR2 
			elsif($line =~ /^CDR2\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+((\d+(\.)?\d+)?)/){
				#~ Passo. Armazene as posições inicial e final da região e obtenha a substring correspondente 
				$start = $1 - 1;
				$end = $2;
				$table2{$id}{cdr2} = substr($table2{$id}{seq},$start,($end - $start));
				
			}
			#~ Passo. Se a linha for da FR3
			elsif($line =~ /^FR3\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+((\d+(\.)?\d+)?)/){
				#~ Passo. Armazene as posições inicial e final da região e obtenha a substring correspondente 
				$start = $1 - 1;
				$end = $2;
				$table2{$id}{fr3} = substr($table2{$id}{seq},$start,($end - $start));
				#~ Passo. Obtenha a posição final da CDR3
				if($table2{$id}{seq} =~ /(.+)WG\wG/ || $table2{$id}{seq} =~ /(.+)FG\wG/){
					$cdr3end = $+[1];
					#~ Passo. Armazene no campo cdr3 a substring entre a posição final anterior ($end) e a posição do regex atual ($cdr3end)
					$table2{$id}{cdr3} = substr($table2{$id}{seq},$end,($cdr3end - $end));
					#~ Passo. Armazene no campo fr4 a substring que vai de $cdr3end até o último caracter da sequência
					$fr4 = length($table2{$id}{seq});			
					$table2{$id}{fr4} = substr($table2{$id}{seq},$cdr3end,($fr4 - $cdr3end));
				}
			}
			else{
				if($flag == 0){
				#~ Passo. Se a linha for dos valores de identidade
					if($line =~ /^V(.+)FOLD-CHANGE:(\d+\.\d+)\s+(\S+)?\s+(\d+(\.\d+)?)?\s+(\S+)/){
						#~ Passo. Armazene o nome da germline com maior valor de identidade
						$table1{$id}{germline} = $3;
						#~ Passo. Armazene o valor da identidade
						$table1{$id}{identity} = $4 + 0;
						#~ Passo. Altere flag para indicar que já foi encontrada a germline com maior identidade
						$flag = 1;
					}
				}
			}
		}
		close IN;
		$lib++;
	}
#~ ----------------------------------------------------------------------------------------------------------------------------------	
					# Leitura do arquivo de links ncbi
#~ ----------------------------------------------------------------------------------------------------------------------------------	
	#~ # Passo. Inicialize o tipo de biblioteca
	#~ $lib = 0;
	#~ # Passo. Enquanto houver arquivos de links ncbi para serem lidos
	#~ while($lib <= 1){
		#~ # Se a biblioteca for VH
		#~ if($lib == 0){
			#~ # Passo. Selecione a saída do IgBlast de VH
			#~ $input = $in5;
		#~ }
		#~ else{
			#~ # Passo. Selecione a saída do IgBlast de VL
			#~ $input = $in6;
		#~ }
		#~ # Passo. Para cada sequência
		#~ foreach $line (keys %table1){
			#~ # Passo. Busque e armazene o link ncbi da germline
			#~ $table1{$line}{link} = `grep -P ">$table1{$line}{germline}</a>" $input`;
			#~ chomp($table1{$line}{link});
		#~ }
		#~ $lib++;
	#~ }

#~ ----------------------------------------------------------------------------------------------------------------------------------	
					# Leitura do arquivo csv
#~ ----------------------------------------------------------------------------------------------------------------------------------	
	# Passo. Inicialize o tipo de biblioteca
	$lib = 0;
	$aux = 0;
	# Passo. Enquanto houver arquivos csv para serem lidos
	while ($lib <= 1){
		# Passo. Se a biblioteca for VH
		if ($lib == 0){
			$input = $csvfilevh;
			$aux = "vh";
		}
		else{
			$input = $csvfilevl;
			$aux = "vl";
		}
		# Passo. Abra o arquivo
		open IN, "<$input" or die "Could not open csv file $input: $!\n";
		while ($line=<IN>){
			chomp($line);
			# Passo. Se a linha for da biblioteca R0 antes da filtragem
			if ($line =~ /^R0,(\d+),raw$/){
				$perda{$aux}{r0}{raw} = $1;
			}
			# Passo. Se a linha for da biblioteca R0 filtrada
			if ($line =~ /^R0,(\d+),filtering$/){
				$perda{$aux}{r0}{filtering} = $1;
			}
			# Passo. Se a linha for da biblioteca RN antes da fitragem
			if ($line =~ /^RN,(\d+),raw$/){
				$perda{$aux}{rn}{raw} = $1;
			}
			# Passo. Se a linha for da biblioteca RN filtrada
			if ($line =~ /^RN,(\d+),filtering$/){
				$perda{$aux}{rn}{filtering} = $1;
			}
		}
		close IN;
		$lib++;
	}

	# Passo. Calcule o percentual de perda de VH R0
	$perdavhr0 = ($perda{vh}{r0}{raw} - $perda{vh}{r0}{filtering}) / $perda{vh}{r0}{raw};
	$perdavhr0 = $perdavhr0 * 100;
	$perdavhr0 = sprintf("%.2f", $perdavhr0);

	# Passo. Calcule o percentual de perda de VH RN
	$perdavhrn = ($perda{vh}{rn}{raw} - $perda{vh}{rn}{filtering}) / $perda{vh}{rn}{raw};
	$perdavhrn = $perdavhrn * 100;
	$perdavhrn = sprintf("%.2f", $perdavhrn);
	
	# Passo. Calcule o percentual de perda de VL R0
	$perdavlr0 = ($perda{vl}{r0}{raw} - $perda{vl}{r0}{filtering}) / $perda{vl}{r0}{raw};
	$perdavlr0 = $perdavlr0 * 100;
	$perdavlr0 = sprintf("%.2f", $perdavlr0);

	# Passo. Calcule o percentual de perda de VL RN
	$perdavlrn = ($perda{vl}{rn}{raw} - $perda{vl}{rn}{filtering}) / $perda{vl}{rn}{raw};
	$perdavlrn = $perdavlrn * 100;
	$perdavlrn = sprintf("%.2f", $perdavlrn);

		

#~ ----------------------------------------------------------------------------------------------------------------------------------	
	#				Leitura dos arquivos dos testes estatísticos
#~ ----------------------------------------------------------------------------------------------------------------------------------	
	# Passo. Enquanto houver arquivos para serem lidos
	$lib = 0;
	$aux = "";
	while ($lib <= 1){
		# Passo. Se a biblioteca for de Vh
		if ($lib == 0){
			$input = $statfilevh;
			$aux = "vh";
		}
		else{
			$input = $statfilevl;
			$aux = "vl";	
		}
		# Passo. Abra o arquivo 
		open IN, "<$input" or die "Could not open input file $input: $!\n ";

		# Passo. Enquanto houver entrada
		while ($line=<IN>){
			chomp($line);
			# Passo. Se a linha for do alfa de Bonferroni
			if ($line =~ /^alfaBonf:(\s+)(\d.\d+)$/){
				$temp = $2;
			}
			# Passo. Se a linha for do id
			if ($line =~ /^>(.+)/){
				$header = $1;
			}
			# Passo. Se a linha for do p-value
			if ($line =~ /^p\-value:(\s+)(\d*\.*\d*e*\-*\d*)$/){
				$stat{$header}{pvalue} = $2;
			}
			# Passo. Se a linha for do infC
			if ($line =~ /^infIC:(\s+)(\-*\d*\.*\d*e*\-*\d*)$/){
				$stat{$header}{infIC} = $2;
			}
			# Passo. Se a linha for do supC
			if ($line =~ /^supIC:(\s+)(\-*\d*\.*\d*e*\-*\d*)$/){
				$stat{$header}{supIC} = $2;
			}
		}
		if($lib == 0){
			$alfavh = $temp;
			$temp = "";
		}
		else{
			$alfavl = $temp;
		}
		close IN;
		$lib++;
	}
	
#~ ----------------------------------------------------------------------------------------------------------------------------------	
					#~ Escrita do arquivo html
#~ ----------------------------------------------------------------------------------------------------------------------------------	

	# Passo. Abra o arquivo html a ser escrito
	open(OUT, '>:encoding(UTF-8)', $out) or die "Could not create output file: $!\n";
	
	# Passo. Imprima instrução DOCTYPE e css da página
	print OUT"<DOCTYPE! html>\n";
	print OUT"<html>\n";
	print OUT"<head>\n";
	print OUT"<meta charset=\"utf-8\">\n";
	print OUT"<link rel=\"stylesheet\"href=\"http://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap.min.css\">\n";
	print OUT"<script src=\"https://ajax.googleapis.com/ajax/libs/jquery/1.12.4/jquery.min.js\"></script>\n";
	print OUT"<script src=\"http://maxcdn.bootstrapcdn.com/bootstrap/3.3.5/js/bootstrap.min.js\"></script>\n";
	print OUT"<title>Selected sequences</title>\n";
	print OUT"<style>\n";
	print OUT"body{\n";
	print OUT"	font-family: \"Verdana\", sans-serif;\n";
	print OUT"	text-align: center;\n";
	print OUT"}\n";
	print OUT".container-fluid{\n";
	print OUT"	max-width: 1600px;\n";
	print OUT"}\n";
	print OUT"header{\n";
	print OUT"	width: 100%;\n";
	print OUT"	height: 125px;\n";
	print OUT"	background-color: #0039a6;\n";
	print OUT"	border-radius: 5px;\n";
	print OUT"	}\n";
	print OUT"h1{\n";
	print OUT"	text-align:left;\n";
	print OUT"	border-bottom-style: solid;\n";
	print OUT"	border-bottom-color: #A9A9A9;\n";
	print OUT"	border-bottom-width: 3px;\n";
	print OUT"	font-size: 20px;\n";
	print OUT"	padding-top: 10px;\n";
	print OUT"	font-weight: normal;\n";
	print OUT"}\n";
	print OUT"h2{\n";
	print OUT"	text-align:left;\n";
	print OUT"	border-bottom-style: solid;\n";
	print OUT"	border-bottom-color: #A9A9A9;\n";
	print OUT"	border-bottom-width: 3px;\n";
	print OUT"	font-size: 20px;\n";
	print OUT"	font-weight: normal;\n";
	print OUT"}\n";
	print OUT"h3{\n";
	print OUT"text-align:left;\n";
	print OUT"border-bottom-style: solid;\n";
	print OUT"border-bottom-color: #A9A9A9;\n";
	print OUT"border-bottom-width: 3px;\n";
	print OUT"font-size: 20px;\n";
	print OUT"font-weight: normal;\n"; 
	print OUT"}\n";
	print OUT"h4{\n";
	print OUT"text-align:left;\n";
	print OUT"font-size: 20px;\n";
	print OUT"font-weight: normal;\n";
	print OUT"color:#EE3B3B;\n";
	print OUT"}\n";
	print OUT".tablegerm{\n";
	print OUT"	width: 70%;\n";
	print OUT"	height: 200px;\n";
	print OUT"	table-layout: fixed;\n";
	print OUT"	font-size: 13px;\n";
	print OUT"	text-align: center;\n";
	print OUT"	margin-left: auto;\n";
	print OUT"	margin-right: auto;\n";
	print OUT"	background-color: #F5F5F5;\n";
	print OUT"}\n";
	print OUT".tabledomain{\n";
	print OUT"	height: 400px;\n";
	print OUT"	width: 100%;\n";
	print OUT"	table-layout: fixed;\n";
	print OUT"	font-size: 10px;\n";
	print OUT"	text-align: center;\n";
	print OUT"}\n";
	print OUT".germlinerow {\n";
	print OUT"	font-weight: bold;\n";
	print OUT"	color: white;\n";
	print OUT"	background-color: #708090; 	\n";
	print OUT"}\n";
	print OUT".domainrow{\n";
	print OUT"	height: 40px;\n";
	print OUT"	font-weight: bold;\n";
	print OUT"	color: white;\n";
	print OUT"	font-size: 13px;\n";
	print OUT"}\n";
	print OUT"#ID{\n";
	print OUT"	width: 10%;\n";
	print OUT"	\n";
	print OUT"}\n";
	print OUT"#fc{\n";
	print OUT"	width: 15%;\n";
	print OUT"	\n";
	print OUT"}\n";
	print OUT"#germ{\n";
	print OUT"	width: 15%;\n";
	print OUT"}\n";
	print OUT"#identity{\n";
	print OUT"	width: 10%;\n";
	print OUT"}\n";
	print OUT"#CDR1{\n";
	print OUT"color: #3399FF;\n";
	print OUT"background-color:  #F5F5F5;\n";
	print OUT"width: 8%;\n";
	print OUT"}\n";
	print OUT"#FR2{\n";
	print OUT"color: #DB7093;\n";
	print OUT"background-color: #F5F5F5;\n";
	print OUT"width: 10%;\n";
	print OUT"}	\n";
	print OUT"#FR1{\n";
	print OUT"color: #FF6347;\n";
	print OUT"background-color: #F5F5F5;\n";
	print OUT"width: 16%;\n";
	print OUT"}\n";
	print OUT"#CDR2{\n";
	print OUT"color: 	#DAA520;\n";
	print OUT"background-color: #F5F5F5;\n";
	print OUT"width: 10%;\n";
	print OUT"}\n";
	print OUT"#FR3{\n";
	print OUT"color: #3CB371;\n";
	print OUT"background-color: #F5F5F5;\n";
	print OUT"width: 20%;\n";
	print OUT"}\n";
	print OUT"#CDR3{\n";
	print OUT"color: #FF8C00;\n";
	print OUT"background-color:  	#F5F5F5;\n";
	print OUT"width: 7%;\n";
	print OUT"}\n";
	print OUT"#FR4{\n";
	print OUT"color: #9999FF;\n";
	print OUT"background-color:  	#F5F5F5;\n";
	print OUT"width: 10%;\n";
	print OUT"}\n";
	print OUT"#BFR1{\n";
	print OUT"background-color: #FF6347;\n";
	print OUT"width: 15%;\n";
	print OUT"}\n";
	print OUT"#BCDR1{\n";
	print OUT"background-color: #3399FF;\n";
	print OUT"width: 11%;\n";
	print OUT"}\n";
	print OUT"#BFR2{\n";
	print OUT"background-color: #DB7093;\n";
	print OUT"width: 10%;\n";
	print OUT"}\n";
	print OUT"#BCDR2{\n";
	print OUT"background-color:  	#DAA520;\n";
	print OUT"width: 10%;\n";
	print OUT"}\n";
	print OUT"#BFR3{\n";
	print OUT"background-color: #3CB371;\n";
	print OUT"width: 20%;\n";
	print OUT"}\n";
	print OUT"#BCDR3{\n";
	print OUT"background-color: #FF8C00;\n";
	print OUT"width: 8%;\n";
	print OUT"}\n";
	print OUT"#BFR4{\n";
	print OUT"background-color: #9999FF;\n";
	print OUT"width: 10%;\n";
	print OUT"}\n";
	print OUT"#idread{\n";
	print OUT"	background-color: #696969;\n";
	print OUT"	width: 5%;\n";
	print OUT"	\n";
	print OUT"}\n";
	print OUT"#idcolumn{\n";
	print OUT"	background-color: #F5F5F5;\n";
	print OUT"	width: 5%;\n";
	print OUT"	\n";
	print OUT"}\n";
	print OUT"#linha{\n";
	print OUT"	border-bottom-style: solid;\n";
	print OUT"	border-bottom-color: #A9A9A9;		\n";
	print OUT"	border-bottom-width: 3px;\n";
	print OUT"	height: 20px;\n";
	print OUT"	\n";
	print OUT"}\n";
	print OUT"footer{\n";
	print OUT"	padding-top: 10px;\n";
	print OUT"	padding-bottom: 10px;\n";
	print OUT"	height: 90px;\n";
	print OUT"	width: 100%;\n";
	print OUT"}\n";
	print OUT"#application{\n";
	print OUT"	width: 100%;\n";
	print OUT"	float: left; \n";
	print OUT"	padding-top: 30px;\n";
	print OUT"	\n";
	print OUT"}\n";
	print OUT"#bioinfo{\n";
	print OUT"width: 7%;\n";
	print OUT"}\n";
	print OUT"#address{\n";
	print OUT"	text-align: left;\n";
	print OUT"	font-size: 15px;\n";
	print OUT"	color: #0039a6;\n";
	print OUT"	width: 40%;\n";
	print OUT"}\n";
	print OUT"#unb{\n";
	print OUT"	float: right;\n";
	print OUT"	width: 60%;\n";
	print OUT"}\n";
	print OUT"#grafico1{\n";
	print OUT"width: 45%;\n";
	print OUT"float: left;\n";
	print OUT"}\n";
	print OUT"#grafico2{\n";
	print OUT"width: 51%;\n";
	print OUT"float: right;\n";
	print OUT"}\n";
	print OUT"#tablefooter{\n";
	print OUT"	table-layout: fixed;\n";
	print OUT"}\n";
	print OUT"#minidiv{\n";
	print OUT"	height: 10px;\n";
	print OUT"}\n";
	print OUT".popover{\n";
	print OUT"	max-width: 100%;\n";
	print OUT"}	\n";
	print OUT"#variabledom{\n";
	print OUT"height: 200 px;\n";
	print OUT"}\n";
	print OUT"#pvalue{\n";
	print OUT"widht: 5%;\n";	
	print OUT"}\n";
	print OUT"#InfIC{\n";
	print OUT"widht: 15%;\n";		
	print OUT"}\n";
	print OUT"#supIC{\n";
	print OUT"widht: 15%;\n";		
	print OUT"}\n";
	print OUT"#lossreads{\n";
	print OUT"height: 820px;\n";	
	print OUT"}\n";
	print OUT"#lossreads1{\n";
	print OUT"height: 820px;\n";	
	print OUT"}\n";
	print OUT"</style>\n";
	print OUT"</head>\n";
	print OUT"<body>\n";
	print OUT"<div class=\"container-fluid\">\n";
	print OUT"<header class=\"container-fluid\">\n";
	print OUT"<div id=\"application\" >\n";
	print OUT"<img src=\"app5.png\"; width=\"70%\"> \n";
	print OUT"</div>\n";
	print OUT"</header>\n";
	print OUT"<script>\n";
	print OUT"\$(document).ready(function(){\n";
	print OUT"\$(\'[data-toggle=\"popover\"]\').popover();   \n";
	print OUT"});\n";
	print OUT"</script>\n";
	print OUT"<div class=\"container-fluid\">\n";
	print OUT"<ul class=\"nav nav-tabs\">\n";
	print OUT"<li class=\"active\"><a data-toggle=\"tab\" href=\"#vh\">VH Library</a></li>\n";
	print OUT"<li><a data-toggle=\"tab\" href=\"#vl\">VL Library</a></li>\n";
	print OUT"</ul>\n";
	  
	  
	
	print OUT"<div class=\"tab-content\">\n";
	

	$lib = 0;
	my $position;
	while ($lib <= 1){	
		if ($lib == 0){
			print OUT"<div id=\"vh\" class=\"tab-pane fade in active\">\n";
			$start = 0;
			$end = $countvh - 1;
			$path1 = $dir1;
			$path2 = $dir2;
			$input = $in1;
			
		}
		else{
			print OUT"<div id=\"vl\" class=\"tab-pane fade\">\n";
			$start = $countvh;
			$end = $countvl - 1;
			$path1 = $dir3;
			$path2 = $dir4;
			$input = $in2;
		}

		# Passo. Caso existam clones candidatos da biblioteca atual		
		$aux = 0;
		$aux = `grep -cP "^>" $input`;
		chomp($aux);

		$counter = 0;

		# Escrita da seção Reads Information
		if($lib == 0){
			print OUT"\t\t\t<h1><a data-toggle=\"collapse\" href=\"#lossreads\">Reads Information</a></h1>\n";		
			print OUT"\t\t\t<div id=\"lossreads\" class=\"panel-collapse collapse in\">\n";						
		}
		else{
			print OUT"\t\t\t<h1><a data-toggle=\"collapse\" href=\"#lossreads1\">Reads Information</a></h1>\n";		
			print OUT"\t\t\t<div id=\"lossreads1\" class=\"panel-collapse collapse in\">\n";						
		}
		print OUT"\t\t\t<div id=\"minidiv\"></div>\n";
		print OUT"\t\t\t<table id=\"tablegerm\" class=\"table table-bordered tablegermline\">\n";
		print OUT"\t\t\t\t<tr class=\"germlinerow\">\n";
		print OUT"\t\t\t\t\t<td>Dataset</td>\n";
		print OUT"\t\t\t\t\t<td>Loss of Reads (Initial Library)</td>\n";
		print OUT"\t\t\t\t\t<td>Loss of Reads (Final Library)</td>\n";
		print OUT"\t\t\t\t</tr>\n";
		print OUT"\t\t\t\t<tr>\n";
		if ($lib == 0){
			print OUT"\t\t\t\t\t<td>VH</td>\n";
			print OUT"\t\t\t\t\t<td>$perdavhr0%</td>\n";
			print OUT"\t\t\t\t\t<td>$perdavhrn%</td>\n";
		}
		else{
			print OUT"\t\t\t\t\t<td>VL</td>\n";
			print OUT"\t\t\t\t\t<td>$perdavlr0%</td>\n";
			print OUT"\t\t\t\t\t<td>$perdavlrn%</td>\n";
		}
		print OUT"\t\t\t\t</tr>\n";
		print OUT"\t\t\t</table>\n";	

		print OUT"\t\t\t<div id=\"minidiv\"></div>\n";
		print OUT"\t\t\t\t<img id=\"grafico1\"src=\"$path1\">\n";
		print OUT"\t\t\t\t<img id=\"grafico2\"src=\"$path2\">\n";
		print OUT"\t\t\t</div>\n";


		
		if ($aux > 0){
		#~ Escrita da seção Candidate Clones

		if ($lib == 0){
			print OUT"<h2><a data-toggle=\"collapse\" href=\"#candidates\">Candidate Clones</a></h2>\n";
			print OUT"<div id=\"candidates\" class=\"panel-collapse collapse in\">\n";		
		}
		else{
			print OUT"<h2><a data-toggle=\"collapse\" href=\"#candidates1\">Candidate Clones</a></h2>\n";
			print OUT"<div id=\"candidates1\" class=\"panel-collapse collapse in\">\n";		
		}
		print OUT"<div id=\"minidiv\"></div>\n";
		print OUT"\t\t\t<table id='tablegerm' class='table table-bordered tablegermline'>\n";
		print OUT"\t\t\t\t<tr class='germlinerow'>\n";
		print OUT"\t\t\t\t\t<td id='ID'>Sequence ID</td>\n";
		print OUT"\t\t\t\t\t<td id='fc'>Fold Change</td>\n";
		if ($lib == 0){
			print OUT"\t\t\t\t\t<td id='pvalue'>p-value (&alpha;Bonf=$alfavh)</td>\n";
		}
		else{
			print OUT"\t\t\t\t\t<td id='pvalue'>p-value (&alpha;Bonf=$alfavl)</td>\n";
		}
		print OUT"\t\t\t\t\t<td id='infic'>InfIC (IC=95%)</td>\n";
		print OUT"\t\t\t\t\t<td id='supic'>SupIC (IC=95%)</td>\n";
		print OUT "\t\t\t\t\t<td id='germ'>Germline</td>\n";
		print OUT "\t\t\t\t\t<td id ='identity'>Identity (%)</td>\n";
		print OUT "\t\t\t\t</tr>\n";

		foreach $line (sort{$table1{$a}{count} <=> $table1{$b}{count}} keys %table1){
			if($table1{$line}{count} >= $start && $table1{$line}{count} <= $end){
				$counter = $counter + 1;
				$position = $counter . "\xB0" ;	
				print OUT "\t\t\t\t\t<td><a href='#tablegerm' data-toggle='popover' data-trigger='hover' data-content='$line'>$position</a></td>\n";
				print OUT "\t\t\t\t\t<td>$table1{$line}{fc}</td>\n";
				print OUT "\t\t\t\t\t<td>$stat{$line}{pvalue}</td>\n";
				print OUT "\t\t\t\t\t<td>$stat{$line}{infIC}</td>\n";
				print OUT "\t\t\t\t\t<td>$stat{$line}{supIC}</td>\n";
				print OUT "\t\t\t\t\t<td>$table1{$line}{germline}</td>\n";
				print OUT "\t\t\t\t\t<td>$table1{$line}{identity}</td>\n";
				print OUT "\t\t\t\t</tr>\n";
			}
			$position = "";
		}
		print OUT "\t\t\t</table>\n";
		print OUT "\t\t</div>\n";

		# Escrita da seção Regions of Variable Domain of Candidate Clones
		print OUT"\t\t<h3><a data-toggle=\"collapse\" href=\"#variabledom\">Regions of Variable Domain of Candidate Clones</a></h3>\n";
		print OUT"\t\t<div id=\"variabledom\" class=\"panel-collapse collapse in\">\n";
		print OUT"\t\t<div id=\"minidiv\"></div>\n";
	
		$counter = 0;	
		print OUT "\t\t\t\t<table id='tabledom' class='table table-bordered tabledomain'>\n";
		print OUT "\t\t\t\t<tr class='domainrow'>\n";	
		print OUT "\t\t\t\t\t<td id='idread'>Sequence ID</td>\n";
		print OUT "\t\t\t\t\t<td id='BFR1'>FR1</td>\n";
		print OUT "\t\t\t\t\t<td id='BCDR1'>CDR1</td>\n";
		print OUT "\t\t\t\t\t<td id='BFR2'>FR2</td>\n";
		print OUT "\t\t\t\t\t<td id='BCDR2'>CDR2</td>\n";
		print OUT "\t\t\t\t\t<td id='BFR3'>FR3</td>\n";
		print OUT "\t\t\t\t\t<td id='BCDR3'>CDR3</td>\n";
		print OUT "\t\t\t\t\t<td id='BFR4'>FR4</td>\n";
		print OUT "\t\t\t\t</tr>\n";
		
				foreach $line (sort{$table2{$a}{count} <=> $table2{$b}{count}} keys %table2){
			if($table2{$line}{count} >= $start && $table2{$line}{count} <= $end){
				$counter = $counter + 1;
				$position = $counter . "\xB0" ;
				print OUT "\t\t\t\t<tr class='rows'>\n";
				print OUT "\t\t\t\t\t<td id=\"idcolumn\"><a href='#tabledom' data-toggle='popover' data-trigger='hover' data-content='$line'>$position</a></td>\n";
				print OUT "\t\t\t\t\t<td id='FR1'>$table2{$line}{fr1}</td>\n";
				print OUT "\t\t\t\t\t<td id='CDR1'>$table2{$line}{cdr1}</td>\n";
				print OUT "\t\t\t\t\t<td id='FR2'>$table2{$line}{fr2}</td>\n";
				print OUT "\t\t\t\t\t<td id='CDR2'>$table2{$line}{cdr2}</td>\n";
				print OUT "\t\t\t\t\t<td id='FR3'>$table2{$line}{fr3}</td>\n";
				print OUT "\t\t\t\t\t<td id='CDR3'>$table2{$line}{cdr3}</td>\n";
				print OUT "\t\t\t\t\t<td id='FR4'>$table2{$line}{fr4}</td>\n";
				print OUT "\t\t\t\t</tr>\n";
			}
			$position = "";
		}
		print OUT "\t\t\t\t</table>\n";
		print OUT "\t\t</div>\n";
		print OUT "\t\t</div>\n";
	
		}
		else{
			print OUT"<h4>No candidate clones found</h4>\n";  
		}
		
		$lib++;
	}
	
	print OUT "\t\t</div>\n";
	print OUT "\t\t</div>\n";
		
	#~ Escrita do footer	
	print OUT"<div id=\"linha\"></div>\n";
	print OUT"<footer >\n";
	print OUT"<table id=\"tablefooter\">\n";
	print OUT"<tr>\n";
	print OUT"<td id=\"bioinfo\">\n";
	print OUT"<img src=\"bioinfo.png\" width=\"80%\">\n";
	print OUT"</td>\n";
	print OUT"<td id=\"address\" >\n";
	print OUT"Universidade de Brasília<br/>\n";
	print OUT"Instituto de Ciências Biológicas<br/>\n";
	print OUT"Laboratório de Bioinformática, Bloco K<br/>\n";
	print OUT"CEP 70.000-00<br/>\n";
	print OUT"Brasília DF\n";
	print OUT"</td>\n";
	print OUT"<td id=\"unb\">\n";
	print OUT"<img id=\"unb\"src=\"unb_logo1.png\" alt=\"UNB\" width=\"85%\">\n";
	print OUT"</td>\n";
	print OUT"</tr>\n";
	print OUT"</table>\n";
	print OUT"</footer>\n";
	print OUT "\t</body>\n";
	print OUT "</html>\n";

	close OUT;
	exit(0);
