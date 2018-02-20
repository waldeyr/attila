#! usr/bin/perl

	use warnings;
	use strict;

#~---------------------------------------------------------------------------------------------
#~ Programa: autoiganalysis.pl
#~ Data: 18/12/2015
#~ Este programa automatiza a execução de uma série de etapas,
#~ de um pipeline de seleção in silico de sequências de imunoglobulinas,
#~ produzidas por tecnologia de phage display. A entrada do programa é 
#~ um arquivo de configuração,contendo argumentos usados nas etapas do 
#~ pipeline. O programa imprime um arquivo em formato csv, contendo o número
#~ de sequências das bibliotecas a cada etapa do pipeline. Além disso, são impressos
#~ dois arquivos de log, um deles contendo erro padrão do próprio programa automatic_selection.pl,
#~ e outro contendo erro padrão dos programas do pipeline. 
#~ Atenção: o arquivo de log dos erros de programas do pipeline conterá informação adicional, 
#~ isto é, mensagens do programa nohup, pois o nohup redireciona suas mensagens para STDERR, uma 
#~ vez que seja redirecionado o erro padrão para um arquivo diferente de nohup.out.
#~---------------------------------------------------------------------------------------------

#~---------------------------------------------------------------------------------------------
#~ 				       Função principal
#~---------------------------------------------------------------------------------------------

	#~ Passo. Declaração de variáveis
	my ($line,$input1dir,$input2dir,$aux,$gencodedir);
	my ($projectdir,$projectname,$libtype,$igblastdir,$packagedir);
	my ($cycle1dir,$cycle3dir,$selectedir,$input,$joined);
	my (@filename,$i,$k,$temp,$listsize,$germclass);
	my ($input1,$input2,$listnt,$numberedf,$germlinedir);
	my ($enriched,$enrichedn,$numbered);
	my (%numberofsequences,$p,$q,$csv);
	my ($filteringout,$pipelinerrorlog);
	my ($configfile,$logfile) = @ARGV;
	my ($input1fasta,$input2fasta);
	my ($pairedend,$minlen,$minqual,$maxlen);
	my ($prinseqcommand1,$prinseqcommand2);
	my ($input1r1dir,$input1r2dir);
	my ($input2r1dir,$input2r2dir);
	my ($afterfiltering1,$afterfiltering2,$beforefilter1,$beforefilter2);
	my ($plot1,$plot2,$workdir,$reportdir, $copy,$temp1,$temp2,$space);
	$minlen = "";
	$maxlen = "";
	$minqual = "";
	my($fastqcommand1,$fastqcommand2);
	my ($rscriptcommand,$sizer0,$sizern,$statscript,$ncandidates);

			
#~---------------------------------------------------------------------------------------------
#~ 				Leitura do arquivo de configuração
#~--------------------------------------------------------------------------------------------
	if(@ARGV < 2){
		die "autoiganalysis: perl autoiganalysis.pl <input.cfg> <error.log>\n";
	}
	open STDERR, ">>$logfile" or die "autoiganalysis.pl: Could not create log file: $!\n";
	open IN, "<$configfile" or die "autoiganalysis.pl: Could not open configuration file: $!\n";
	
	while($line = <IN>){
		chomp($line);
		#~ Passo. Se a linha contiver o nome do arquivo do ciclo 1
		if($line =~ /^input1dir:\s(.+)/){
			$input1dir = $1;	
		}
		elsif($line =~ /^input2dir:\s(.+)/){
			$input2dir = $1;				
		}
		elsif($line =~ /^input1r1dir:\s(.+)/){
			$input1r1dir = $1;
		}
		elsif($line =~ /^input1r2dir:\s(.+)/){
			$input1r2dir = $1;
		}
		elsif($line =~ /^input2r1dir:\s(.+)/){
			$input2r1dir = $1;
		}
		elsif($line =~ /^input2r2dir:\s(.+)/){
			$input2r2dir = $1;
		}
		elsif($line =~ /^projectdir:\s(.+)/){
			$projectdir = $1;		
		}
		elsif($line =~ /^projectname:\s(.+)/){
			$projectname = $1;
		}
		elsif($line =~ /^packagedir:\s(.+)/){
			$packagedir = $1;
		}
		elsif($line =~ /^libtype:\s(.+)/){
			$libtype = $1;				
		}
		elsif($line =~ /^listsize:\s(.+)/){
			$listsize = $1;
		}
		elsif($line =~ /^pairedend:\s(.+)/){
			$pairedend = $1;
		}
		elsif($line =~ /^igblastdir:\s(.+)/){
			$igblastdir = $1;
		}
		elsif($line =~ /^minlen:\s(.+)/){
			$minlen = $1;
		}
		else {
			if($line =~ /^minqual:\s(.+)/){
				$minqual = $1;
			}
		}
	}
	
	close IN;
#~---------------------------------------------------------------------------------------------
#~		 	Criação de diretórios e nomes de arquivos do projeto
#~---------------------------------------------------------------------------------------------
	#~ Passo. Crie o diretório do projeto
	if($libtype == 0){
		$aux = $projectdir . $projectname . "/VH" ;
	}
	else{
		$aux = $projectdir . $projectname . "/VL" ;
	}
	system("mkdir $aux");
	#~ Passo. Crie os subdiretórios do ciclo 1 e do ciclo 3
	$cycle1dir = $aux . "/InitialRound";
	system("mkdir $cycle1dir");
	$cycle3dir = $aux . "/FinalRound";
	system("mkdir $cycle3dir");
	$selectedir = $aux . "/SelectedSequences";
	system("mkdir $selectedir");
	$reportdir = $projectdir . $projectname . "/Report" ;
	$copy = "cp " . $packagedir . "data/*.png " . $reportdir ;
	system("$copy");
	#~ Passo. Obtenha os nomes dos arquivos de entrada dos reads single-end
	if($pairedend == 0){
		if($input1dir =~ /(.+)\/(.+)\.fq$/ || $input1dir =~ /(.+)\/(.+)\.fastq$/){
			$input1 = $2;
		}
		if($input2dir =~ /(.+)\/(.+)\.fq$/ || $input2dir =~ /(.+)\/(.+)\.fastq$/){
			$input2 = $2;
		}
	}	
	#Passo. Obtenha os nomes dos arquivos de entrada dos reads paired-end	
	if($pairedend == 1){
		if($input1r1dir =~ /(.+)\/(.+)\.fq$/ || $input1r1dir =~ /(.+)\/(.+)\.fastq$/){
			$input1 = $2;
		}
		if($input2r1dir =~ /(.+)\/(.+)\.fq$/ || $input2r1dir =~ /(.+)\/(.+)\.fastq$/){
			$input2 = $2;
		}
	}	
	
	#~ Passo. Crie nomes para os arquivos de entrada e saída de todas as etapas do pipeline
	$input1fasta = $cycle1dir . "/" . $input1 ;
	$input2fasta = $cycle3dir . "/" . $input2 ;
	$aux = "";
	$i = 0;
	while($i < 2){
		if($i == 0){
			$k = 0;
			$aux = $input1fasta;
		}
		else{
			$k = 3;
			$aux = $input2fasta;
		}
		
		$filename[$k] = $aux . "aa.fasta";
		$filename[$k+1] = $aux . "nt.fasta";
		$filename[$k+2] = $filename[$k];
		$filename[$k+2] =~ s/.fasta//;
		$filename[$k+2] .= "freq.txt";
		$i++;
	}

	if($libtype == 0){
		$aux = "vh";
	}
	else{
		$aux = "vl";
	}
   	$enriched = $selectedir . "/" . $aux . "enriched" . ".fasta";
	$enrichedn = $selectedir . "/" . $aux . "list" . $listsize . ".fasta";
	$numbered = $selectedir . "/" . $aux . "list" . $listsize . "numbered.txt";
	$numberedf = $selectedir . "/" . $aux . "list" . $listsize . "numbered.fasta";
	$listnt = $selectedir . "/" . $aux . "list" . $listsize . "numberednt.fasta";
	$germclass = $selectedir . "/" . $aux . "list" . $listsize . "numbered" . "germlineclassification.txt" ;
	$pipelinerrorlog = $projectdir . $projectname . "/" . $aux . "PipelineError.log";		
	$gencodedir = $packagedir . "data/genetic_code.fasta";
	$plot1 = $reportdir . "/" . "plot1" . $aux . ".r" ;
	$plot2 = $reportdir . "/" . "plot2" . $aux . ".r" ;
	$statscript = $reportdir . "/" . $aux . "diffproportion.R " ;

#~---------------------------------------------------------------------------------------------
#~ 				Etapas do pipeline
#~---------------------------------------------------------------------------------------------

	# Passo. Entre no diretório de links simbólicos
	chdir("ATTILASymLinks/");
	#~ Passo. Controle de qualidade
	if($pairedend == 0){
		system("nohup fastqc -q -o $cycle1dir $input1dir >> $pipelinerrorlog 2>&1");
		system("nohup fastqc -q -o $cycle3dir $input2dir >> $pipelinerrorlog 2>&1");
	}
	else{
		system("nohup fastqc -q -o $cycle1dir $input1r1dir >> $pipelinerrorlog 2>&1");
		system("nohup fastqc -q -o $cycle1dir $input1r2dir >> $pipelinerrorlog 2>&1");
		system("nohup fastqc -q -o $cycle3dir $input2r1dir >> $pipelinerrorlog 2>&1");
		system("nohup fastqc -q -o $cycle3dir $input2r2dir >> $pipelinerrorlog 2>&1");
		#~ Passo. Montagem de reads
		$fastqcommand1 = "nohup fastq-join " . $input1r1dir . " " . $input1r2dir  . " -o " . $input1fasta . ".fq" . " >> " .  $pipelinerrorlog . " 2>&1";
		$fastqcommand2 = "nohup fastq-join " . $input2r1dir . " " . $input2r2dir  . " -o " . $input2fasta . ".fq" . " >> " .  $pipelinerrorlog . " 2>&1";
		system($fastqcommand1);
		system($fastqcommand2);
		#~ Passo. Atualize os nomes dos arquivos de entrada para a etapa de filtragem
		$input1dir = $input1fasta . ".fqjoin";
		$input2dir = $input2fasta . ".fqjoin";
       }

	# Passo. Modifique os identificadores dos reads caso contenham espaço em branco
	$space = `egrep -c "^@(.+)\\s(.+)" $input1dir` ;
	chomp($space);
	if ($space != 0){
		$temp1 = $input1dir;
		$temp2 = $input2dir;
		$input1dir = "";
		$input2dir = "";
		$input1dir = $input1fasta . "newid.fq" ;
		$input2dir = $input2fasta . "newid.fq" ;
		system("nohup perl parserid.pl $temp1 $input1dir >> $pipelinerrorlog 2>&1");
		system("nohup perl parserid.pl $temp2 $input2dir >> $pipelinerrorlog 2>&1");
	}

	#~ Passo. Filtragem
	if($minqual ne ""){
		$aux = $minqual;
		$minqual = " -min_qual_mean " . $aux;
	}
	if($minlen ne ""){
		$aux = $minlen;
		$minlen = " -min_len " . $aux;
	}
	if($maxlen ne ""){
		$aux = $maxlen;
		$maxlen = " -max_len " . $aux ;
	}

	# Passo. Crie arquivos fasta a partir dos arquivos fastq originais, sem filtar
	$beforefilter1 = $input1fasta . "beforefiltering";
	$beforefilter2 = $input2fasta . "beforefiltering";
	$prinseqcommand1 =  "nohup prinseq-lite -fastq " . $input1dir . " -out_format 1 -out_bad null -out_good " . $beforefilter1 . " >> $pipelinerrorlog 2>&1" ;
	$prinseqcommand2 =  "nohup prinseq-lite -fastq " . $input2dir . " -out_format 1 -out_bad null -out_good " . $beforefilter2 . " >> $pipelinerrorlog 2>&1" ;
	
	system($prinseqcommand2);
	system($prinseqcommand1);
	$prinseqcommand1 = "";
	$prinseqcommand2 = "";

	$afterfiltering1 = $input1fasta . "filtered";
	$prinseqcommand1 = "nohup prinseq-lite -fastq " . $input1dir . $minqual . $minlen . $maxlen . " -out_format 5 -out_bad null -out_good " . $afterfiltering1 . " >> $pipelinerrorlog 2>&1" ;
	$afterfiltering2 = $input2fasta . "filtered" ;
	$prinseqcommand2 = "nohup prinseq-lite -fastq " . $input2dir . $minqual . $minlen . $maxlen . " -out_format 5 -out_bad null -out_good " . $afterfiltering2 . " >> $pipelinerrorlog 2>&1";
	
	system($prinseqcommand1);
	system($prinseqcommand2);
	$afterfiltering1 .= ".fastq";
	$afterfiltering2 .= ".fastq";

	# Passo. Controle de qualidade após filtragem
	system("nohup fastqc -q -o $cycle1dir $afterfiltering1 >> $pipelinerrorlog 2>&1");
	system("nohup fastqc -q -o $cycle3dir $afterfiltering2 >> $pipelinerrorlog 2>&1");

	$afterfiltering1 =~ s/\.fastq/\.fasta/;
	$afterfiltering2 =~ s/\.fastq/\.fasta/;

	#~ Passo. Tradução
	system("nohup ./translateab9 $afterfiltering1 $gencodedir $filename[0] $filename[1] $libtype >> $pipelinerrorlog 2>&1");
	system("nohup ./translateab9 $afterfiltering2 $gencodedir $filename[3] $filename[4] $libtype >> $pipelinerrorlog 2>&1");

	#~ Passo. Cálculo de frequência relativa
	system("nohup perl frequency_counter3.pl $filename[0] $afterfiltering1 $filename[2] >> $pipelinerrorlog 2>&1");
	system("nohup perl frequency_counter3.pl $filename[3] $afterfiltering2 $filename[5] >> $pipelinerrorlog 2>&1");

	#~ Passo. Cálculo de fold change
	system("nohup perl find_duplicates7.pl $filename[2] $filename[5] $enriched >> $pipelinerrorlog 2>&1");
	system("nohup perl get_nsequences.pl $enriched $enrichedn $listsize >> $pipelinerrorlog 2>&1");
	
	#~ Passo. Numeração de resíduos
	system("nohup perl numberab.pl $enrichedn $numbered >> $pipelinerrorlog 2>&1");
	system("nohup perl convertofasta.pl $numbered $numberedf >> $pipelinerrorlog 2>&1");
	
	#~ Passo. Recuperação de sequências de nucleotídeos
	system("nohup perl get_ntsequence2.pl $numberedf $enrichedn $filename[4] $listnt >> $pipelinerrorlog 2>&1");
	

	#~ Passo. Teste estatístico de diferença de proporção, usando correção de Bonferroni
	$sizer0 = `grep -c -P "^>" $afterfiltering1`;
	chomp($sizer0);
	$sizern = `grep -c -P "^>" $afterfiltering2`;
	chomp($sizern);
	$ncandidates = `grep -c -P "^>" $numberedf`;
	chomp($ncandidates);
	system("nohup perl statscript_creator.pl $enrichedn $numberedf $reportdir $statscript $sizer0 $sizern $libtype $ncandidates >> $pipelinerrorlog 2>&1");
	system("chmod 777 $statscript");
	system("$statscript");


	# Passo. Armazene o diretório atual
	$workdir = `pwd` ;
	chomp($workdir);
	
	# Passo. Classificação de germlines
	chdir($igblastdir);
	if($libtype == 0){
		$germlinedir = $packagedir . "data/germline_human/VH_germline_human";
	}
	else{
		$germlinedir = $packagedir . "data/germline_human/VL_germline_human";
	}
	system("bin/igblastp -germline_db_V $germlinedir -query $numberedf -organism human -domain_system kabat -outfmt 7 > $germclass");
	
	chdir($workdir) ;
#~---------------------------------------------------------------------------------------------
#~ 				Processamento dos dados produzidos e formatação da saída
#~---------------------------------------------------------------------------------------------
	
	#~ Passo. Contagem de sequências das bibliotecas de cada etapa do pipeline
	$numberofsequences{raw}{number} = 0;
	if($pairedend == 1){
		$numberofsequences{joining}{number} = 1;
	}
	$numberofsequences{filtering}{number} = 3;
	$numberofsequences{translation}{number} = 4;
	$numberofsequences{frequency}{number} = 5;
	$numberofsequences{enrichment}{number} = 6;
	$numberofsequences{numeration}{number} = 7;

	$i = 0;
	$aux = "";
	while($i < 2){
		if($i == 0){
			$filteringout = $afterfiltering1;
			$k = 'input1';
			$p = $i + 1;
			$q = $p + 1;
			if($pairedend == 0){
				$input = $input1dir;
			}
			else{
				$joined = $input1dir;
				$input = $input1r1dir;
			}
		}
		else{
			$filteringout = $afterfiltering2;
			$k = 'input2';
			$p = $i + 3;
			$q = $p + 1;
			if($pairedend == 0){
				$input = $input2dir;
			}
			else{
				$joined = $input2dir;
				$input = $input2r1dir;
			}
		}
		$aux = `wc -l $input | grep -woP "(\\d*)"`;
		chomp($aux);
		$aux = $aux / 4;
		$numberofsequences{raw}{$k} = $aux;
		if($pairedend == 1){
			$aux = `wc -l $joined | grep -woP "(\\d*)"`;
			chomp($aux);
			$aux = $aux / 4;
			$numberofsequences{joining}{$k} = $aux;
		}

		$aux = `grep -c -P "^>" $filteringout`;
		chomp($aux);
		$numberofsequences{filtering}{$k} = $aux;
		$aux = `grep -c -P "^>" $filename[$p]`;
		chomp($aux);
		$numberofsequences{translation}{$k} = $aux;
		$aux = `grep -c -P "^#" $filename[$q]`;
		chomp($aux);
		$numberofsequences{frequency}{$k} = $aux;
		$i++;
	}	
	$aux = `grep -c -P "^>" $enriched`;
	chomp($aux);
	$numberofsequences{enrichment}{sequences} = $aux;
	$aux = `grep -c -P "^>" $numbered`; 
	chomp($aux);
	$numberofsequences{numeration}{sequences} = $aux;
	$aux = "";
	if ($libtype == 0){
		$csv = $projectdir . $projectname . "/VH/vhSequenceCounting.csv";
	}
	else{
		$csv = $projectdir . $projectname . "/VL/vlSequenceCounting.csv";
	}

	#~ Passo. Impressão dos tamanhos das bibliotecas em um arquivo csv
	open OUT, ">$csv" or die "autoiganalysis.pl: Could not create csv file: $!\n";
	print OUT "library,reads,step\n";
	foreach $aux(sort{$numberofsequences{$a}{number} <=> $numberofsequences{$b}{number}} keys %numberofsequences){
		if($aux ne 'enrichment' && $aux ne 'numeration' && $aux ne 'joining'){	
			print OUT "R0,$numberofsequences{$aux}{input1},$aux\n";
			print OUT "RN,$numberofsequences{$aux}{input2},$aux\n";
		}
		else{
			if($aux eq 'enrichment' or $aux eq 'numeration'){
				print OUT "Selected,$numberofsequences{$aux}{sequences},$aux\n";
			}
			else{
				if($pairedend == 1){
					print OUT "R0,$numberofsequences{$aux}{input1},$aux\n";
					print OUT "RN,$numberofsequences{$aux}{input2},$aux\n";
				}
			}
		}
	}
	close OUT;
	$aux = "";
	if ($libtype == 0){
		$aux = "vh";
	}
	else{
		$aux = "vl";
	}
	$minlen =~ s/-min_len//g ;
 $rscriptcommand = "nohup perl rscript_creator.pl " . $beforefilter1 . ".fasta " . $beforefilter2 . ".fasta " . $csv .  " " . $plot1 . " " . $plot2 . " " . $reportdir . "/" . $minlen . " " . $aux ;
	$aux = "";
	$aux = " >> " .  $pipelinerrorlog . " 2>&1";
	$rscriptcommand .= $aux;
	system($rscriptcommand);
	system("chmod 777 $plot1");
	system("chmod 777 $plot2"); 
	system("$plot1");		
	system("$plot2");
	$aux = "";
	$aux = $beforefilter1 . ".fasta" ;
	system("rm $aux");
	$aux = "";
	$aux = $beforefilter2 . ".fasta" ;
	system("rm $aux");

	close STDERR;
	exit(0);
	
