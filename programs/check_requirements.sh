#! /bin/bash


#~ Script: check_requirements.sh
#~ Este script verifica os requisitos para rodar o método, e informa ao 
#~ usuário se ele pode rodar a análise ou se precisa atender a algum dos requisitos ainda.
#~ Este script assume que as bibliotecas ggplot2 e scales estão no diretório default de bibliotecas do R, isto é,
#~ no primeiro diretório listado pela função libPaths() do R. Caso todos os requisitos estejam instalados, este script
# cria um link simbólico para o attilacli2.sh

	echo "Checking requirements ..."
	flag=0
	check=$(which fastqc)
	if [ -z "$check" ]
		then
			flag=`expr $flag + 1`
			install[1]="FASTQC"
	fi
	
	check=$(which prinseq-lite)
	if [ -z "$check" ]
		then
			flag=`expr $flag + 1`
			install[2]="Prinseq-lite"
	fi
	
	check=$(which fastq-join)
	if [ -z "$check" ]
		then
			flag=`expr $flag + 1`
			install[3]="FastqJoin"
	fi
	
	check=$(which perl)
	if [ -z "$check" ]
		then
			flag=`expr $flag + 1`
			install[4]="Perl"
	fi

	check=$(which R)
	if [ -z "$check" ]
		then
			flag=`expr $flag + 1`
			install[5]="R"
		else
			r=$(find `pwd` -name "findrlib.r")
			dir=$(Rscript $r | grep "[1]" | cut -d '"' -f2)
	fi

	check=$(find $dir -name "ggplot2")
	if [ -z "$check" ]
		then
			flag=`expr $flag + 1`
			install[6]="ggplot2"
	fi
	
	check=$(find $dir -name "scales")
	if [ -z "$check" ]
		then
			flag=`expr $flag + 1`
			install[7]="scales"
	fi
	
	check=$(find `pwd` -name "igblastp")
	if [ -z "$check" ]
		then
			flag=`expr $flag + 1`
			install[8]="IgBlast"
		else
			echo "$check" >> paths_attila.txt
	fi
	
	
	check=$(find `pwd` -name "attilacli.sh")
	echo "$check" >> paths_attila.txt
	
	if [ $flag -ne 0 ]
	then
		echo "Please install the following tools/packages"
		for i in "${install[@]}"
		do
			if [ $i == "FASTQC" ]
			then 
				echo "$i: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/ "
			elif [ $i == "Prinseq-lite" ]
			then
				echo "$i: http://prinseq.sourceforge.net/ "
			elif [ $i == "FastqJoin" ]
			then
				echo "$i: https://code.google.com/archive/p/ea-utils/ "
			elif [ $i == "Perl" ]
			then
				echo "$i: https://www.perl.org/get.html "
			elif [ $i == "IgBlast" ]
			then
				echo "$i: ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/ "
			elif [ $i == "R" ]
			then
				echo "$i: https://cran.r-project.org/ "
				echo "$i: Use R function 'install.packages(\"ggplot2\") "
				echo "$i: Use R function 'install.packages(\"scales\") "
			elif [ $i == "ggplot2" ]
			then
				echo "$i: Use R function 'install.packages(\"ggplot2\") "
			else
				echo "$i: Use R function 'install.packages(\"scales\") "
			fi
		done
	else
		ln -s $check ./attilacli.sh
		echo "Type "./attilacli.sh" to run ATTILA "
	fi
	
