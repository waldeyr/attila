#! /bin/bash

# ---------------------------------------------------------------------------
# Script: iganalyzer.sh
# Este script lê informações dadas pelo usuário para definir
# as configurações da automatização da análise de 
# sequências de imunoglobulinas, desenvolvida pelo grupo de 
# Bioinformática da UnB. O programa imprime ajuda caso o usuário
# entre com o argumento --help. Após imprimir as configurações num
# arquivo, são criados links simbólicos para todos os programas 
# pertencentes ao pacote IgAnalyzer, no diretório atual. Finalmente,
# este script shell executa o script perl de automatização da análise.
# ----------------------------------------------------------------------------
	echo "***********************************************************************************************************************************"
	echo "  				IgAnalyzer: A tool for immunoglobulin sequence analysis"
	echo "***********************************************************************************************************************************"
	echo " "
# --------------------------------------------------------------------------------------------
					# Analysis settings
# --------------------------------------------------------------------------------------------
	# settings[1]			Project name
	# settings[2]			Project path
	# settings[3]			IgAnalyzer package path
	# settings[4]			Reads are paired-end (0/1)
	# settings[5]			VH R0 reads r1 path
	# settings[6]			VH R0 reads r2 path
	# settings[7]			VH Rs reads r1 path
	# settings[8]			VH Rs reads r2 path
	# settings[9]			VL R0 reads r1 path
	# settings[10]			VL R0 reads r2 path
	# settings[11]			VL Rs reads r1 path
	# settings[12]			VL Rs reads r2 path
	# settings[13]			VH R0 path
	# settings[14]			VH Rs path
	# settings[15]			VL R0 path
	# settings[16]			VL Rs path
	# settings[17]			IgBlast package path
	# settings[18]			Minimum read length 
	# settings[19]			Minimum base quality
	# settings[20]			Number of candidates to rank

	echo -n "Settings file exists (y or n): "
	read "config"
	if [[ $config == "" || $config == "n" ]]
	then
		echo -n "Enter project name: "
		read -e "settings[1]"
		flag=0
		while [ $flag -eq 0 ]
		do		
			echo -n "Enter project path: "
			read -e "settings[2]"
			if [ ! -d ${settings[2]} ]
			then
				echo "Error: directory does not exist or is not a directory"
			else
				flag=1
			fi		
		done
		
		flag=0
		while [ $flag -eq 0 ]
		do
			echo -n "Enter IgAnalyzer package path: "
			read -e "settings[3]"	
			if [ ! -d ${settings[3]} ]
			then
				echo "Directory does not exist or is not a directory"
			else
				flag=1
			fi
		done
	
		echo -n "Reads are paired-end (y or n): "
		read  "reads"
		if [[ $reads == "y"  || $reads == "" ]]
		then
			settings[4]=1
		else
			settings[4]=0
		fi
		if [[ ${settings[4]} =~ ^[0,1]$ ]]
		then
			if [ ${settings[4]} -eq 1 ]
			then
				flag=0
				while [ $flag -eq 0 ]
				do
					echo -n "Enter VH R0 reads r1 path: "
					read -e "settings[5]"
					if [ ! -f ${settings[5]} ]
					then
						echo "Error: file does not exist or is not a regular file"
					else
						if [[ ${settings[5]} =~ (.+)\.fastq$ ]]
						then
							flag=1
						fi
						if [[ ${settings[5]} =~ (.+)\.fq$ ]]
						then
							flag=1
						fi
						if [ $flag -eq 0 ]
						then
							echo "Input format is incorrect. Please use fastq format"
						fi
					fi
				done
				flag=0
				while [ $flag -eq 0 ]
				do
					echo -n "Enter VH R0 reads r2 path: "
					read -e "settings[6]"
					if [ ! -f ${settings[6]} ]
					then
						echo "Error: file does not exist or is not a regular file"
					else
						if [[ ${settings[6]} =~ (.+)\.fastq$ ]]
						then
							flag=1
						fi
						if [[ ${settings[6]} =~ (.+)\.fq$ ]]
						then
							flag=1
						fi
						if [ $flag -eq 0 ]
						then
							echo "Input format is incorrect. Please use fastq format"
						fi
					fi
				done
				flag=0
				while [ $flag -eq 0 ]
				do
					echo -n "Enter VH Rs reads r1 path: "
					read -e "settings[7]"
					if [ ! -f ${settings[7]} ]
					then
						echo "Error: file does not exist or is not a regular file"
					else
						if [[ ${settings[7]} =~ (.+)\.fastq$ ]]
						then
							flag=1
						fi
						if [[ ${settings[7]} =~ (.+)\.fq$ ]]
						then
							flag=1
						fi
						if [ $flag -eq 0 ]
						then
							echo "Input format is incorrect. Please use fastq format"
						fi
					fi
				done
				flag=0
				while [ $flag -eq 0 ]
				do
					echo -n "Enter VH Rs reads r2 path: "
					read -e "settings[8]"
					if [ ! -f ${settings[8]} ]
					then
						echo "Error: file does not exist or is not a regular file"
					else
						if [[ ${settings[8]} =~ (.+)\.fastq$ ]]
						then
							flag=1
						fi
						if [[ ${settings[8]} =~ (.+)\.fq$ ]]
						then
							flag=1
						fi
						if [ $flag -eq 0 ]
						then
							echo "Input format is incorrect. Please use fastq format"
						fi
					fi
				done
				flag=0
				while [ $flag -eq 0 ]
				do
					echo -n "Enter VL R0 reads r1 path: "
					read -e "settings[9]"
					if [ ! -f ${settings[9]} ]
					then
						echo "Error: file does not exist or is not a regular file"
					else
						if [[ ${settings[9]} =~ (.+)\.fastq$ ]]
						then
							flag=1
						fi
						if [[ ${settings[9]} =~ (.+)\.fq$ ]]
						then
							flag=1
						fi
						if [ $flag -eq 0 ]
						then
							echo "Input format is incorrect. Please use fastq format"
						fi
					fi
				done
				flag=0
				while [ $flag -eq 0 ]
				do
					echo -n "Enter VL R0 reads r2 path: "
					read -e "settings[10]"
					if [ ! -f ${settings[10]} ]
					then
						echo "Error: file does not exist or is not a regular file"
					else
						if [[ ${settings[10]} =~ (.+)\.fastq$ ]]
						then
							flag=1
						fi
						if [[ ${settings[10]} =~ (.+)\.fq$ ]]
						then
							flag=1
						fi
						if [ $flag -eq 0 ]
						then
							echo "Input format is incorrect. Please use fastq format"
						fi
					fi
				done
				flag=0
				while [ $flag -eq 0 ]
				do
					echo -n "Enter VL Rs reads r1 path: "
					read -e "settings[11]"
					if [ ! -f ${settings[11]} ]
					then
						echo "Error: file does not exist or is not a regular file"
					else
						if [[ ${settings[11]} =~ (.+)\.fastq$ ]]
						then
							flag=1
						fi
						if [[ ${settings[11]} =~ (.+)\.fq$ ]]
						then
							flag=1
						fi
						if [ $flag -eq 0 ]
						then
							echo "Input format is incorrect. Please use fastq format"
						fi	
					fi
				done
				flag=0
				while [ $flag -eq 0 ]
				do
					echo -n "Enter VL Rs reads r2 path: "
					read -e "settings[12]"
					if [ ! -f ${settings[12]} ]
					then
						echo "Error: file does not exist or is not a regular file"
					else
						if [[ ${settings[12]} =~ (.+)\.fastq$ ]]
						then
							flag=1
						fi
						if [[ ${settings[12]} =~ (.+)\.fq$ ]]
						then
							flag=1
						fi
						if [ $flag -eq 0 ]
						then
							echo "Input format is incorrect. Please use fastq format"
						fi
					fi
				done
				flag=0
				while [ $flag -eq 0 ]
				do
					echo -n "Enter VH R0 path: "
					read -e "settings[13]"
					if [ ! -f ${settings[13]} ]
					then
						echo "Error: file does not exist or is not a regular file"
					else
						if [[ ${settings[13]} =~ (.+)\.fastq$ ]]
						then
							flag=1
						fi
						if [[ ${settings[13]} =~ (.+)\.fq$ ]]
						then
							flag=1
						fi
						if [ $flag -eq 0 ]
						then
							echo "Input format is incorrect. Please use fastq format"
						fi
					fi
				done
				flag=0
				while [ $flag -eq 0 ]
				do
					echo -n "Enter VH Rs path: "
					read -e "settings[14]"
					if [ ! -f ${settings[14]} ]
					then
						echo "Error: file does not exist or is not a regular file"
					else
						if [[ ${settings[14]} =~ (.+)\.fastq$ ]]
						then
							flag=1
						fi
						if [[ ${settings[14]} =~ (.+)\.fq$ ]]
						then
							flag=1
						fi
						if [ $flag -eq 0 ]
						then
							echo "Input format is incorrect. Please use fastq format"
						fi	
					fi
				done
				flag=0
				while [ $flag -eq 0 ]
				do
					echo -n "Enter VL R0 path: "
					read -e "settings[15]"
					if [ ! -f ${settings[15]} ]
					then
						echo "Error: file does not exist or is not a regular file"
					else
						if [[ ${settings[15]} =~ (.+)\.fastq$ ]]
						then
							flag=1
						fi
						if [[ ${settings[15]} =~ (.+)\.fq$ ]]
						then
							flag=1
						fi
						if [ $flag -eq 0 ]
						then
							echo "Input format is incorrect. Please use fastq format"
						fi
					fi
				done
				flag=0
				while [ $flag -eq 0 ]
				do		
					echo -n "Enter VL Rs path: "
					read -e "settings[16]"
					if [ ! -f ${settings[16]} ]
					then
						echo "Error: file does not exist or is not a regular file"
					else
						if [[ ${settings[16]} =~ (.+)\.fastq$ ]]
						then
							flag=1
						fi
						if [[ ${settings[16]} =~ (.+)\.fq$ ]]
						then
							flag=1
						fi
						if [ $flag -eq 0 ]
						then
							echo "Input format is incorrect. Please use fastq format"
						fi
					fi
				done
			fi
		fi
		
		flag=0
		while [ $flag -eq 0 ]
		do
			echo -n "Enter IgBlast package path: "
			read -e "settings[17]"
			if [ ! -d ${settings[17]} ]
			then
				echo "Directory does not exist or is not a directory"
			else
				flag=1
			fi
		done
		
		echo -n "Minimum read length (default=300) Change (y or n): "
		read "chminlen"
		if [[ $chminlen == "" || $chminlen == "n" ]]
		then
			settings[18]=300
		else
			flag=0
			while [ $flag -eq 0 ]
			do
				echo -n "Enter minimum read length: "
				read "setting[18]"
				if [[ ! ${settings[18]} =~ ^[0-9]+$ ]]
				then
					echo "Invalid value"
				else
					flag=1
				fi
			done
		fi
	
		echo -n "Minimum base quality (default=20) Change (y or n): "
		read "chminqual"
		if [[ $chminqual == "" || $chminqual == "n" ]]
		then
			settings[19]=20
		else
			flag=0
			while [ $flag -eq 0 ]
			do
				echo -n "Enter minimum base quality: "
				read "settings[19]"
				if [[ ! ${settings[19]} =~ ^[0-9]+$ ]]
				then
					echo "Invalid value"
				else
					flag=1
				fi
			done
		fi
	
		echo -n "Enter number of candidates to rank: "
		flag=0
		while [ $flag -eq 0 ]
		do
			read "settings[20]"
			if [[ ! ${settings[20]} =~ ^[0-9]+$ ]]
			then
				echo "Invalid value"
			else
				flag=1
			fi
		done

		else
		flag=0
		while [ $flag -eq 0 ]
		do	
			echo -n "Enter VH setting file path: "
			read -e "vhfilecfg"
			if [ ! -f $vhfilecfg ]
			then
				echo "File does not exist or is not a regular file"
			else
				format1=`egrep "*.fastq" $vhfilecfg`
				format2=`egrep "*.fq" $vhfilecfg`
				if [[ ${format1} =~ (.+)\.fastq$ ]]
				then
					flag=1
				fi
				if [[ ${format2} =~ (.+)\.fq$ ]]
				then
					flag=1
				fi
				if [ $flag -eq 0 ]
				then
					echo "Input format is incorrect. Please use fastq format."
				fi
			fi
		done
		flag=0
		while [ $flag -eq 0 ]
		do
			echo -n "Enter VL setting file path: "
			read -e "vlfilecfg"
			if [ ! -f $vlfilecfg ]
			then
				echo "File does not exist or is not a regular file"
			else
				format3=`egrep "*.fastq" $vlfilecfg`
				format4=`egrep "*.fq" $vlfilecfg`
				if [[ ${format3} =~ (.+)\.fastq$ ]]
				then
					flag=1
				fi
				if [[ ${format4} =~ (.+)\.fq$ ]]
				then
					flag=1
				fi
				if [ $flag -eq 0 ]
				then
					echo "Input format is incorrect. Please use fastq format."
				fi
			fi
		done
		settings[1]=`grep "projectname:" $vhfilecfg | sed 's/projectname: //g'`
		settings[2]=`grep "projectdir:" $vhfilecfg | sed 's/projectdir: //g'`
		settings[3]=`grep "packagedir:" $vhfilecfg | sed 's/packagedir: //g'`
	fi
	
#--------------------------------------------------------------------------------------------
				# Check settings
# -------------------------------------------------------------------------------------------
	if [ $config == "n" ]
	then
		echo "------------------------------------------------------------------------------------------------------------------------------------"
		echo "---------------------------------------------Settings for current analysis----------------------------------------------------------"
		echo " "

		echo "Project name: ${settings[1]}"
		echo "Project path: ${settings[2]}"
		echo "IgAnalyzer path: ${settings[3]}"
		if [ ${settings[4]} -eq 1 ]
		then
			echo "Reads are paired-end: yes"
			echo "VH R0 reads r1: ${settings[5]}"
			echo "VH R0 reads r2: ${settings[6]}"
			echo "VH Rs reads r1: ${settings[7]}"
			echo "VH Rs reads r2: ${settings[8]}"
			echo "VL R0 reads r1: ${settings[9]}"
			echo "VL R0 reads r2: ${settings[10]}"
			echo "VL Rs reads r1: ${settings[11]}"
			echo "VL Rs reads r2: ${settings[12]}"
		else
			echo "Reads are paired-end: no"
			echo "VH R0: ${settings[13]}"
			echo "VH Rs: ${settings[14]}"
			echo "VL R0: ${settings[15]}"
			echo "VL Rs: ${settings[16]}"
		fi
		echo "IgBlast package path: ${settings[17]}"
		echo "Minimum read length: ${settings[18]}"
		echo "Mininum base quality: ${settings[19]}"
		echo "Number of candidates: ${settings[20]}"
		echo "------------------------------------------------------------------------------------------------------------------------------------"
		echo "Configuration is correct (y or n): "
		c=21
		read "start"
		if [[ $start != "" && $start != "y" ]] 
		then
			echo "--------------------------------Menu to correct settings----------------------------------------------------------------------------"
			echo " "
			echo "Project Name (1)			Project Path (2)		IgAnalyzer Package Path (3)	    Paired-end Reads (4)"
			if [ ${settings[4]} -eq 1 ]
			then
				echo "VH R0 reads r1 (5)	         	VH R0 reads r2 (6)		VH Rs reads r1 (7)		    VH Rs reads r2 (8)"
				echo "VL R0 reads r1 (9)	         	VL R0 reads r2 (10)		VL Rs reads r1 (11)		    VL Rs reads r2 (12)"
			else
				echo "VH R0 (13)				VH Rs (14)			VL R0 (15)		            VL Rs (16)"
			fi
			echo "IgBlast Package Path (17)		Minimum Read Length (18)	Mininum Base Quality (19)"
			echo "Number of Candidates (20)		Save and exit settings (0)"
			echo "------------------------------------------------------------------------------------------------------------------------------------"
			while [ $c -ne 0 ] 	
			do 
				echo -n "Enter corresponding integer to correct settings: "
				read "c"
				if [ $c -ne 0 ]
				then
					echo -n "Enter new setting: "
					read -e "set"
					settings[$c]=$set
				fi
			done
		fi
	fi

# -------------------------------------------------------------------------------------------
		# Write settings to file
# -------------------------------------------------------------------------------------------
	if [[ $config == "" || $config == "n" ]]
	then
		for i in 0 1
		do
			if [ $i -eq 0 ]
			then
				filecfg="${settings[1]}set_iganalysisvh.cfg"
				vhfilecfg=$filecfg
			else
				filecfg="${settings[1]}set_iganalysisvl.cfg"
				vlfilecfg=$filecfg
			fi
			echo "----------------------------------------------------------------------" >> $filecfg
			echo "# Settings for immunoglobulin sequence analysis" >> $filecfg
			echo "# [ Section: files and directories ]" >> $filecfg
			echo "projectname: ${settings[1]}" >> $filecfg
			echo "projectdir: ${settings[2]}"	>> $filecfg
			echo "packagedir: ${settings[3]}" >> $filecfg	
			echo "igblastdir: ${settings[17]}" >> $filecfg
			if [ $i -eq 0 -a ${settings[4]} -eq 0 ]
			then
				echo "input1dir: ${settings[13]}" >> $filecfg
				echo "input2dir: ${settings[14]}" >> $filecfg
			elif [ $i -eq 0 -a ${settings[4]} -eq 1 ]
			then
				echo "input1r1dir: ${settings[5]}" >> $filecfg
				echo "input1r2dir: ${settings[6]}" >> $filecfg
				echo "input2r1dir: ${settings[7]}" >> $filecfg
				echo "input2r2dir: ${settings[8]}" >> $filecfg
			elif [ $i -eq 1 -a ${settings[4]} -eq 0 ] 	
			then
				echo "input1dir: ${settings[15]}" >> $filecfg
				echo "input2dir: ${settings[16]}" >> $filecfg
			else
				echo "input1r1dir: ${settings[9]}" >> $filecfg
		                echo "input1r2dir: ${settings[10]}" >> $filecfg
			        echo "input2r1dir: ${settings[11]}" >> $filecfg
				echo "input2r2dir: ${settings[12]}" >> $filecfg
			fi 
			echo "# [ Section: analysis arguments ]" >> $filecfg
			echo "libtype: $i" >> $filecfg
			echo "listsize: ${settings[20]}"	>> $filecfg
			echo "pairedend: ${settings[4]}" >> $filecfg
			echo "minlen: ${settings[18]}" >> $filecfg
			echo "minqual: ${settings[19]}" >> $filecfg
		done
	fi
# --------------------------------------------------------------------------------------------
		# Create symbolic links
# --------------------------------------------------------------------------------------------

	echo -n "Create symbolic links (y or n): "
	read "links"
	if [[ $links == "" || $links == "y" ]]
	then
		autoiganalysis="${settings[3]}programs/autoiganalysis3.pl"
		translateab9="${settings[3]}programs/translateab9"
		frequencycounter="${settings[3]}programs/frequency_counter3.pl"
		finduplicates="${settings[3]}programs/find_duplicates7.pl"
		getnsequences="${settings[3]}programs/get_nsequences.pl"
		numberab="${settings[3]}programs/numberab.pl"
		convertofasta="${settings[3]}programs/convertofasta.pl"
		getntsequence="${settings[3]}programs/get_ntsequence2.pl"
		rscriptcreator="${settings[3]}programs/rscript_creator.pl"
		htmlcreator="${settings[3]}programs/html_creator.pl"
		parserid="${settings[3]}programs/parserid.pl"



		ln -s $autoiganalysis autoiganalysis3.pl
		ln -s $translateab9 translateab9
		ln -s $frequencycounter frequency_counter3.pl
		ln -s $finduplicates find_duplicates7.pl
		ln -s $getnsequences get_nsequences.pl
		ln -s $getntsequence get_ntsequence2.pl
		ln -s $numberab numberab.pl
		ln -s $convertofasta convertofasta.pl	
		ln -s $rscriptcreator rscript_creator.pl
		ln -s $htmlcreator html_creator.pl
		ln -s $parserid parserid.pl
	fi
# --------------------------------------------------------------------------------------------
		# Run analysis
# --------------------------------------------------------------------------------------------
	
	echo "Creating project directory"
	project="${settings[2]}${settings[1]}"
	mkdir $project 
	reportdir="${project}/Report"
	mkdir $reportdir
	echo "Running VH analysis ..."
	vherrorlog="${project}/vherror.log"
	time perl autoiganalysis3.pl $vhfilecfg $vherrorlog
	echo "---------------------------------------------------------------------------------"
	echo "VH Analysis Completed"
	echo "---------------------------------------------------------------------------------"
	vlerrorlog="${project}/vlerror.log"
	echo "Running VL analysis ..."
	time perl autoiganalysis3.pl $vlfilecfg $vlerrorlog
	echo "---------------------------------------------------------------------------------"
	echo "VL Analysis Completed"
	echo "---------------------------------------------------------------------------------"
	
# --------------------------------------------------------------------------------------------------
		# Create web page
# --------------------------------------------------------------------------------------------------
	oldir=`pwd`
	selectedirvh="${project}/VH/SelectedSequences/"
	cd $selectedirvh
	numberedvh=`ls | egrep "*numbered.fasta"`
	germlinevh=`ls | egrep "*germlineclassification.txt"`
	cd $reportdir
	plot1vh=`ls | egrep "*length_vh.png"`
	plot2vh=`ls | egrep "task_vh.png"`
	selectedirvl="${project}/VL/SelectedSequences/"
	cd $selectedirvl
	numberedvl=`ls | egrep "*numbered.fasta"`
	germlinevl=`ls | egrep "*germlineclassification.txt"`
	cd $reportdir
	plot1vl=`ls | egrep "*length_vl.png"`
	plot2vl=`ls | egrep "task_vl.png"`
	cd $oldir
	
	perl html_creator.pl "${project}/VH/SelectedSequences/${numberedvh}" "${project}/VL/SelectedSequences/${numberedvl}" "${project}/VH/SelectedSequences/${germlinevh}" "${project}/VL/SelectedSequences/${germlinevl}" "${reportdir}/Report.html" "$plot1vh" "$plot2vh" "$plot1vl" "$plot2vl" > "${reportdir}/webpage.log" 2>&1
	
	webpagelog=`du -b "${reportdir}/webpage.log" | cut -f1`
	if [ $webpagelog == "0" ]
	then
		echo "Web page is ready !"
	else
		echo "Could not create web page"
	fi
	exit 0
