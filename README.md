# ATTILA - AutomaTed Tool For Immunoglobulin Analysis

ATTILA searches for candidate immunoglobulin sequences in phage display libraries, generating as main output a list of sequences of heavy and light chain, which were selected by phage display experiment, and code for antibody fragments that can probably bind to the target molecule. ATTILA package has programs developed in C, Perl and Shell script to execute eight steps of a completely automated analysis.

ATTILA can analyse human sequences coding for variable domain of heavy (VH) and light chain (VL), produced by phage display technology. Therefore, the input for the method must be VH and VL libraries, from initial and final rounds. Considering that our approach uses distances between canonical aminoacid residues of variable domain based on human sequences, the analysis may also be performed on mouse sequences, since their distances are similar to those of human.


## Requirements

In order to install and run ATTILA, you must have:

Linux system

Internet

FASTQC http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Prinseq-lite http://prinseq.sourceforge.net/

FastqJoin https://code.google.com/archive/p/ea-utils/

Perl https://www.perl.org/get.html

IgBlast ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/

R package https://cran.r-project.org/

ggplot2 Use R function install.packages(“ggplot2”)

scales Use R function install.packages(“scales”)


## Installation

After installing all requirements, perform the following steps:

1. Download ATTILA package at ...

2. Uncompress the tar.gz file using command line:

$ tar -vzxf attila-1.0.tar.gz

3. Go to the directory where you want to install ATTILA, using cd command. Note that ATTILA and IgBlast packages must be subdirectories of the installation directory.

4. Type the following command line:

$ ln -s <path to check_requirements.sh> check_requirements.sh
  
Example: ln -s home/Attila/programs/check_requirements.sh check_requirements.sh

5. Run the following command line:

$ ./check_requirements.sh

If check requirements.sh prints "Type ./attilacli.sh to run ATTILA", then ATTILA is ready to run ! 
If not, check requirements.sh prints a list of requirements that still need to be installed.

