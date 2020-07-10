# Processing MS raw data <!-- omit in toc -->

This document describes the initial processing of the raw mass spectrometry data in order to obtain peptide identification results. There are many ways you can go about this and the described method is just one of them. I like this pipeline because it is sensitive, accurate, and open-source. This pipeline is also based on command line execution from a Linux terminal. However, for all of the tools described there are GUI implementations available on Windows, so the analysis can easily be ported over to this OS. This document is more of a high-level explanation, so for very detailed explanations of different aspects I would refer you to the individual pages linked below for the specific tools.

## System and software information

These data were processed on a CentOS Linux release 7.4.1708 system with dual Intel(R) Xeon(R) CPU E5-2690 processors (16 total cores, 32 total threads) and 125GB of ram. For software, this script makes use of:

* RawTools [(version 2.0.2)](https://github.com/kevinkovalchik/RawTools)
* SearchGUI [(version 3.3.16)](https://github.com/compomics/searchgui)
* PeptideShaker [(version 1.16.42)](https://github.com/compomics/peptide-shaker)
* R language [(version 3.6.3)](https://www.r-project.org/)
* R package 'Biostrings' (I used version 2.54.0)
* R package 'tidyverse' (I used version 1.3.0)
* R package 'Peptides' (I used version 2.4.2)
* R package 'OrgMassSpecR' (I used version 0.5-3)

All R packages were installed using BiocManager::install('library').

## Running the analysis

The analysis is already coded into two shell scripts, provided below. Although the below text can seem complex, it is really just changing file paths based on your system. In addition, you should only need to do it once. In your linux terminal, navigate to the directory where you would like to store your scripts and prepare to run the main shell script.

```bash
cd /directory/where/you/want/to/store/scripts
vim nameOfScript.sh
```

Paste the code for the shell script provided below in the opened vim file and edit as required (press 'i' while vim is open to be able to change text, hit Esc when you are done editing, then type ':wq' and hit Enter to save and exit). There are a couple of spots where you will need to edit the file, mostly to change the file paths for your own system or change the modifications for your own analysis. For example, you will need to change:

* sampleTag="10Oct2016_ADROIT_PrepCompare_TMT10_hph" - this is the shared name of your file(s) that you want to process.
* desiredBaseLocation="/projects/ptx_results/2016/FusionData/10-Oct/Samples/CH/" - this is where the data analysis will be output.
* dataStorageLocation="/projects/ptx_results/2016/FusionData/10-Oct/Samples/CH/" - this is where the raw data are currently stored.

Below are the locations of the individual assets we need on my machine. You should change based on your machine.

* crapDatabase="/projects/ptx_analysis/chughes/databases/crap_jan2019.fasta" - this can be obtained from [here](ftp://ftp.thegpm.org/fasta/crap/).
* rawtoolsLocation="/projects/ptx_analysis/chughes/software/RawTools/rawtools202/RawTools.exe"
* searchGuiLocation="/projects/ptx_analysis/chughes/software/SearchGUI-3.3.16/SearchGUI-3.3.16.jar"
* peptideshakerLocation="/projects/ptx_analysis/chughes/software/PeptideShaker-1.16.42/PeptideShaker-1.16.42.jar"
* eval "export PATH=/gsc/software/linux-x86_64-centos7/R-3.6.3/bin:$PATH" - this is where your R distribution is installed.

For modifications, it is going to depend on how you prepared your sample. Below, the modifications are set up for a standard TMT analysis with carbamidomethylation of cysteine. If you need to add a modification and don't know the name, use the command `java -cp /path/to/searchGui.jar eu.isas.searchgui.cmd.IdentificationParametersCLI`.

* fixedModifications='"Carbamidomethylation of C, TMT 11-plex of peptide N-term"'
* variableModifications='"Oxidation of M, TMT 11-plex of K"'

The last thing to pay attention to is the memory size calls. For example:

* peptideshakerExecution="java -Xmx600g

Here we are allowing java to use 600gigs of RAM, which your system may not have. PeptideShaker can be extremely memory thirsty, especially when you are analyzing a large collection of files, so scale this based on your needs and system capabilities.

The last thing you will need to do is create an R script and tell the shell script where it is. The R script is provided below, so simply create a new file with vim as before, paste in the code, and save it with a '.R' extension. Once you have created the R script, you will need to change the following line in the shell script so it will be able to see it:

* rAnnotationIndex="Rscript /projects/ptx_analysis/chughes/projectsWorkspace/globalScripts/annotateUniprotDatabase.R

You are now ready to go. Execute the script as below.

```bash
chmod +x nameOfScript.sh
nameOfScript.sh |& tee dataProcessing.txt
```

## Shell script for data processing

```bash
#! /bin/bash

##############################################################
#this first part calls the help if the user triggers it
usage(){
printf "
This script will process proteomic data to generate PeptideShaker Reports.\n
The sample tag and desired base location parameters should be edited before running the script.\n
You sould execute this script in the shell scripting directory.\n
All results will be written to the folder created based on the sample tag and base location.\n\n"
}

if [[ $1 == -h ]] || [[ $1 == --help ]];
then
  usage
  exit
fi


##############################################################
##############################################################
#users must edit the below variables
#this is the text that will identify your files
sampleTag="10Oct2016_ADROIT_PrepCompare_TMT10_hph"
#this is the desired location for the output of the data processing process
desiredBaseLocation="/projects/ptx_results/2016/FusionData/10-Oct/Samples/CH/"
#this is the base location where your raw data is stored
dataStorageLocation="/projects/ptx_results/2016/FusionData/10-Oct/Samples/CH/"
#this is text that will be appended to your output folder that is created for this analysis
folderAdapter="dataProcessing_"
#users do not need to edit the statement below
folderToCreate="$desiredBaseLocation$folderAdapter$sampleTag"
#starting directory
startingDirectory="${PWD}"


##############################################################
##############################################################
#locations of software tools
crapDatabase="/projects/ptx_analysis/chughes/databases/crap_jan2019.fasta"
rawtoolsLocation="/projects/ptx_analysis/chughes/software/RawTools/rawtools202/RawTools.exe"
searchGuiLocation="/projects/ptx_analysis/chughes/software/SearchGUI-3.3.16/SearchGUI-3.3.16.jar"
peptideshakerLocation="/projects/ptx_analysis/chughes/software/PeptideShaker-1.16.42/PeptideShaker-1.16.42.jar"



##############################################################
##############################################################
####first we need to make the location for the data storage
makingDirectory="mkdir $folderToCreate"

#check if directory does not exist and create it if it doesn't
if [ ! -d "$folderToCreate" ];
then
  printf "Creating the processing directory.\n"
  eval $makingDirectory
fi

#move into the processing directory
eval "cd $folderToCreate"
printf "\nData processing proceeding in the directory: $PWD\n\n"



##############################################################
##############################################################
####set up the database to use
uniprotVersion=$( date +%b%Y )
databaseExtension="TargetDecoy"

#check if the database exists, and if not, get it
if [ ! -f "./uniprotHumanCrap$databaseExtension$uniprotVersion.fasta" ]; then
    printf "Database not found!\n\n"
    eval "wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640_9606.fasta.gz"
    eval "gunzip UP000005640_9606.fasta.gz"
    eval "mv UP000005640_9606.fasta ./uniprotHuman$uniprotVersion.fasta"
    eval "cat ./uniprotHuman$uniprotVersion.fasta $crapDatabase > ./uniprotHumanCrap$uniprotVersion.fasta"
    eval "java -cp $searchGuiLocation eu.isas.searchgui.cmd.FastaCLI -in ./uniprotHumanCrap$uniprotVersion.fasta -decoy"
    eval "mv ./uniprotHumanCrap${uniprotVersion}_concatenated_target_decoy.fasta ./uniprotHumanCrap$databaseExtension$uniprotVersion.fasta"
    eval "scp ./uniprotHuman$uniprotVersion.fasta /projects/ptx_analysis/chughes/databases/uniprot/"
    eval "scp ./uniprotHumanCrap$uniprotVersion.fasta /projects/ptx_analysis/chughes/databases/uniprot/"
    eval "scp ./uniprotHumanCrap$databaseExtension$uniprotVersion.fasta /projects/ptx_analysis/chughes/databases/uniprot/"
    printf "Database created and ready to use.\n\n"
else
  printf "\nDatabase file already exists: ./uniprotHumanCrap$databaseExtension$uniprotVersion.fasta\n\n"
fi

####process the database to make an annotation set with R
#just in case you are in a screen session, tell the screen session where R is
eval "export PATH=/gsc/software/linux-x86_64-centos7/R-3.6.3/bin:$PATH"

#check if the the annotation index exists, and if not, create it
if [ ! -f "uniprotHuman$uniprotVersion.fasta.annotated.rds" ]; then
  targetId="uniprotHuman$uniprotVersion.fasta"
  printf "Annotation index for $targetId does not exist, so it will be created. This can be slow but will only happen once per index.\n\n"
  rAnnotationIndex="Rscript /projects/ptx_analysis/chughes/projectsWorkspace/globalScripts/annotateUniprotDatabase.R ${PWD}/uniprotHuman$uniprotVersion.fasta"
  eval $rAnnotationIndex
  eval "scp ${PWD}/uniprotHuman$uniprotVersion.fasta.annotated.rds /projects/ptx_analysis/chughes/databases/uniprot/"
fi


##############################################################
##############################################################
####set up the parameters file
fixedModifications='"Carbamidomethylation of C, TMT 11-plex of peptide N-term"'
variableModifications='"Oxidation of M, TMT 11-plex of K"'

#####call the parameters file creation
callingParameters='java -cp "$searchGuiLocation" eu.isas.searchgui.cmd.IdentificationParametersCLI -out ./databaseSearchParameters.par -db ./uniprotHumanCrap"$databaseExtension""$uniprotVersion".fasta -prec_tol 20 -frag_tol 0.5 -db_pi ./uniprotHumanCrap"$databaseExtension""$uniprotVersion".fasta -fixed_mods "$fixedModifications" -variable_mods "$variableModifications" -msgf_instrument 0 -msgf_protocol 4 -msgf_fragmentation 1'

####call the command
if [ ! -f "./databaseSearchParameters.par" ]; then
  printf "\nParameters file does not exist. Creating it.\n\n"
  eval "$callingParameters"
else
  printf "Parameters file exists: databaseSearchParameters.par.\n\n"
fi



##############################################################
##############################################################
####run rawtools to get MGF files
##first find the raw files you want to process
fileList=()
  while IFS= read -d $'\0' -r file ; do
    fileList=("${fileList[@]}" "$file")
  done < <(find $dataStorageLocation -name "*$sampleTag*.raw" -print0 | sort -z)

###now execute rawtools on each of these
for i in "${fileList[@]}"; do
  targetId=$(basename "$i" ".raw")
  if [ -f "$targetId.raw.mgf" ]; then
    printf "MGF output for $targetId already exists, skipping file.\n\n"
  else
    rawtoolsExecution="mono $rawtoolsLocation -f $i -q -r TMT11 -uxmR -o $PWD"
    eval $rawtoolsExecution
  fi
done


##############################################################
##############################################################
####run the searchGUI analysis
fileList=()
  while IFS= read -d $'\0' -r file ; do
    fileList=("${fileList[@]}" "$file")
  done < <(find $PWD -name "*.mgf" -print0 | sort -z)

###now execute on each of these
for i in "${fileList[@]}"; do
  if [ -f "${i::-4}_searchgui_out.zip" ]; then
    targetId=$(basename "$i" ".zip")
    printf "Search output for $targetId already exists, skipping file.\n\n"
  else
    searchGuiExecution="java -Xmx120g -cp $searchGuiLocation eu.isas.searchgui.cmd.SearchCLI -spectrum_files $i -output_folder $PWD -id_params ./databaseSearchParameters.par -output_option 1 -xtandem 1 -msgf 1 -comet 1"
    eval $searchGuiExecution
  fi
done


##############################################################
##############################################################
####run peptideshaker and write the reports
####first move the quant files
if [ ! -d "quantFiles" ]; then
  eval "mkdir quantFiles"
  eval "mv *.txt ./quantFiles"
fi

###now execute on each of these
while ! [ -f "n_${sampleTag}_1_Default_PSM_Report.txt" ]; do
  printf "\nNo peptide shaker reports output exists, so it will be created.\n"
  #peptide shaker execution
  peptideshakerExecution="java -Xmx600g -cp $peptideshakerLocation eu.isas.peptideshaker.cmd.PeptideShakerCLI -experiment n -sample $sampleTag -replicate 1 -identification_files $PWD -spectrum_files $PWD -id_params ./databaseSearchParameters.par -out $PWD/$sampleTag.out.cpsx"
  eval $peptideshakerExecution
  #reports output
  reportsExecution="java -Xmx600g -cp $peptideshakerLocation eu.isas.peptideshaker.cmd.ReportCLI -in $sampleTag.out.cpsx -out_reports $PWD -reports 0,3"
  eval $reportsExecution
done
printf "\nFinished all exporting reports for $sampleTag data set.\n\n"


####final output
eval "rm *.html"
eval "mv ${startingDirectory}/dataProcessing.txt $PWD"
printf "\n\nFinished processing a total of ${#fileList[@]} files!\n\n"

```

## R script for database annotation

```R
#!/usr/bin/env Rscript

##libraries
suppressMessages(require(Biostrings))
suppressMessages(require(tidyverse))
suppressMessages(require(Peptides))
suppressMessages(require(OrgMassSpecR))

##get the input from the command line
commandLineInputs = commandArgs(trailingOnly=TRUE)

##process the fasta database
fastaDb = readAAStringSet(commandLineInputs[1])
#fastaDb = readAAStringSet('C:/Users/chris/OneDrive/Documents/bccrc/databases/uniprot/uniprotHuman202006.fasta')
fastaIndex = tibble('metadata' = names(fastaDb)) %>%
  mutate(accession = sub(".*[sptr]\\|(.*)\\|.*$", "\\1", metadata)) %>%
  mutate(length = width(fastaDb))
fastaIndex$gene = ifelse(grepl('GN=', fastaIndex$metadata), sub(".*GN=(.*) [PE].*$", "\\1", fastaIndex$metadata), NA)
fastaIndex$species = ifelse(grepl('sapiens', fastaIndex$metadata), 'human', 
                             ifelse(grepl('\\=9606', fastaIndex$metadata), 'human',
                                    ifelse(grepl('musculus', fastaIndex$metadata), 'mouse', 'other')))
fastaIndexFinal = dplyr::select(fastaIndex, accession, gene, species, length)

##get info about tryptic peptides
trypticPeptides = vector()
detectableLength = vector()
for (i in 1:length(fastaDb)){
  aaSeq = as.character(fastaDb[[i]])
  seqDigest = OrgMassSpecR::Digest(aaSeq, enzyme = 'trypsin', missed = 0, custom = list(code = c('X','U','Z','B'), mass = c(50, 60, 70, 80)))
  seqDigest$pepLength = (seqDigest$stop - seqDigest$start) + 1
  seqDigestSub = subset(seqDigest, (seqDigest$pepLength > 5) & (seqDigest$pepLength < 31))
  trypticPeptides = c(trypticPeptides, nrow(seqDigestSub))
  detectableLength = c(detectableLength, sum(seqDigestSub$pepLength, na.rm = TRUE))
  message(paste('Finished ', i, ' proteins.', sep = ''))
}

#now add to the database
fastaIndexFinal$detectablePeptides = trypticPeptides
fastaIndexFinal$detectableLength = detectableLength

##save the index
saveRDS(fastaIndexFinal, paste(commandLineInputs[1],'.annotated.rds',sep = ''))

```
