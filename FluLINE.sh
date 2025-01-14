##############################################################
#  Script     : FluLINE.sh
#  Author     : Uma Sangumathi
#  Date       : 28/04/2015
#  Last Edited: 09/04/2017, uma 
#  Description: Flu Analysis pipeline: guided consensus genome and variants
##############################################################
# Purpose: 
#  Create guided consensus genome from the reads and nearest reference genome     
#  Map reads using bowtie2    
#  Snps detected using Lofreq2 
#  Plot coverage graph    
#  moves the result files 
# Requirements:
#  1. Biopython module    
#  2. Bwa
#  3. Bowtie2
#  4. components of Vipr pipeline  
#  5. Samtools 
#  6. NCBI blast commandline tool and database
#  7. picard-tools-1
#  8. Circos software for graph
#############################################################

# Path to softwares used : Please edit this if "bin" and "src" folders are not in the current working directory
SCRIPTS_DIR=./FluLINE/bin
SOFTWARE_DIR=./FluLINE/src
PYTHONPATH=./FluLINE/src/dist-packages
# Location of the data , the output


=======
hmdir=/home/administrator/Desktop/Package_FluAnalysis
fqdir=$hmdir/Fastq    #Enter the fastq files directory
refdir=$hmdir/Reference  
blastn=/home/administrator/Downloads/ncbi-blast-2.2.31+/bin/blastn # installed blastn location
nt_db=/home/administrator/Downloads/ncbi-blast-2.2.31+/db/nt/nt # NCBI database location
out_dir=$hmdir/Analysis  #Enter the Output directory
info_file=$hmdir/Sample_info.csv  # Sample Key
name_replace=.fastq
email_id=emailid@gmail.com #Enter your emailid

Nextension=0


export PATH=$SCRIPTS_DIR:$PATH
export PATH=$SOFTWARE_DIR:$PATH
export PICARDDIR=$SOFTWARE_DIR/picard-tools-1.105/
export PATH=$PATH:/home/administrator/Downloads/ncbi-blast-2.2.31+/bin  # Edit to installed blastn location

############## Do not edit below this unless you want to change the pipeline flow ##########################
cd $fqdir

echo "Run QC on the raw fastq files.."
run_QC.sh  $fqdir $name_replace  

for x in *.fastq
do 
  echo $x
  partial=$(echo $x | sed 's/.fastq//g' ) 

  echo "Generate consensus genome and map reads using Bowtie..."
  GenerateConsensusGenome_withBlast.py -r $refdir  -q $fqdir -f $x  -b $blastn -d $nt_db -o $out_dir -i $info_file -e $email_id
 
  echo "Generate files required to make the graphs...."
  createGraphfiles-full.py  $out_dir/$partial

  echo "Improve the 5' and 3' N by replacing the nts from softclipping..."
  N_consensus.py $Nextension $out_dir/$partial 

done 

echo "Generate Graphs..."
mkdir $out_dir/circos 
cp  $out_dir/*/*txt   $out_dir/circos
cp $SCRIPTS_DIR/*conf  $out_dir/circos
generate_covplot.py  $info_file  $out_dir/circos $SCRIPTS_DIR

echo "Copy results to separate folders.."
MoveToFolders.py $out_dir $fqdir $info_file
  
echo "Combine all the annotation files of the consensus genome..."
combine.annot.py  $out_dir/Summary_result  $info_file

###### The END ########



