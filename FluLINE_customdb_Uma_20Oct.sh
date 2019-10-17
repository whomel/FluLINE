###############################################################
#  Script     : FluAnalysisPipeline.sh
#  Author     : Uma Sangumathi
#  Date       : 28/04/2015
#  Last Edited: 17/10/2019, miguel 
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

#main folder
hmdir=/home/miguelg/bm14_scratch/miguel/who/Uma_FluLINE

# Path to softwares used : Please edit this if "bin" and "src" folders are not in the current working directory
SCRIPTS_DIR=$hmdir/bin_Custom_2018
SOFTWARE_DIR=$hmdir/src
#PYTHONPATH=$hmdir/src/dist-packages

=======
fqdir=/home/miguelg/bm14_scratch/miguel/who/flub/DanRunLib1912-16 #Enter the fastq files directory
refdir=$hmdir/Reference  
blastn=/home/miguelg/bm14/apps/ncbi-blast-2.7.1+/bin/blastn # installed blastn location
nt_db=/scratch/blastdb/FASTA/nt # NCBI database location
out_dir=$hmdir/out_DanRunLib1912-16 #Enter the Output directory
info_file=/home/miguelg/bm14_scratch/miguel/who/flub/Lib1912_16info.csv  # Sample Key
name_replace=.DanLib1912_16.fastq
email_id=miguel.graulopez@monash.edu #Enter your emailid

Nextension=0


export PATH=$SCRIPTS_DIR:$PATH
export PATH=$SOFTWARE_DIR:$PATH
#export PICARDDIR=/home/user/Desktop/FluLINE/src/picard/build/libs
#module load picard/2.9.2
#export PATH=$PATH:/home/user/Downloads/ncbi-blast-2.7.1+/bin  # Edit to installed blastn location


## set DATABASE 
# Make custom DB
# DOwnload from GISAID, DNA fasta = Type|Segment|DNA Accession no.|
# /home/administrator/Deskto2018_p/Package_FluAnalysis/fasta_tools-master/bin/fasta_unique ./Custom_DB/Fludb_19Jun2018.fasta > ./Custom_DB/Fludb_19Jun2018-unique.fasta
# makeblastdb -in ./Custom_DB/Fludb_19Jun2018-unique.fasta -parse_seqids -dbtype nucl

database_mode=NCBI #Custom ## takes either "NCBI" or "Custom" 
BLASTDB=/scratch/blastdb/FASTA/nt  #$hmdir/TEST/Btest_1.fasta #/home/administrator/Downloads/ncbi-blast-2.2.31+/db/nt/nt
#genome_coords=/home/administrator/Desktop/Package_FluAnalysis/RSV/RSV_coords.txt


############## Do not edit below this unless you want to change the pipeline flow ##########################


cd $fqdir

echo "Run QC on the raw fastq files.."
run_QC.sh  $fqdir $name_replace  

for x in *-Qc.fastq
do 
  echo $x
  partial=$(echo $x | sed 's/-Qc.fastq//g' ) 

  echo "Generate consensus genome and map reads using Bowtie..."
  GenerateConsensusGenome_withBlast_custom.py -r $refdir  -q $fqdir -f $x  -b $blastn -d $nt_db -o $out_dir -i $info_file -e $email_id -t $database_mode
 
  echo "Generate files required to make the graphs...."
  createGraphfiles_Full_17Nov.py  $out_dir/$partial
  #R -q --no-save --no-restore --slave -f $SCRIPTS_DIR/plot_linear_cov.R --args $genome_coords $out_dir/$partial/*-sort.coverage
  #createGraphfiles_Full_17Nov.py for full genome graph

  echo "Improve the 5' and 3' N by replacing the nts from softclipping..."
  N_consensus.py $Nextension $out_dir/$partial 

done 

echo "Generate Graphs..."
mkdir $out_dir/circos 
cp  $out_dir/*/*txt   $out_dir/circos
cp $SCRIPTS_DIR/*conf  $out_dir/circos
cp -r $hmdir/src/etc $out_dir/circos/
generate_covplot.py  $info_file  $out_dir/circos $SCRIPTS_DIR

echo "Copy results to separate folders.."
MoveToFolders.py $out_dir $fqdir $info_file
  
echo "Combine all the annotation files of the consensus genome..."
combine.annot.py  $out_dir/Summary_result  $info_file

###### The END ########



