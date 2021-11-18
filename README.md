# RNA-Seq_analysis
Here I describe the main steps for T. atroviride RNA-seq data analysis using data obtained from ncbi geo website

# Objetive
# The aim of this repository is to distribute the pipeline used for T. atroviride RNA-seq data analysis.

# Biological relevance
# Analyze the response of the fungus Trichoderma atroviride (wild type and RNAi-machinery-mutant strains) to fungal preys under three stages

# Sequencing Technology used and libraries obtained
# 90 libraries were sequenced by Illumina TruSeq 1X100 single-end.

#Quality analysis and reads mapping
#The quality of the RNA-seq libraries were analyzed by FastQC version 0.11.8. Around 9 millions of high-quality read per library were obtained. Cleaned reads were mapped to the new genome reference of T. atroviride IMI206040 (Atriztán-Hernández et al., in prep) using HISAT2 version 2.1.0. We use the code as next:

module load FastQC/0.11.8
fastqc /path-to-.fastq.gz* --outdir /path-to-outputdirectory

#Trimmomatic
module load Trimmomatic/0.32
java -jar $TRIMMOMATIC SE -phred33 “path.to-.fastq.gz” “outputdirectory” ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

#HISAT2
module load  hisat2/2.1.0
#create index
hisat2-build -p 8 -f “genome.fasta”   “Output-directory”
#aligning
hisat2 -p 4 -t “path to index” -U “path to trimed-.fq.gz” --dta -S “path-to output-.sam”

#Read quantification and Differential expressed Analysis (DE)

#For mapping reads counting to each gene Rsubread package was used:

counts<-featureCounts("path_to_bam_file",
                         annot.ext = "path_to_annotation_file_in_.gff3",
                         isGTFAnnotationFile = TRUE,
                         GTF.featureType = "mRNA",
                         GTF.attrType = "ID",
                         countMultiMappingReads = TRUE,
                         fraction = TRUE,
                         allowMultiOverlap = TRUE)
matrix_counts_file<-cbind(counts$counts,counts2$counts,countsN$counts)                         
write.csv(matrix_counts_file, file = "path_to_save_the_matrix_in_.csv")




