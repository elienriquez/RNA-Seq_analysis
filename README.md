# RNA-Seq_analysis
Here I describe the main steps for T. atroviride RNA-seq data analysis using data obtained from ncbi geo website

# Objetive
The aim of this repository is to distribute the pipeline used for T. atroviride RNA-seq data analysis.

# Biological relevance
Analyze the response of the fungus Trichoderma atroviride (wild type and RNAi-machinery-mutant strains) to fungal preys under three stages

# Sequencing Technology used and libraries obtained
90 libraries were sequenced by Illumina TruSeq 1X100 single-end.

# Quality analysis and reads mapping
The quality of the RNA-seq libraries were analyzed by FastQC version 0.11.8. Around 9 millions of high-quality read per library were obtained. Cleaned reads were mapped to the new genome reference of T. atroviride IMI206040 (Atriztán-Hernández et al., in prep) using HISAT2 version 2.1.0. We use the code as next:

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


# Read quantification and Differential expressed Analysis (DE)
For mapping reads counting to each gene Rsubread package was used:

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

# To performe the DE analysis the EdgeR and Limma packages were used 

1. Filtering genes with more than five reads per libraries and present in at least two different libraries 
 
counts = counts[rowSums(cpm(counts) >= 2) >=5,]

2. From filtered data, DGElist object was made, the replicates are introduced by factor function

group <- factor( c(rep("1",3),rep("2",3),rep("3",3))) 
design <-model.matrix(~ group)
set1<-DGEList(counts=counts,group=group)

3. Calculate dispersion and normalization by using calcNormFactors and estimateDisp functions

set1 <- calcNormFactors(set1)
set1 <- estimateDisp(set1,design)

4. Getting DE genes using Likelihood ratio test

fit     <- glmFit(set1, design)
lrt_fit  <- glmLRT(fit, coef=2)
lrt_fit2 <- glmLRT(fit, coef=3)
 
5. To retrieve the list of DE genes by using decideTests function with p.value = 0.05 and LogFoldChange = 1

res<-decideTests(lrt_fit,p.value=0.05,lfc=1)

6. Plotting Venn Diagrams for upregulated and downregulated genes 

res1<-decideTests(lrt_fit,p.value=0.05,lfc=1)
res2<-decideTests(lrt_fit2,p.value=0.05,lfc=1)
 
x_object           <- cbind(res1[,1], res2[,1])
colnames(x_object) <-c("condition_name_1","condiction_name_2")
#upregulated
vennDiagram(x_object [,1:2] == 1,circle.col = c("red","blue"),cex = c(1.5))
#downregulated
vennDiagram(x_object [,1:2] == -1,circle.col = c("red","blue"),cex = c(1.5))


# Enrichment analysis by using topGO package

#Merging a list of DE genes to a functional annotation file containing GO terms for each gene of the organism used, in this case T. atroviride

genes_to_GO <-merge("file_with_DE_genes","functional_annotation_file",by.x = "gene_ID",by.y = "gene_ID")
write.csv(genes_to_GO  , file = "path_to_save_file_in_.csv")

#Adding the annotations by "readMappings" function

genes2GO<-readMappings(file = "path_to_genes_to_GO",sep = "\t", IDsep = ";")

#Assigning pvalues to a geneList

GeneList <-file_with_DE_genes[,"number_of_column_containing_pvalues"]
names(GeneList) <-file_with_DE_genes[,"number_of_column_containing_geneIDs"]

#Building the topGO object specifying the type of ontology (MF, BP or CC)

GOdata  <- new("topGOdata",
               description = "Simple session", ontology = "MF",
               allGenes =GeneList ,geneSel = topDiffGenes,
               nodeSize = 10,
               annot=annFUN.gene2GO,gene2GO = genes2GO)
               
#Performing the Fisher Test, KS or KS_elim by using "runTest" function
 
Fisher  <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
KS      <- runTest(GOdata, algorithm = "classic", statistic = "ks")
KS_elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")

#Retrieve the results of the tests performed by "GenTable" function, in this case for the top 10 significant GO terms

Genetable<-GenTable(GOdata, classicFisher = Fisher,classicKS = KS, elimKS = KS_elim,orderBy = "elimKS", ranksOf = "classicFisher", topNodes =10)

# Upset graphs 

Here we describe the procedure to create upset graphs. The files needed for this section were provided as "contras_up.txt" and "contras_down.txt". The lists of upregulated and downregulated differentially expressed genes (DEG) with
a Log2FC  1 and p &lt; 0.05 were taken for each interaction with respect to the control (wild strain growing without interaction). These lists were loaded to create
the UpSet charts using the following code:

BiocManager::install(&quot;UpSetR&quot;)
library(&quot;UpSetR&quot;)
library(&quot;ggplot2&quot;)
library(&quot;plyr&quot;)
library(&quot;gridExtra&quot;)
library(&quot;grid&quot;)

For upregulated genes:

data_inter_up=fromList(lista_contras$contras_up.txt)
UpSet plot of genes that are upregulated at the before contact and contact stage:
upset(data_inter_up, sets = c(&quot;WTvsR.solani_AG2_BC&quot;,
&quot;WTvsR.solani_AG5_BC&quot;, &quot;WTvsA.alternata_BC&quot;,&quot;WTvsR.solani_AG2_C&quot;,
&quot;WTvsR.solani_AG5_C&quot;, &quot;WTvsA.alternata_C&quot;),
main.bar.color = &quot;navyblue&quot;, sets.bar.color = &quot;red&quot;, point.size = 4.5,
number.angles = 10, text.scale = c(rep(1.5, 5), 1.5),
set_size.show = TRUE, order.by = &quot;freq&quot;, keep.order = TRUE,
set_size.scale_max = 1300, matrix.color = &quot;SteelBlue&quot;, nintersects = 15,
mb.ratio = c(0.55, 0.45),
queries = list(list(query = intersects, params = list(&quot;WTvsR.solani_AG2_BC&quot;,
&quot;WTvsR.solani_AG5_BC&quot;, &quot;WTvsA.alternata_BC&quot;),
color =&quot;chartreuse3&quot;, active = T),list(query = intersects, params =
list(&quot;WTvsR.solani_AG2_C&quot;, &quot;WTvsR.solani_AG5_C&quot;, &quot;WTvsA.alternata_C&quot;),color
= &quot;chartreuse3&quot;, active = T)))
UpSet plot of genes that are upregulated at the after contact stage:
upset(data_inter_up, sets = c(&quot;WTvsR.solani_AG2_AC&quot;,
&quot;WTvsR.solani_AG5_AC&quot;, &quot;WTvsA.alternata_AC&quot;),
main.bar.color = &quot;navyblue&quot;, sets.bar.color = &quot;red&quot;, point.size = 4.5,
number.angles = 10, text.scale = c(rep(1.5, 5), 1.5),
set_size.show = TRUE, order.by = &quot;freq&quot;, keep.order = TRUE,
set_size.scale_max = 100, matrix.color = &quot;SteelBlue&quot;,
mb.ratio = c(0.55, 0.45),
queries = list(list(query = intersects, params = list(&quot;WTvsR.solani_AG2_AC&quot;,
&quot;WTvsR.solani_AG5_AC&quot;, &quot;WTvsA.alternata_AC&quot;),
color =&quot;chartreuse3&quot;, active = T)))

For downregulated genes:

data_inter_down=fromList(lista_contras$contras_down.txt)
UpSet plot of genes that are downregulated at the before contact and contact
stage:
upset(data_inter_down, sets = c(&quot;WTvsR.solani_AG2_BC&quot;,
&quot;WTvsR.solani_AG5_BC&quot;, &quot;WTvsA.alternata_BC&quot;,&quot;WTvsR.solani_AG2_C&quot;,
&quot;WTvsR.solani_AG5_C&quot;, &quot;WTvsA.alternata_C&quot;),
main.bar.color = &quot;navyblue&quot;, sets.bar.color = &quot;red&quot;, point.size = 4.5,
number.angles = 10, text.scale = c(rep(1.5, 5), 1.5),
set_size.show = TRUE, order.by = &quot;freq&quot;, keep.order = TRUE,
set_size.scale_max = 1300, matrix.color = &quot;SteelBlue&quot;, nintersects = 15,
mb.ratio = c(0.55, 0.45), queries = list(list(query = intersects, params =
list(&quot;WTvsR.solani_AG2_BC&quot;, &quot;WTvsR.solani_AG5_BC&quot;, &quot;WTvsA.alternata_BC&quot;),
color =&quot;chartreuse3&quot;, active = T)))
UpSet plot of genes that are downregulated at the after contact stage:
upset(data_inter_down, sets = c(&quot;WTvsR.solani_AG2_AC&quot;,
&quot;WTvsR.solani_AG5_AC&quot;, &quot;WTvsA.alternata_AC&quot;),
main.bar.color = &quot;navyblue&quot;, sets.bar.color = &quot;red&quot;, point.size = 4.5,
number.angles = 10, text.scale = c(rep(1.5, 5), 1.5),
set_size.show = TRUE, order.by = &quot;freq&quot;, keep.order = TRUE,
set_size.scale_max = 100, matrix.color = &quot;SteelBlue&quot;, nintersects = 15,
mb.ratio = c(0.55, 0.45),
queries = list(list(query = intersects, params = list(&quot;WTvsR.solani_AG2_AC&quot;,
&quot;WTvsR.solani_AG5_AC&quot;, &quot;WTvsA.alternata_AC&quot;),
color =&quot;chartreuse3&quot;, active = T)))
 
