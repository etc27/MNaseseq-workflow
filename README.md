# Workflow for generating nucleosome occupancy plots from MNase-seq data
Emma Tung Corcoran (10/13/2020)

## Introduction
This document covers my basic workflow for processing paired-end MNase-seq (micrococcal nuclease digestion with deep sequencing) samples and assaying nucleosome occupancy at different genomic loci. I used the [Ruddle HPC cluster at the Yale Center for Research Computing](https://docs.ycrc.yale.edu/clusters-at-yale/clusters/ruddle/) for my HPC environment.

## 1. Setup

### Set up Miniconda environment
I followed the directions on the [YCRC Conda Documentation](https://docs.ycrc.yale.edu/clusters-at-yale/guides/conda/) to set up Miniconda as a module on my HPC account. First, I generated a conda environment (name=env_name) containing all of the packages I need for MNase-seq that are not included with the default conda installation.

```
module load miniconda
conda create -n env_name fastqc trim-galore subread multiqc samtools bowtie2 bedops
```

Then, I am able to load the conda environment containing all of the required packages with the following code.

```
module load miniconda
conda activate env_name
```
### Set up the folder structure
In order to organize all of the files generated from processing the RNA-seq raw data, I utilized the following folder structure adapted from [this RNA-seq workflow](https://github.com/twbattaglia/RNAseq-workflow).

```
── MNaseseq_data/
  │   └── annotation/               <- Genome annotation file (.GTF/.GFF)
  │  
  │   └── genome/                   <- Reference genome file (.FASTA)
  │  
  │   └── input/                    <- Location of input MNase-seq data
  │  
  │   └── results/                  <- Data generated during processing steps
  │       ├── 1_initial_qc/         <- Quality check of input files
  │       ├── 2_trimmed_output/     <- Trimmed read files and quality check of trimmed reads
  │       ├── 3_aligned_sequences/  <- Main alignment files for each sample
  │           ├── aligned_bam/      <- Alignment files generated from bowtie2 (.BAM)
  │           ├── aligned_logs/     <- Log from running bowtie2 alignment step
  │           ├── aligned_sam/      <- Alignment files generated from bowtie2 (.SAM)
  │       ├── 4_multiQC/            <- Overall report of logs for each step
  │  
  │   └── bowtie2/                  <- Folder to store the indexed genome files from bowtie2
  │    
  │   └── plot2DO/                  <- Folder to store plot2DO script and output
  │       ├── annotations/          <- Annotations folder for plot2D
  │       ├── config/               <- Config folder for plot2DO (don't edit)
  │       ├── misc/                 <- Examples of plot2DO output
  │       ├── output/               <- Output from plot2D
  │       ├── sample_BAM_files/     <- Sample BAM files provided for plot2DO workflow
  │       ├── source/               <- Source code for plot2DO (don't edit)
  ```
  
### Download the reference genome and annotation
I downloaded the *Arabidopsis thaliana* reference genome (Araport 11) from the [JGI Genome Porta](https://genome.jgi.doe.gov/portal/pages/dynamicOrganismDownload.jsf?organism=Athaliana) to the `genome/` folder. The genome assembly was called `Athaliana_447_TAIR10.fa.gz`
I downloaded the *Arabidopsis thaliana* annotation (Araport 11) from [TAIR](https://www.arabidopsis.org/download/index-auto.jsp?dir=%2Fdownload_files%2FGenes%2FAraport11_genome_release) to the `annotation/` folder. The annotation was called `Araport11_GFF3_genes_transposons.201606.gtf`. I also downloaded `Araport11_gene_type`, a list of the gene types corresponding to each gene.

### Download raw sequencing data
In order to access the sequencing data from the MNase-seq experiments, I followed the directions on the [Ruddle documentation](https://docs.ycrc.yale.edu/clusters-at-yale/clusters/ruddle/#access-sequencing-data). Briefly, the Yale Center for Genome Analysis sent me a url that looks like this: 
`http://fcb.ycga.yale.edu:3010/randomstring/sample_dir_001` and I used the ycgaFastq tool to make a soft link to the data with the following command.
```
/home/bioinfo/software/knightlab/bin_Mar2018/ycgaFastq  fcb.ycga.yale.edu:3010/randomstring/sample_dir_001
```
My raw sequencing data contains paired-end Illumina sequencing reads. Note that the fastq files corresponding to pair1 and pair2 for one sample are labeled `sample_R1_001.fastq.gz` and `sample_R2_001.fastq.gz`

## 2. Analyze sequence quality with FastQC

### Description
[FastQC: A quality control tool for high throughput sequence data.](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
"FastQC aims to provide a simple way to do some quality control checks on raw sequence data coming from high throughput sequencing pipelines. It provides a modular set of analyses which you can use to give a quick impression of whether your data has any problems of which you should be aware before doing any further analysis."

### Command
Note: it is not necessary to unzip the fastq file for any of the following processes.
```
# Run FastQC
# -o: output directory
# --noextract: do not uncompress the output file after creating it

fastqc -o results/1_initial_qc/ --noextract input/sample.fastq.gz
```

### Output
```
── results/1_initial_qc/
    └──  sample_fastqc.html   <- HTML file of FastQC quality analysis figures
    └──  sample_fastqc.zip    <- FastQC report data
```

### Explanation of FastQC figures
The [FastQC Manual](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) provides detailed explanations of each figure generated. Additionally, [this webpage](https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon/lessons/qc_fastqc_assessment.html) gives a good explanation of the files generated by FastQC and their significance for RNA-seq.

#### Sequence Quality Histograms and Per Sequence Quality Scores
Two of the most important figures for quality analysis are the **Sequence Quality Histograms**, which provide the mean quality value across each base position in the read, and the **Per Sequence Quality Scores**, which display the average quality score on the x-axis by the number of sequences with that average on the y-axis.

#### Per Base Sequence Content
The **Per Base Sequence Content** shows the proportion of each base position for which each of the four normal DNA bases has been called.

#### Per Sequence GC Content and Overrepresented Sequences
The **Per Sequence GC Content** measures the average GC content of reads and compares it to a modelled normal distribution of GC content. In a normal random library you should expect to see a roughly normal distribution of GC content where the central peak corresponds to the expected GC content for the organism. The **Overrepresented Sequences** shows the total amount of overrepresented sequences found in each library (if amount is over 1%). If the **Per Sequence GC Content** in the previous module looked suspect, this table can help identify the source. If the Possible Source column shows No Hit, you can BLAST the sequence to determine the identity.

#### Sequence Length Distribution
The **Sequence Length Distribution** shows the distribution of fragment sizes (read lengths). Note that this module may show **FAIL** after trimming and adapter cleaning if you removed reads shorter than a certain length during that process, but the warning is safe to ignore in that case.

#### Sequence Duplication Levels
Finally, the **Sequence Duplication Levels** displays the relative level of duplication found for every sequence. [This page](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/8%20Duplicate%20Sequences.html) provides a detailed explanation of how to interpret these data.

## 3. Perform quality control with Trim Galore

### Description
[Trim Galore: A wrapper tool around Cutadapt and FastQC to consistently apply quality and adapter trimming to FastQ files.](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
Trim Galore performs adapter trimming (the default is the first 13 bp of Illumina standard adapters ('AGATCGGAAGAGC'), but it is able to autodetect the adapter sequence). It also removes sequences that become too short during the trimming process. With the `--paired` option, Trim Galore removes both reads in a pair if at least one of the two sequences becomes shorter than the threshold. Additionally, it can run FastQC on the output files to assess quality once trimming has been completed. I kept the default options for quality and length.

### Command
```
# Run Trim Galore! (input is both paired end sequencing files for a sample)
#--paired: remove both reads in a pair if at least one of the two sequences becomes shorter than the threshold
#--fastqc: run FastQC in the default mode on the FastQ file once trimming is complete
#--output_dir: output directory

trim_galore --paired --fastqc --output_dir results/2_trimmed_output/ input/sample_R1_001.fastq.gz input/sample_R2_001.fastq.gz
```

### Output
```
── results/2_trimmed_output/
     └──  sample_R1_001_val_1.fq.gz                          <- Compressed trimmed sequencing file (for read1)
     └──  sample_R1_001_val_1_fastqc.html                    <- HTML file of FastQC quality analysis figures (for read1)
     └──  sample_R1_001_val_1_fastqc.zip                     <- FastQC report data (for read1)
     └──  sample_R1_001.fastq.gz_trimming_report.txt         <- Cutadapt trimming report (for read1)
     └──  sample_R2_001_val_2.fq.gz                          <- Compressed trimmed sequencing file (for read2)
     └──  sample_R2_001_val_2_fastqc.html                    <- HTML file of FastQC quality analysis figures (for read2)
     └──  sample_R2_001_val_2_fastqc.zip                     <- FastQC report data (for read2)
     └──  sample_R2_001.fastq.gz_trimming_report.txt         <- Cutadapt trimming report (for read2)
```

## 4. Align reads to genome with bowtie2

### Description
[Bowtie 2: Fast gapped-read alignment.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3322381/)
[Bowtie 2](https://github.com/BenLangmead/bowtie2) is an ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences.

### Generating index
Before running bowtie2, it is necessary to generate an index of your reference genome. I stored this index in the `bowtie2` folder. It is only necessary to run this step once.
```
# Build Arabidopsis thaliana genome index for bowtie2
# first argument specifies path to reference genome (FASTA)
# second argument specifies output prefix

bowtie-build genome/Athaliana_447_TAIR10.fa.gz bowtie2/Athaliana_447_TAIR10
```

### Run bowtie2
```
# Run bowtie2
#--very-sensitive: preset parameters that are optimal for MNase-seq (more sensitive, but slower)
#-x: specifies path to the directory where the genome indices are stored
#-1: path to sample read1 file
#-2: path to sample read2 file
#-S: path to output SAM file
#-p: specified number of parallel search threads (4)

bowtie2 --very-sensitive -x ../../bowtie2/Athaliana_447_TAIR10 -1 results/2_trimmed_output/sample_R1_001_val_1.fq.gz \
-2 results/2_trimmed_output/sample_R2_001_val_2.fq.gz -S results/4_aligned_sequences/sample.sam -p 4
```

### Generate log of alignment (summary mapping statistics) and .bam file
```
# Generate .bam file
samtools view -bS sample.sam > sample.bam

# Generate log of alignment
samtools stats sample.sam > sample.txt
```

### Move output files to correct subdirectories
```
# Move the BAM file into the correct folder
mv sample.bam aligned_bam

# Move the log file into the correct folder
mv sample.txt aligned_logs

# Move the SAM file into the correct folder
mv sample.sam aligned_sam
```

### Output
```
── results/3_aligned_sequences/
    └── aligned_bam/sample.bam           <- BAM alignment
    └── aligned_logs/sample.txt          <- Summary mapping statistics
    └── aligned_sam/sample.sam           <- SAM alignment
```

## 5. Generate analysis report with MultiQC

### Description
[MultiQC: summarize analysis results for multiple tools and samples in a single report.](https://pubmed.ncbi.nlm.nih.gov/27312411/)
MultiQC is a tool to create a single report visualizing output from multiple tools across many samples, including the output of FastQC, Trim_Galore, and Bowtie 2.

### Command
```
# Run multiQC
#--outdir: output directory
multiqc results --outdir results/4_multiQC
```

### Output
```
── results/4_multiQC/
    └── multiqc_report.html     <- Figures representing the logs from each step
    └── multiqc_data/           <- Folder of data that multiqc found from various log files
```

### Explanation of MultiQC Figures
MultiQC summarizes results from FastQC before and after trimming, as well as logs from the alignment and gene counts steps. Many of the figures in MultiQC are the same or similar to those produced from FastQC (see FastQC section above). One new figure in the FastQC section is **Sequence Counts**, displaying the sequence counts for each sample divided into unique reads and estimated duplicate reads. The ratio of unique to duplicate reads gives information about library complexity vs. sequencing depth.

## 6. Generate nucleosome occupancy plots using Plot2DO

### Description
[Plot2DO: Creating 2D Occupancy Plots](https://link.springer.com/protocol/10.1007/978-1-0716-0301-7_5). 
[Plot2DO](https://github.com/rchereji/plot2DO), is an open source package written in R for evaluating the quality of MNase-seq and MNase-ChIP-seq data and for visualizing nucleosome distributions.

### Download plot2DO
Follow the instructions on the [plot2DO GitHub page](https://github.com/rchereji/plot2DO) to install plot2DO. You need to have R (version >= 3.5.0) installed on your computer. Briefly, run the following command to download the plot2DO package:
```
git clone https://github.com/rchereji/plot2DO.git
```
Then, run the following command to install other dependencies:
```
Rscript plot2DO_setup.R
```

### Command
Run the following script within the plot2DO folder to execute plot2DO. By default, plot2DO aligns to the TSS as the reference point and considers DNA fragments between lengths 50-200 bp.
```
# Run plot2DO
#-f: input bam file
#-g: genome (tair10 for Arabidopsis thaliana)
Rscript plot2DO.R -f ../results/4_aligned_sequences/aligned_bam/sample.bam -g tair10
```

### Output
```
── plot2DO/output/
    └── OCC_matrix.TSS.50_200.sample.pdf             <- The three panels generated by plot2DO: (1) 2D occupancy plot (2) One-dimensional occupancy plot (3) DNA fragment length histogram
    └── OCC_matrix.TSS.50_200.sample.RData           <- R Data generated from plot2DO execution
```

## 7. Using plot2DO to generate plots for custom lists of genes
By default, plot2DO computes the coverage of all sites across the genome, but it is possible to have plot2DO compute the coverage for custom lists of sites. For *Arabidopsis thaliana*, I am interested in computing the coverage for all protein-coding genes.

### Generating custom list of sites in bed format
Plot2DO accepts a list of sites in bed format as a parameter for alignment. First, I generated this list of sites using the command line and R. The following example demonstrates using this workflow for a list of all protein-coding genes in *Arabidopsis thaliana*. R files and annotation files can be found within the annotation/ folder.

### Commands run in command line (part 1)
Navigate to annotation/ folder. Save the rows of the Araport11 annotation corresponding to transcripts as Araport_transcripts.gtf
```
awk '$3 == "transcript" {print}' Araport11_GFF3_genes_transposons.201606.gtf > Araport_transcripts.gtf
```
Save the gene names of Araport11_gene_type.txt corresponding to protein-coding genes as proteincoding.txt (list of genes)
```
awk '$2 == "protein_coding" {print $1}' Araport11_gene_type.txt > proteincoding.txt
```

### Commands run in R
Run **nonredundantGenes.R** to remove redundant genes from Araport_transcripts.gtf
```
#first parameter: input gtf file
#second parameter: output file name
Rscript nonredundantGenes.R "Araport_transcripts.gtf" "Araport_transcripts_nonredundant.gtf"
```
Run **geneSublist.R** to generate custom list of genes
```
#first parameter: input gtf file
#second parameter: list of gene names
#third parameter: output file name
Rscript geneSublist.R "Araport_transcripts_nonredundant.gtf" "proteincoding.txt" "proteincoding_sublist.gtf"
```

### Commands run in command line (part 2)
Copy custom list of genes to plot2DO/annotations/ folder.
```
cp annotation/proteincoding_sublist.gtf plot2DO/annotations
```
Convert list of genes to bed format.
```
gtf2bed < proteincoding_sublist.gtf > proteincoding_full.bed
```
Keep only the first 6 columns of bed file so that it is compatible with plot2DO.
```
awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6}' proteincoding_full.bed > proteincoding.bed
```

### Run plot2DO with custom list of sites
```
# Run plot2DO
#-f: input bam file
#-g: genome (tair10 for Arabidopsis thaliana)
#--sites: user-provided sites to be aligned (BED file)
#--align: points of the provided intervals to be aligned
#--minLength: the smallest DNA fragment to be considered [default = 50]
#--maxLength: the largest DNA fragment to be considered [default = 200]
#--siteLabel: label for the aligned sites [default = Sites]
Rscript plot2DO.R -f ../results/4_aligned_sequences/aligned_bam/sample.bam -g tair10 \ 
--sites=annotations/proteincoding.bed --align=fivePrime --minLength=140 --maxLength=160 --siteLabel=proteincoding
```

### Output
```
── plot2DO/output/
    └── OCC_matrix.proteincoding.140_160.sample.pdf             <- The three panels generated by plot2DO: (1) 2D occupancy plot (2) One-dimensional occupancy plot (3) DNA fragment length histogram
    └── OCC_matrix.proteincoding.140_160.sample.RData           <- R Data generated from plot2DO execution
```
