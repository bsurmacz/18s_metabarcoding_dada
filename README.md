# 18 metabarcoding pipeline using Ares - for Daniel  
This repository contains a pipeline for the fast analysis of  18 metabarcoding.
The pipeline employs *dada2* R package (Callahan et al. 2016) and is based onthe official *dada2* tutorial (benjjneb.github.io/dada2/tutorial.html).
All calculations are performed using R 3.6.3. Additionally the pipeline uses *Cutadapt* 4.6 (Martin 2011) to remove primer sequences

OUTPUT:  
**<RUN_NAME>**_zOTU_TABLE.csv - zOTU table with classification (with taxonomy classification by dada2 AssignTaxonomy on custom database)    
**<RUN_NAME>**_zOTUs.fasta - fasta file with zOTU sequences

## log in to ares using your login in terminal
```
ssh <yourname>@ares.cyfronet.pl
```
## start screen session to make the script run even after logging out 
```
screen -S meta
```
(to go back to the session use screen -r meta)
## start interactive session in slurm using resources from the computing grant (for now: 10 hours of using 16 cores and 50G of RAM)
```
srun -N 1 -n 16 -p plgrid -A plgtardigrada-cpu -t 10:0:0 --mem=50G --pty /bin/bash -l
```

## load the modules
```
ml cutadapt/3.5-gcccore-11.2.0
ml r/4.4.2-foss-2023a

```

## go to your working directory or wherever you want
```
cd $SCRATCH

```
## start R:
```
R
```
In interactive R session: 
Loading required libraries
```
library(dada2)  
library(dplyr)
library(data.table)
library(seqinr)
library(tibble)
library(stringr)
```
 
### *INPUT_DIRECTORY* is a path to folder containing demultiplexed files with raw sequencing reads  
```
RUN_NAME <- "SUPER_NICE_DATASET"
MARKER_NAME<-"18S"
INPUT_DIRECTORY<-"PATH/TO_FOLDER_WITH/DATA"
DATABASE_PATH <- "/net/pr2/projects/plgrid/plggiecoteam/REFERENCE_DATABASES/18S/pr2_version_5.1.1_SSU_dada2_WITH_PEATBOGS_WITH_JEZIORKA.fasta.gz" # for taxonomy classification

### trimming of raw reads to discard low-quality read ends  (OPTIONAL!!!)
TRIM_LENGTH = 300

```
### Primer sequences used for metabarcoding:
```
  FWD<-  "GCTTGTCTCAAAGATTAAGCC" # check this
  REV<- "GCCTGCTGCCTTCCTTGGA"   # check this
 
```

### File input and N-filtering (filtering reads with uncalled bases)  
```
path <-paste0(INPUT_DIRECTORY)
list.files(path)   # list files to check if the directory is correct
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))  ## check if this pattern matches your files and adjust it!!!!
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE)) 	 

fnFs.filtN <- file.path(path, "filtN", basename(fnFs))		 
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, truncLen=TRIM_LENGTH,  multithread = TRUE)
```
### Setting flags for Cutadapt and trimming primers, filtering reads without primers
```
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))
	
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

CUTADAPT_FLAGS<-paste0("cutadapt -a ^",FWD,"...",REV.RC," -A ^",REV,"...",FWD.RC," --discard-untrimmed -o")

for(i in seq_along(fnFs)) {
	system(paste(CUTADAPT_FLAGS, fnFs.cut[i],"-p",fnRs.cut[i],fnFs[i], fnRs[i]))
}
```	
 
### Quality-filtering of reads
```
cutFs <- sort(list.files(path.cut, pattern = "_R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2_001.fastq.gz", full.names = TRUE))

filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, 
                     maxN=0, truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=T,minLen=190 ) 
```
### Error-rate estimation
```	
## calculating errors - on a subset of data to save time (https://benjjneb.github.io/dada2/bigdata.html)
set.seed(100)
errF <- learnErrors(filtFs, multithread=TRUE, nbases = 1e8, randomize=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE, nbases = 1e8, randomize=TRUE)
```
### Sample inference
```	
dadaFs <- dada(filtFs, err=errF, multithread=T)
dadaRs <- dada(filtRs, err=errR, multithread=T)
```
### Merging parired reads
```	
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs)
```
### Creating sequence table
```	 
seqtab <- makeSequenceTable(mergers)
```
### Removing chimeras
```	 
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
 #saving seqtab to file
fwrite( seqtab.nochim,  file=paste0(RUN_NAME,"_seqtab.nochim.csv"),sep=",",row.names=TRUE)
```
### Checking distribution of merged sequences' lengths
```
sums <- apply(seqtab.nochim,2, sum)
lengths <- nchar(names(sums))
lenght_distribution <- data.table( len=lengths, sum=sums) %>% group_by(len) %>% summarise( sum_of_given_len=sum(sum)) %>% data.frame()

lenght_distribution[order(lenght_distribution$sum_of_given_len),]
```

### Abundance filtering, length filtering  
 
proper length for 18S: **~350-380**  
This step removes lots of erronneous sequences and saves a lot of RAM to the next steps. According to the COI database, all tardigrade sequences have the same length (it was also checked on the actual metabarcoding data).
```	 
# if you want you can remove rare variants to save memory for the next steps 
seqtab.nochim.filtered_length.filtered_min_abundance<-seqtab.nochim.filtered_length[ , apply(seqtab.nochim.filtered_length,2,sum) >= 10]   #filtering reads winth counts < 10 
# sorting by abundance
seqtab.nochim.filtered_length.filtered_min_abundance<-seqtab.nochim.filtered_length.filtered_min_abundance[, order(-colSums(seqtab.nochim.filtered_length.filtered_min_abundance))]  

rm(seqtab)
rm(seqtab.nochim)
```

 
### Taxonomy classification with dada2 for 18S
```	 
 sequences <- colnames(seqtab.nochim.filtered_length.filtered_min_abundance)
 TAXA <- assignTaxonomy(sequences,DATABASE_PATH,minBoot=80,tryRC=TRUE, multithread = TRUE,outputBootstraps=TRUE)
  TAXA_DF<-data.frame(TAXA)
  seqtab_with_taxonomy<-cbind( TAXA_DF, t(seqtab.nochim.filtered_length.filtered_min_abundance.filtered_translation))
  fwrite( seqtab_with_taxonomy, file=paste0(RUN_NAME,"_zOTU_TABLE.csv"),sep=",",row.names=TRUE)  
```
### Exporting FASTA file with zOTU sequences
```	 
names<-paste0("zotu",1:ncol(seqtab.nochim.filtered_length.filtered_min_abundance.filtered_translation))
sequences<-colnames(seqtab.nochim.filtered_length.filtered_min_abundance.filtered_translation)
ZOTUS_FILENAME<-paste0("zotus_",RUN_NAME,"_",MARKER_NAME,".fasta")
write.fasta(as.list(sequences),names,paste0(RUN_NAME,"_zOTUS.fasta"))
```

