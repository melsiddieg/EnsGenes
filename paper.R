library(data.table)
library(tidyverse)
library(janitor)
library(UpSetR)
library(gprofiler2)
mrt37 <- fread('data/grch37.tsv')
mrt37 <- janitor::clean_names(mrt37)
setindexv(mrt37,c("gene_stable_id","gene_name"))
mrt37 <- clean_names(mrt37)
mrt38 <- fread('data/grch38.txt.gz')
mrt38 <- janitor::clean_names(mrt38)
setindexv(mrt38,c("gene_stable_id","gene_name"))
mrt37 %>% filter(gene_type=='protein_coding')  %>%select(gene_name) %>% unique() %>% nrow()
mrt38 %>% filter(gene_type=='protein_coding')  %>%select(gene_name) %>% unique() %>% nrow()

sharedgenes <- inner_join(mrt37,mrt38,by='gene_stable_id')
# grch37 will have x suffix and grch38 will have y suffix
names(sharedgenes)
##
####
# fwrite(DiscrepantGenes,'Results/discrpant_genes.tsv', sep = '\t')
#### find all genes that changed names and biotype between assemblies
mrtype <- sharedgenes %>%filter(gene_type.y=='protein_coding'&gene_type.x!='protein_coding') %>% 
  group_by(gene_stable_id) %>% 
  summarise(number=n())
## locate genes that changes between assembly
idx2 <- match(mrtype$gene_stable_id,sharedgenes$gene_stable_id)
## final set of genes
ChangedGenes <- sharedgenes[idx2,]
setDT(ChangedGenes)
## basic stats
### genes that changed names
ChangedGenes[gene_name.x!=gene_name.y,.N]
### how many phenotypes in the grch37 subset
ChangedGenes[nchar(phenotype_description.x)>1,.N]
### how many phenotypes in the grch38 subset
ChangedGenes[nchar(phenotype_description.y)>1,.N]
# ### with phenotypes
# pheno <- ChangedGenes[nchar(phenotype_description.y)>1,]
# phenoTable <- pheno %>% select(gene_stable_id,gene_name.x,gene_name.y, phenotype_description.x,phenotype_description.y, ref_seq_match_transcript)
# fwrite(phenoTable,"phenotypes2.tsv", sep = "\t")
### with refseq
# ChangedGenes[nchar(ref_seq_match_transcript)>1,.N]
# ChangedGenes[,`:=`(inRefSeq=nchar(ref_seq_match_transcript)>1,havePhenotype=nchar(phenotype_description.y)>1, inClinGen=nchar(source_name.y)>1)]
# ChangedGenes[gene_name.x!=gene_name.y,haveNewHGNC:=TRUE]
# ChangedGenes[,haveNewBiotype:=TRUE]
# final <- ChangedGenes %>% select(gene_stable_id,gene_type.x,starts_with("have"), starts_with("in"))
# names(final)[2] <- "Grch37BioType"
# final2 <- final %>% mutate_if(is_logical, as.numeric)
# final3 <- final2 %>% mutate_all(~replace(., is.na(.), 0))
###
# upset(final3,sets = c("havePhenotype","haveNewHGNC","haveNewBiotype","inRefSeq","inClinGen"), order.by =  "freq")


### double checking refseq
#ref <- fread("Results/gProfiler_hsapiens_8-8-2020_9-11-46 PM.csv")
ref <- gconvert(query = unique(ChangedGenes$gene_stable_id), organism = "hsapiens", target="REFSEQ_MRNA")
names(ref)[c(2,4)] <- c("gene_stable_id","converted_alias")
names(ref)
setDT(ref)
###  sanity check
ref[,length(unique(gene_stable_id))]
###
ref2 <- ref %>% group_nest(gene_stable_id)
pluckRefSeq <- function(data){
  if(data$converted_alias=="None"){
    NA
  }else{
    res = data$converted_alias[[1]]
    res
  }
}

ref2$refseq <- ref2$data %>% map_chr(pluckRefSeq)
setDT(ref2)
ref2[,data:=NULL]
###
ft <- ChangedGenes %>% left_join(ref2)
setDT(ft)
ft[,`:=`(inRefSeq=!is.na("refseq"),havePhenotype=nchar(phenotype_description.y)>1, inClinGen=nchar(source_name.y)>1)]
ft[gene_name.x!=gene_name.y,haveNewHGNC:=TRUE]
ft[,haveNewBiotype:=TRUE]
final <- ft %>% select(gene_stable_id,gene_type.x,starts_with("have"), starts_with("in"))
names(final)[2] <- "Grch37BioType"
final2 <- final %>% mutate_if(is_logical, as.numeric)
final3 <- final2 %>% mutate_all(~replace(., is.na(.), 0))
upset(final3,sets = c("havePhenotype","haveNewHGNC","haveNewBiotype","inRefSeq","inClinGen"), order.by =  "freq")
####
## basic stats
### genes that changed names
ft[gene_name.x!=gene_name.y,.N]
### how many phenotypes in the grch37 subset
ft[nchar(phenotype_description.x)>1,.N]
### how many phenotypes in the grch38 subset
ft[nchar(phenotype_description.y)>1,.N]
### how many have refseq transcripts
ft[nchar(refseq)>1,.N]
###
pheno <- ft[nchar(phenotype_description.y)>1,]
phenoTable <- pheno %>% select(gene_stable_id,gene_name.x,gene_name.y, phenotype_description.x,phenotype_description.y, refseq)
fwrite(phenoTable,"phenotypes2.tsv", sep = "\t")
