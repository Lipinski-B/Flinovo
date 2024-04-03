library(Seurat)
library(readr)
library(readxl)
library(Matrix)
library(ggplot2)
library(patchwork)
library(stringr)

# Download the list of dissociation genes from the Genome biol paper 
diss.genes.GenBiol <- read_csv("/Users/Manon/Documents/Manon/Bioinfo_scRNAseq/Scripts_local/Laurenti_gene dissoc/ListegènedissociationGenome biol.csv")

# Rename by gene name in "gene symbol" column
diss.genes.GenBiol <- diss.genes.GenBiol$gene_symbol

# Load Seurat Object

diss<- list(diss.genes.GenBiol)
DefaultAssay(singlet) <- "RNA"
singlet = AddModuleScore(singlet, features = diss, name = 'stress_response_score')

Idents(singlet) <- "Condition"
new <- subset(singlet, subset = Condition == c('collagenase + actinomycin D', 'collagenase', 'subtilisin'))

VlnPlot(object = new, features = 'stress_response_score1', group.by = "Condition") +
  ggtitle('Dissociation signature') +
  theme(axis.title.x = element_blank()) +
  theme(text=element_text(size=17, family="Arial")) +
  ggplot2::scale_fill_manual(values = c('#86BBD8', '#00B159','#9CB3F6'))


