##############################################################################
# Nombre del Script: metabarcoding_bioinfo_Gerald.R
# Descripción: "Metabarcoding 16S"
# Autor: Gerald Moreno Morales
# Fecha: 07/07/2024
# Referencias:
#"https://benjjneb.github.io/dada2/tutorial.html"
#"https://www.bioconductor.org/packages/3.3/bioc/manuals/dada2/man/dada2.pdf"
#"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4927377/"
##############################################################################

###Instalar BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(BiocManager); packageVersion ("BiocManager")
####Instalar Dada2
BiocManager::install("dada2", force= TRUE)
library(dada2); packageVersion("dada2")
###Como verifico cual es mi path???
getwd()
path <- "C:/bioinfo/MiSeq_SOP"
list.files(path)
#Ahora leemos los nombres de los archivos fastq y realizamos algunas manipulaciones de cadenas para obtener listas coincidentes de los archivos fastq directos e inversos.
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq andSAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
###Nota: Consideraciones para sus propios datos: Es posible que sea necesario modificar las manipulaciones de cadenas si el formato de su nombre de archivo es diferente.
##Empezamos visualizando los perfiles de calidad de las lecturas directas:
plotQualityProfile(fnFs[5:6])
##Ahora visualizamos el perfil de calidad de las lecturas inversas:
plotQualityProfile(fnRs[1:4])

#Asigne los nombres de archivo para los archivos fastq.gz filtrados.
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# On Windows set multithread=FALSE
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE)
head(out)

##Aprender errores de Forwards y Reverse
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
##Siempre vale la pena visualizar las tasas de error estimadas:
plotErrors(errF, nominalQ=TRUE)

#####Ahora estamos listos para aplicar el algoritmo de inferencia de muestra central a los datos de secuencia filtrados y recortados.
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
## Inspeccionando el dada-class objeto de vuelto:
dadaFs[[2]]

###Unir forward y reverse
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

##Ahora podemos construir una tabla de variantes de secuencia de amplicones (ASV), una versión de mayor resolución de la tabla OTU producida por métodos tradicionales.
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

### Inspeccionar la distribución de longitudes de secuencia.
table(nchar(getSequences(seqtab)))

###Remover chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus",
                                    multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
###Como verificación final de nuestro progreso, veremos la cantidad de lecturas que superaron cada paso del proceso:
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers,
                                                                       getN), rowSums(seqtab.nochim))
# Si procesa una sola muestra, elimine las llamadas de sapply: p. reemplazar sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged",
                     "nonchim")
rownames(track) <- sample.names
head(track)

###Guardar nuestro avance, aquí solo van los archivos generados anteriormente. 
save(seqtab.nochim, mergers, dadaRs,dadaFs, file = "microbioma.RData")
load("microbioma.RData")


#NOTA: DESCARGAR BASE DE DATOS DE SILVA: (https://zenodo.org/records/4587955)
taxa <-
  assignTaxonomy(seqtab.nochim,"C:/bioinfo/MiSeq_SOP/silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread=TRUE)
#Inspeccionemos las asignaciones taxonómicas:
  # Removing sequence rownames for display only
  taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)
      
##Instalar phyloseq
BiocManager::install("phyloseq")
library(phyloseq); packageVersion("phyloseq")
##Instalar Biostrings
BiocManager::install("Biostrings")
library(Biostrings); packageVersion("Biostrings")
###Instalar ggplot2
library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw())
BiocManager::install("GenomeInfoDb")
library(GenomeInfoDb)
BiocManager::install("DECIPHER")
library(DECIPHER)
library(phangorn)
library(dada2)
####Crear la metada y el objeto phylose
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)


day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out

###Otra opción, realizar un metadata en drive o excel y descargarlo como csv. El archivo csv es metadata.csv
df2 <- read.csv(file = "metadata2.csv", row.names = 1)


###Introducir arbol filogenético usando NJ
#Extraer secuencias de Dada2 output
seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs
###Correr el alineamiento de secuencias usando DECIPHER
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)
##Change sequence alignment output into a phyDat structure
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
###Create distance matrix
dm <- dist.ml(phangAlign)
##Perform Neighbor joining
treeNJ <- NJ(dm) 
###Internal maximum likelihood
fit = pml(treeNJ, data=phangAlign)
##Este código se utiliza para ajustar un modelo de máxima verosimilitud para filogenias, específicamente un modelo GTR (General Time Reversible).
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
# Import phyloseq and Remove mock sample
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
               sample_data(df2),
               tax_table(taxa), phy_tree(fitGTR$tree))
ps2 <- prune_samples(sample_names(ps) != "Mock", ps)

ps2

# Computamos prevalencia para cada feature y la guardamos en un data frame
prevdf = apply(X = otu_table(ps2),
               MARGIN = ifelse(taxa_are_rows(ps2), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Le agregamos la taxonomía
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps2),
                    tax_table(ps2))

# Seleccionamos las taxa de interés
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps2, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps2),color=Phylum)) +
  # Agregamos una línea para nuestro umbral
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

# Definimos el umbral de prevalencia a un 5%
(prevalenceThreshold = 0.05 * nsamples(ps2))

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
(ps3 = prune_taxa(keepTaxa, ps2))

# Reemplazamos las secuencias por un nombre genérico
taxa_names(ps3) <- paste0("ASV", seq(ntaxa(ps3)))




