############################################################
# Nombre del Script: metabarcoding_bioinfo_Gerald.R
# Descripción: "Metabarcoding 16S"
# Autor: Gerald Moreno Morales
# Fecha: 07/07/2024
# Referencias:
#"https://benjjneb.github.io/dada2/tutorial.html"
#"https://www.bioconductor.org/packages/3.3/bioc/manuals/dada2/man/dada2.pdf"
#"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4927377/"


############################################################


rm(list =ls())
setwd("/home/gerald/Documentos/maestria/metabarcoding/metagenomica/archivos_ubigem")
getwd()


################################################## 
#    Flujograma de Metabarcoding por denoising   #
##################################################
##Paquetes necesarios
###Instalar BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(BiocManager); packageVersion ("BiocManager")
####Instalar Dada2
BiocManager::install("dada2", force= TRUE)
library(dada2); packageVersion("dada2")
##Ejemplo solo con 1 par de reads
#Listar archivos Fastq con los que vot a trabajar
getwd()
path <- "/home/gerald/Documentos/maestria/metabarcoding/metagenomica/archivos_ubigem/fasqs"
list.files(path)

##Ahora leemos los nombres de los archivos fastq y realizamos algunas manipulaciones de cadenas para obtener listas coincidentes de los archivos fastq directos e inversos.
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq andSAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
print(fnFs)
print(fnRs)

# Extraer los nombres de las muestras
sample.namesF <- sapply(strsplit(basename(fnFs), "_"), function(x) paste(x[1], x[2], x[3], "R1", sep="_"))
sample.namesR <- sapply(strsplit(basename(fnRs), "_"), function(x) paste(x[1], x[2], x[3], "R2", sep="_"))
print("Forward sample names:")
print(sample.namesF)
print("Reverse sample names:")
print(sample.namesR)

###Nota: Consideraciones para sus propios datos: Es posible que sea necesario modificar las manipulaciones de cadenas si el formato de su nombre de archivo es diferente.
##Empezamos visualizando los perfiles de calidad de las lecturas directas:
plotQualityProfile(fnFs)
##Ahora visualizamos el perfil de calidad de las lecturas inversas:
plotQualityProfile(fnRs)

#Asigne los nombres de archivo para los archivos fastq.gz filtrados.
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.namesF, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.namesR, "_R_filt.fastq.gz"))
names(filtFs) <- sample.namesF
names(filtRs) <- sample.namesR

# Asegúrate de que los nombres de las muestras estén correctamente asignados
print(filtFs)
print(filtRs)


#En la región V3-V4 del gen 16S rRNA, la longitud típica es de 
#aproximadamente 460 pares de bases (pb). Aquí hay algunos pasos y recomendaciones 
#para ajustar tu código:
#Estoy cortando, los forwards en 200 y los reverse en 150:
#La longitud total después de la truncación debe ser aproximadamente 
#igual a la longitud del amplicón (460 pb) más la longitud de la superposición deseada. 
#Por ejemplo, si decides tener una superposición de 30 pb, las longitudes truncadas de 
#las secuencias forward y reverse deben sumar aproximadamente 490 pb 
#(460 pb del amplicón + 30 pb de superposición).


# On Linux set multithread=TRUE, en Windows multithread=False
# Aquí se asume truncLen de forward = 250 y reverse = 240 para obtener una buena superposición
truncLen_forward <- 250
truncLen_reverse <- 240
print(fnFs)
print(fnRs)
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(truncLen_forward, truncLen_reverse),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)

#otra manera
fastqPairedFilter
print(out)

###Aprender errores de Forwards y Reverse
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

##Siempre vale la pena visualizar las tasas de error estimadas:
plotErrors(errF, nominalQ=TRUE)



