library(dplyr)
library(minfi)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(GEOquery)


cwd <-"/Users/Asus/Documents/GSEs/GSE97362/idat"
basedir <- "/Users/Asus/Documents/GSEs/GSE97362"

file_name <- "97362_sample_sheet.csv"
full_path <- file.path(basedir, file_name)
gse <-getGEO(filename = "/Users/Asus/Documents/GSEs/GSE97362/GSE97362_series_matrix.txt" )
pd <- pData(gse)

filenames <- list.files(cwd)
chunked_names <- lapply(strsplit(filenames, "_"),function (x) x[1:3])
chunked_names <-as.data.frame(chunked_names)
chunked_names <-t(chunked_names)
chunked_names <-as.data.frame(chunked_names)

colnames(chunked_names)<-c("Basename","Sentrix_ID","Sentrix_Position")
rownames(chunked_names) <-NULL
char_values <- as.character(pd$geo_accession)
matching <- match(char_values, chunked_names$Basename)

sampled <- chunked_names[matching,]
sample <- cbind(pd, sampled)
View(sample)

colnames(sample)[colnames(sample)=="geo_accession"] <-"Sample_ID"

write.csv(sample,file =full_path,row.names = FALSE)

rm(gse, basedir, char_values, file_name, cwd, matching, full_path, sample, sampled, pd, chunked_names)
