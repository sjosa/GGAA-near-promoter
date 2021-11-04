### Load repeats dataset, format each dataset to take chr, start, end and give them an individual ID. Each dataset requires individual formatting

### install library GenomicRanges
library(valr)

### load GGAA repeats data. Format data to chr, start, end, ID. CHANGE PATH TO YOUR REPEATS DATASET
# repeatmasker = read.delim("//ATLAS/IIER_TumoresInfantiles/RESULTADOS LABORATORIO/PROYECTO LIPI/4-Datasets/Tandem Repeats from USCS Table Browser - 2021.abr.15/Repeats-RepeatMasker-GRCh37.txt", header= TRUE, sep = "\t", dec = ".") #load repeats from RepeatMasker
# repeatmasker_GGAA = repeatmasker[repeatmasker$repName == "(GGAA)n" | repeatmasker$repName == "(TTCC)n",] #take only those with GGAA in sense (GGAA)n or antisense (TTCC)n. Other combinations GGAA where in strange chr
# GGAA_mSat = read.delim("//ATLAS/IIER_TumoresInfantiles/RESULTADOS LABORATORIO/PROYECTO LIPI/4-Datasets/Orth et al_ESCLA/Orth_et_al_Suppl1_GGAA-mSat_A673.txt") #A673 (GGAA)n-mSat from Orth et al, 2021
orth_chipseq_core = read.delim("//ATLAS/IIER_TumoresInfantiles/RESULTADOS LABORATORIO/PROYECTO LIPI/4-Datasets/Orth et al_ESCLA/Orth_et_al_Suppl2_Core-ChIPseq_ONLY_mSat.txt", header= TRUE, sep = "\t", dec = ".") #load localization of Binding Sites for EF by ChIPseq from Orth et al, 2021
# riggi_chipseq_common = read.delim("//ATLAS/IIER_TumoresInfantiles/RESULTADOS LABORATORIO/PROYECTO LIPI/4-Datasets/Riggi et al/NIHMS639914-supplement-1_common.txt") #A673/SKNMC FLI1 ChIPseq from Riggi et al, 2014

### format Repeatmasker: chr, start, end, id
# names(repeatmasker_GGAA) = tolower(names(repeatmasker_GGAA))
# repeatmasker_GGAA = rename(repeatmasker_GGAA, "chr" = "genoname", "start" = "genostart", "end" = "genoend")
# repeatmasker_GGAA = repeatmasker_GGAA[c("chr", "start", "end")] #use only those columns
# repeatmasker_GGAA$id = paste("RM_",0001:nrow(repeatmasker_GGAA), sep = "") #give ID to each repetition
# cat("track name=\"RepeatMasker GGAA\"\n", file = "results/browser/repeatmasker_GGAA_for_browser.txt") #header for GenomeBrowser
# write.table(repeatmasker_GGAA[c("chr", "start", "end", "id")], file = "results/browser/repeatmasker_GGAA_for_browser.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE) #store in txt for GenomeBrowser
# repeatmasker_GGAA$chr = substring(repeatmasker_GGAA$chr, 4) #remove word chr

### format GGAA_mSat: chr, start, end, id
# names(GGAA_mSat) = tolower(names(GGAA_mSat))
# cat("track name=\"GGAA-mSat\"\n", file = "results/browser/GGAA_mSat_for_browser.txt") #header for GenomeBrowser
# write.table(GGAA_mSat[c("chr", "start", "end", "id")], file = "results/browser/GGAA_mSat_for_browser.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE) #store in txt for GenomeBrowser
# GGAA_mSat$chr = substring(GGAA_mSat$chr, 4) #remove word chr

### format orth_chipseq_core: chr, start, end, id
names(orth_chipseq_core) = tolower(names(orth_chipseq_core)) #headings in lowercase
orth_chipseq_core$id = paste("orthCScore_",0001:nrow(orth_chipseq_core), sep = "") #give ID to each repetition
cat("track name=\"Core GGAA Orth ChIPseq\"\n", file = "results/browser/orth_chipseq_core_for_browser.txt") #header for GenomeBrowser
write.table(orth_chipseq_core[c("chr", "start", "end", "id")], file = "results/browser/orth_chipseq_core_for_browser.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE) #store in txt for GenomeBrowser
orth_chipseq_core$chr = substring(orth_chipseq_core$chr, 4) #remove word chr

### format riggi_chipseq_common: chr, start, end, id
# names(riggi_chipseq_common) = tolower(names(riggi_chipseq_common)) #headings in lowercase
# riggi_chipseq_common$id = paste("riggiCScommon_",0001:nrow(riggi_chipseq_common), sep = "") #give ID to each repetition
# cat("track name=\"Common GGAA Riggi ChIPseq\"\n", file = "results/browser/riggi_chipseq_common_for_browser.txt") #header for GenomeBrowser
# write.table(riggi_chipseq_common[c("chr", "start", "end", "id")], file = "results/browser/riggi_chipseq_common_for_browser.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE) #store in txt for GenomeBrowser
# riggi_chipseq_common$chr = substring(riggi_chipseq_common$chr, 4) #remove word chr

# all_ef1_bs = rbind(repeatmasker_GGAA, GGAA_mSat, orth_chipseq_core, riggi_chipseq_common) #join all GGAA cordinates in a unique variable
all_ef1_bs = orth_chipseq_core

### remove non valid chr
row_to_keep = TRUE
deleted_rows = all_ef1_bs[FALSE,]
counter1 = 1
for (counter2 in 1:dim(all_ef1_bs)[1]){ 
  if (all_ef1_bs$chr[counter2] %in% c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, "X", "Y")) {
    row_to_keep[counter2] = TRUE
  } else {
    row_to_keep[counter2] = FALSE
    deleted_rows[counter1,] = all_ef1_bs[counter2,]
    counter1 = counter1 + 1
  }
}
all_ef1_bs = all_ef1_bs[row_to_keep,]

cat("track name=\"All EF1 BS\"\n", file = "results/browser/all_ef1_bs.txt") #header for GenomeBrowser
write.table(all_ef1_bs[c("chr", "start", "end", "id")], file = "results/browser/all_ef1_bs.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE) #store in txt for GenomeBrowser

### calculate middle position of binding sites
average_position = vector() 
for (element in 1:dim(all_ef1_bs)[1]) {
  average_position[element] = round(mean(c(all_ef1_bs[element,]$start, all_ef1_bs[element,]$end)))
}
all_ef1_bs$average_position = average_position

write.table(all_ef1_bs,"results/tables/all_ef1_bs.txt",sep="\t",row.names=FALSE, quote=FALSE) #output with all repeats
print(head(all_ef1_bs))