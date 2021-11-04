### format genes_dataset
### INPUT: hgnc_symbol log2fc padj (can be changed)
### OUTPUT: hgnc_symbol hgnc_ID chr start end strand log2fc padj (can be changed)

### load BioMart. Load dataset of genes of interest. Assign position info of genes of interest. Remove genes without relevant info

library("biomaRt") #load bioMart, from where extract gene info

### select database and dataset from BioMart
# listMarts() #list of databases
ensembl = useEnsembl("genes", GRCh = 37) #select database hg19
# ensembl = useEnsembl("ENSEMBL_MART_ENSEMBL") #select database hg38
# datasets = listDatasets(ensembl) #list of datasets
ensembl = useDataset("hsapiens_gene_ensembl", mart = ensembl) #select dataset of human and update variable

### load genes_datasets
genes_dataset = read.table("//ATLAS/IIER_TumoresInfantiles/RESULTADOS LABORATORIO/PROYECTO LIPI/4-Datasets/FEZF1/genes_dif_reg_EF1.txt", header = TRUE, dec = ",") #load file by path. CHANGE PATH
names(genes_dataset) = tolower(names(genes_dataset)) #headers in lowercase
# genes_dataset$log2fc = as.numeric(genes_dataset$log2fc) #column as numbers. change column name if needed
# genes_dataset$padj = as.numeric(genes_dataset$padj) #column as numbers. change column name if needed

### getBM(). three possible arguments:  attributes (output), filters and values for filters
# genes_info = getBM(attributes = c("hgnc_symbol", "hgnc_id", "chromosome_name", "start_position", "end_position", "strand"), mart = ensembl) #info of all human genes
genes_info = getBM(attributes = c("hgnc_symbol", "hgnc_id", "chromosome_name", "start_position", "end_position", "strand"), filters = "hgnc_symbol", values = genes_dataset, mart = ensembl) #info of a list of genes of interest. Other info columns can be selected

### merge genes_dataset with needed info
# genes_dataset = genes_info #when all genes are extracted
genes_dataset_backup = genes_dataset #backup just in case
# genes_dataset = genes_dataset_backup #recovery of dataset, keep commented
genes_dataset = merge(x = genes_dataset, y = genes_info, by.x = "hgnc_symbol", by.y = "hgnc_symbol", all.x = FALSE) #add info to genes dataset
genes_dataset = rename(genes_dataset, "chr" = "chromosome_name", "start" = "start_position", "end" = "end_position") #change col names for next steps
genes_dataset$start = as.numeric(genes_dataset$start) #column as numbers.
genes_dataset$end= as.numeric(genes_dataset$end) #column as numbers.

### remove non valid chr
row_to_keep = TRUE
deleted_rows = genes_dataset[FALSE,]
counter1 = 1
for (counter2 in 1:dim(genes_dataset)[1]){ 
  if (genes_dataset$chr[counter2] %in% c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, "X", "Y")) {
    row_to_keep[counter2] = TRUE
  } else {
    row_to_keep[counter2] = FALSE
    deleted_rows[counter1,] = genes_dataset[counter2,]
    counter1 = counter1 + 1
  }
}
genes_dataset = genes_dataset[row_to_keep,]
write.table(deleted_rows,"results/tables/deleted_rows.txt",sep="\t",row.names=FALSE, quote=FALSE) #create folder results/tables/ previously
cat(nrow(deleted_rows),"rows deleted. Check file at results/tables/deleted_rows\n") #to check if an interesting gene has been removed

### remove rows with blank symbol
row_to_keep = TRUE
deleted_rows = genes_dataset[FALSE,]
counter3 = 1
for (counter4 in 1:dim(genes_dataset)[1]){ 
  if (genes_dataset$hgnc_symbol[counter4] != "") {
    row_to_keep[counter4] = TRUE
  } else {
    row_to_keep[counter4] = FALSE
    deleted_rows[counter3,] = genes_dataset[counter4,]
    counter3 = counter3 + 1
  }
}
genes_dataset = genes_dataset[row_to_keep,]

### determine TSS (BioMart gives TSS of each transcript, we need the further one)
TSS = vector() 
for (counter in 1:dim(genes_dataset)[1]) { #wanders each gene, if sense, TSS is start; if antisense, TSS is end
  gene = genes_dataset[counter,]
  if (genes_dataset[counter,]$strand == "1") {
    new_var = genes_dataset[counter,]$start
  } else {
    new_var = genes_dataset[counter,]$end
  }
  TSS[counter] = new_var
}
genes_dataset$TSS = TSS

names(genes_dataset) = tolower(names(genes_dataset))
print(head(genes_dataset))
