rm(list=ls()) #remove all variables
cat("\014") #clean window

### create folders results/tables/ and results/browser/
### modify paths in scripts to use your own data
### execute each script in the proposed order

### load libraries
library(tidyverse)
library(valr)
library(dplyr)
library(writexl)

### format databases
# load  formatted gene data: symbol, id, coordinates, strand, FC...
source("scripts/format_genes.R")
genes_dataset_all = genes_dataset

# load formated GGAA repeats data: ID, chr, start, end...
source("scripts/format_repeats.R")
tandem_repeats = all_ef1_bs
names_for_columns = c("left_repeat_ID", "left_overlapping_ggaa_msat", "left_repeat_distance", "left_start", "left_end", "right_repeat_ID", "right_overlapping_ggaa_msat", "right_repeat_distance", "right_start", "right_exp_end")
source("scripts/genes_vs_repeats.R")

# source("scripts/format_repeats.R")
# tandem_repeats = all_ef1_bs[grepl('^orthCScore|^riggiCScommon', all_ef1_bs$id),]
# names_for_columns = c("left_exp_repeat_ID", "left_exp_repeat_distance", "left_exp_start", "left_exp_end", "right_exp_repeat_ID", "right_exp_repeat_distance", "right_exp_start", "right_exp_end")
# source("scripts/genes_vs_repeats.R")
# 
# tandem_repeats = all_ef1_bs[grepl('^RM|^mSat', all_ef1_bs$id),]
# names_for_columns = c("left_inf_repeat_ID", "left_inf_repeat_distance", "left_inf_start", "left_inf_end", "right_inf_repeat_ID", "right_inf_repeat_distance", "right_inf_start", "right_inf_end")
# source("scripts/genes_vs_repeats.R")

### determine the closer GGAA, weather informatic or experimental
closest_left_repeat_ID = vector()
closest_left_distance = vector()
closest_right_repeat_ID = vector()
closest_right_distance = vector() #create vectors to store data
for (counter6 in 1:dim(genes_dataset_all)[1]) {
  row_evaluated = genes_dataset_all[counter6,]
  if (is.na(row_evaluated$left_exp_repeat_distance)) { # if a value is NA and the script crashes
    row_evaluated$left_exp_repeat_distance = 1000000000000
  } 
  if (is.na(row_evaluated$left_inf_repeat_distance)) {
    row_evaluated$left_inf_repeat_distance = 1000000000000
  }
  if (is.na(row_evaluated$right_exp_repeat_distance)) { # if a value is NA and the script crashes
    row_evaluated$right_exp_repeat_distance = 1000000000000
  } 
  if (is.na(row_evaluated$right_inf_repeat_distance)) {
    row_evaluated$right_inf_repeat_distance = 1000000000000
  }
  
  if (row_evaluated$left_exp_repeat_distance <  row_evaluated$left_inf_repeat_distance) {
    closest_left_repeat_ID[counter6] =  row_evaluated$left_exp_repeat_ID
    closest_left_distance[counter6] =  row_evaluated$left_exp_repeat_distance
  } else {
    closest_left_repeat_ID[counter6] =  row_evaluated$left_inf_repeat_ID
    closest_left_distance[counter6] =  row_evaluated$left_inf_repeat_distance
  }
  if (row_evaluated$right_exp_repeat_distance <  row_evaluated$right_inf_repeat_distance) {
    closest_right_repeat_ID[counter6] =  row_evaluated$right_exp_repeat_ID
    closest_right_distance[counter6] =  row_evaluated$right_exp_repeat_distance
  } else {
    closest_right_repeat_ID[counter6] =  row_evaluated$right_inf_repeat_ID
    closest_right_distance[counter6] =  row_evaluated$right_inf_repeat_distance
  }
}

genes_dataset_all$closest_left_repeat_ID = closest_left_repeat_ID
genes_dataset_all$closest_left_distance = as.integer(closest_left_distance)
genes_dataset_all$closest_right_repeat_ID = closest_right_repeat_ID
genes_dataset_all$closest_right_distance = as.integer(closest_right_distance)

write_xlsx(genes_dataset_all,"results/tables/genes_near_promoter_TF_up_AsiEF.xlsx") #store in a txt file. CHANGE FILE TO STORE OUTPUT
print(head(genes_dataset_all))

closest_distance = vector()
closest_repeat_ID = vector()
for (counter7 in 1:dim(genes_dataset_all)[1]) {
  row_evaluated_2 = genes_dataset_all[counter7,]
  
  if (is.na(row_evaluated_2$closest_left_distance)) { # if a value is NA and the script crashes
    row_evaluated_2$closest_left_distance = 1000000000000
  } 
  if (is.na(row_evaluated_2$closest_right_distance)) {
    row_evaluated_2$closest_right_distance = 1000000000000
  }

  if (row_evaluated_2$closest_left_distance < row_evaluated_2$closest_right_distance) {
    closest_repeat_ID[counter7] = row_evaluated_2$closest_left_repeat_ID
        closest_distance[counter7] = row_evaluated_2$closest_left_distance
  } else {
    closest_repeat_ID[counter7] = row_evaluated_2$closest_right_repeat_ID
    closest_distance[counter7] = row_evaluated_2$closest_right_distance
  }
}

genes_dataset_all$closest_repeat_ID = closest_repeat_ID
genes_dataset_all$closest_distance = as.integer(closest_distance)

print(head(genes_dataset_all))

write_xlsx(genes_dataset_all,"results/tables/genes_near_promoter_genes_dif_reg_EF1.xlsx") #store in a txt file. CHANGE FILE TO STORE OUTPUT
