### Wanders list of genes; for each genes, wanders list of repeats
all_left_repeat_IDs = vector()
all_left_overlapping = vector()
all_left_distances = vector()
all_left_start = vector()
all_left_end = vector()
all_right_repeat_IDs = vector()
all_right_overlapping = vector()
all_right_distances = vector() #create vectors to store data of repeats and distance
all_right_start = vector()
all_right_end = vector()

for (counter5 in 1:dim(genes_dataset_all)[1]) { #wanders each row of genes of interest
  current_gene = genes_dataset_all[counter5,]$hgnc_symbol
  current_chr = genes_dataset_all[counter5,]$chr #select only repeats of the same chromosome
  current_TSS = genes_dataset_all[counter5,]$tss
  
  distance_left = 1000000000000 #initialize with a high value, to look for the minimum
  distance_right = 1000000000000
  selected_left_repeat_ID = "" #to store the repeat
  selected_right_repeat_ID = ""
  
  repeats_current_chr = tandem_repeats[tandem_repeats$chr == current_chr,] #subset of repeats of the chromosome of the gene
  if (dim(repeats_current_chr)[1] == 0) { #when a chromosome does not have repeats
    all_left_repeat_IDs[counter5] = NA
    all_left_overlapping[counter5] = NA
    all_left_distances[counter5] = NA
    all_left_start[counter5] = NA
    all_left_end[counter5] = NA    
    all_right_repeat_IDs[counter5] = NA
    all_right_overlapping[counter5] = NA
    all_right_distances[counter5] = NA
    all_right_start[counter5] = NA
    all_right_end[counter5] = NA
    next
  }
  for (counter2 in 1:dim(repeats_current_chr)[1]) { #wanders each entry of repeats
    current_repeat_ID = repeats_current_chr[counter2,]$id #repeat to evaluate
    current_repeat_position = repeats_current_chr[counter2,]$average_position #position of repeat to evaluate
    if (current_repeat_position < current_TSS) { #if repeat position is lower than TSS, it is on the left
      posible_distance = current_TSS - current_repeat_position #distance repeat to TSS
      if (posible_distance < distance_left) { #if distance is lower than the stored distance, change stored ID and distance. If not, stays the same
        distance_left = posible_distance
        selected_left_repeat_ID = current_repeat_ID
      } else {
        # print("not changed")
      }
    } else { #if position of repeat is higher than TSS, it is on the right
      posible_distance = current_repeat_position - current_TSS
      if (posible_distance < distance_right) {
        distance_right = posible_distance
        selected_right_repeat_ID = current_repeat_ID
      } else {
        # print("not changed")
      }
    }
  }
  # print(paste(selected_left_repeat_ID, distance_left, selected_right_repeat_ID, distance_right))
  if (distance_left == 1000000000000) {
    all_left_repeat_IDs[counter5] = NA 
    all_left_distances[counter5] = NA
    all_left_overlapping[counter5] = NA
    all_left_start[counter5] = NA
    all_left_end[counter5] = NA
  } else {
    all_left_repeat_IDs[counter5] = selected_left_repeat_ID
    all_left_distances[counter5] = distance_left
    all_left_overlapping[counter5] = repeats_current_chr$overlapping_ggaa_msat[repeats_current_chr$id == selected_left_repeat_ID]
    all_left_start[counter5] = repeats_current_chr$start[repeats_current_chr$id == selected_left_repeat_ID]
    all_left_end[counter5] = repeats_current_chr$end[repeats_current_chr$id == selected_left_repeat_ID]
  }
  if (distance_right == 1000000000000) {
    all_right_repeat_IDs[counter5] = NA
    all_right_distances[counter5] = NA
    all_right_overlapping[counter5] = NA
    all_right_start[counter5] = NA
    all_right_end[counter5] = NA
  } else {
    all_right_repeat_IDs[counter5] = selected_right_repeat_ID
    all_right_distances[counter5] = distance_right
    all_right_overlapping[counter5] = repeats_current_chr$overlapping_ggaa_msat[repeats_current_chr$id == selected_right_repeat_ID]
    all_right_start[counter5] = repeats_current_chr$start[repeats_current_chr$id == selected_right_repeat_ID]
    all_right_end[counter5] = repeats_current_chr$end[repeats_current_chr$id == selected_right_repeat_ID]
  }
}
genes_dataset_all[names_for_columns[1]] = all_left_repeat_IDs
genes_dataset_all[names_for_columns[2]] = all_left_overlapping
genes_dataset_all[names_for_columns[3]] = all_left_distances
genes_dataset_all[names_for_columns[4]] = all_left_start
genes_dataset_all[names_for_columns[5]] = all_left_end
genes_dataset_all[names_for_columns[6]] = all_right_repeat_IDs
genes_dataset_all[names_for_columns[7]] = all_right_overlapping
genes_dataset_all[names_for_columns[8]] = all_right_distances
genes_dataset_all[names_for_columns[9]] = all_right_start
genes_dataset_all[names_for_columns[10]] = all_right_end

print(head(genes_dataset_all))
