## Comparison of the pairwise distances in the 16s gene sequences with respect to the ones calculated on hyper variable parts of the gene


### Introduction

For this project I will be performing an analysis on the 16S rRNA gene, which is about 1500 bp and highly conserved among different species of bacteria. The conserved portions are separated by nine variable regions (called V regions) located in its sequence. These variable regions, numbered from 1 to 9, allow to discriminate bacterial taxa from one another, to tell different bacterial strains apart, and to help asses the genetic diversity in environmental samples. It’s frequently used in phylogenetic studies and microbial taxonomy.

Since I can use 16S to tell different bacteria apart, I want to see if it’s possible to just use a small part of the entirety of the gene (one of the aforementioned hyper variable V regions) to also be able to tell bacteria apart. If that were to be true, it would mean that just analyzing one of the V regions would be enough to assess with a good degree of certainty from which bacteria that region comes from, thus cutting computational times since the region to analyze would be just a small fraction of the entire 16S gene.

The workflow is then set as follows: first a distance matrix must be computed for each full 16S aligned sequence pair. Then, after selecting the hyper variable regions of our choice and extracting them from the gene sequences, a distance matrix must be computed for each region aligned sequence pair. If the regression model built using the two matrices is linear and after seeing if the hyper variable regions perfectly correlate to the 16S full length gene sequence, we can suggest that just using the hyper variable region’s sequence instead of the full 16S sequence may be sufficient to discriminate bacteria from each other in taxonomy.

### Full lenght sequences

I used a FASTA file containing sequences of the 16S gene of many different bacteria. The file, called “input_sequences”, was aligned using MAFFT (Multiple Alignment using Fast Fourier Transform).
The MAFFT output FASTA file was uploaded on R. I calculated the pairwise distances for each of the sequences in the aligned file and put the information into a distance matrix. Each cell of the distance matrix contains the squared root of the pairwise distance of each pair of sequences: for this reason, to get the actual percentage identity values, I squared the matrix, converted the distances into similarities and multiplied by 100 to get the percentages. I then replaced the similarity values of sequences compared with themselves with “NA”.

```r
alignments <- read.alignment("/Users/chiaraalbertini/Desktop/brilli/aligned_sequences.fasta")
pairwise_distances <- dist.alignment(alignments)
distance_matrix <- as.matrix(pairwise_distances)
distance_matrix_squared <- distance_matrixˆ2
similarity_scores <- (1 - distance_matrix_squared) * 100
diag(similarity_scores) <- NA
preview <- similarity_scores[1:3, 1:3]
print(preview)
```

### V3-V4 sequences

The next step is to find the coordinates of the hyper variable regions of interest: the first one I analyzed is the region spanning over the regions V3 and V4. I found the forward and reverse primers for the V3-V4 region in the paper “Comparison of Two 16S rRNA Primers (V3–V4 and V4–V5) for Studies of Arctic Microbial Communities” (https://doi.org/10.3389/fmicb.2021.637526) and with those I could identify and “cut” the V3-V4 region out of the sequences.

```r
file <- "/Users/chiaraalbertini/Desktop/brilli/input_sequences.fasta"
all_sequences <- readDNAStringSet(file, format = "fasta")
fw_primer <- DNAString("CCTACGGGNGGCWGCAG")
rv_primer <-DNAString("GACTACHVGGGTATCTAATCC")
fw_matches <- vmatchPattern(fw_primer, all_sequences,max.mismatch=2)
rev_matches <- vmatchPattern(reverseComplement(rv_primer),all_sequences,max.mismatch=2)
starts <- start(fw_matches) %>% as.numeric()
ends <- end(rev_matches) %>% as.numeric()
valid_indices <- !is.na(starts) & !is.na(ends)
V3_V4_ranges <- IRanges(start=starts[valid_indices], end=ends[valid_indices])
#I check how many sequences I've lost
length(ends) - length(V3_V4_ranges)
```

```r
#now I extract the sequences
all_sequences_filtered <- all_sequences[valid_indices]
V3_V4_sequences <- subseq(all_sequences_filtered, V3_V4_ranges)
V3_V4_sequences
```
I exported the output file and aligned it using MAFFT.

```r
output_v3_v4 <- "/Users/chiaraalbertini/Desktop/brilli/V3-V4.fasta"
writeXStringSet(V3_V4_sequences, filepath = output_v3_v4, format = "fasta")
```
I uploaded the aligned file on R and calculated the pairwise distances; then I built the distance matrix as I did before (this time the file only contained the V3-V4 region sequences).

```r
alignments_v3_v4 <- read.alignment("/Users/chiaraalbertini/Desktop/brilli/aligned_V3-V4.fasta")
pairwise_distances_v3_v4 <- dist.alignment(alignments_v3_v4)
distance_matrix_v3_v4 <- as.matrix(pairwise_distances_v3_v4)
distance_matrix_v3_v4_squared <- distance_matrix_v3_v4ˆ2
similarity_scores_v3_v4 <- (1 - distance_matrix_v3_v4_squared) * 100
diag(similarity_scores_v3_v4) <- NA
```

```r
preview34 <- similarity_scores_v3_v4[1:3, 1:3]
print(preview34)
```

I built a regression model to see if the comparison of the pairwise distances in the full 16S gene and the ones
calculated only on the V3-V4 region was linear.

```r
ds <- similarity_scores[valid_indices,valid_indices]
dist_vector_full <- ds[lower.tri(ds)]
dist_vector_v3_v4 <- similarity_scores_v3_v4[lower.tri(similarity_scores_v3_v4)]
regression <- lm(dist_vector_v3_v4 ~ dist_vector_full)
summary(regression)
```

The high R-squared value (0.9203), the significant p-value and the beta coefficient of 1.039 suggests that the linear model fits the data extremely well (even more impressive when considering the big dataset we are working with) and that it exists an almost proportional relationship between the similarity values calculated on the full 16S sequences and the ones calculated on the V3-V4 region sequences. This is not enough: to be certain that I can use hyper variable sequences in place of the full 16s gene when performing bacterial taxonomy, I must be sure that the hyper variable sequence are in a perfect proportional relationship with the full length sequences. If that was not the case, I may end up with hyper variable sequences that leads to incorrectly classification of bacteria. This can be checked in two ways: with a linear hypothesis test and by looking at the Log2 ratio of the full 16S similarity / V3-V4 similarity.

```r
log2_ratios34 <- log2((ds) / (similarity_scores_v3_v4))
log2_ratios_vector34 <- as.vector(log2_ratios34)[lower.tri(log2_ratios34)]
hist(log2_ratios_vector34,
  main = "Log2 ratios (16s similarity / V3-V4 similarity)",
  xlab = "Log2 ratio",
  ylab = "Frequency",
  col = "blue",
  breaks = 100)
```
The histogram of the Log2 full 16S similarity / V3-V4 similarity is, as expected, showing a Gaussian distribution; but the distribution is not perfectly centered on 1, it’s slightly shifted toward a negative value. The mean value of the Log2 ratios vector is very close to 1, but once again to see if it’s significantly different from 1 or not I have to perform a linear hypothesis test.

```r
linearHypothesis(regression, "dist_vector_full = 1")
```
The very high F statistic value combined with the significant p-value show that, as we saw before, there is indeed a strong correlation between the identities of the two groups, but that relationship isn’t perfectly proportional x=y. This means that there is a possibility that two 16S sequences which have more than 97% identity may have their identity value drop below that threshold if we only compare their respective V3-V4 regions, leading to the belief that the two sequences belong to two distinct bacteria.

### V7-V9 regions
We have seen this behavior between the pairwise distances of the full 16S sequences and the ones calculated on the hyper variable region V3-V4. I did the same analysis picking a different region to see if results were consistent. This time I wasn’t able to find primers for only V8 and V9 regions; however, the paper “Primer, Pipelines, Parameters: Issues in 16S rRNA Gene Sequencing” (doi: 10.1128/mSphere.01202-20) had primer sequences for the region spanning from V7 to V9. I therefore selected these primers to proceed with the analysis.

```r
file <- "/Users/chiaraalbertini/Desktop/brilli/input_sequences.fasta"
all_sequences <- readDNAStringSet(file, format = "fasta")

fw_primer <- DNAStringSet("CAACGAGCGCAACCCT")
rv_primer <- DNAString("TACGGYTACCTTGTTACGACTT")

fw_matches <- vmatchPattern("CAACGAGCGCAACCCT", all_sequences, max.mismatch = 2)
rev_matches <- vmatchPattern(reverseComplement(rv_primer), all_sequences, max.mismatch = 2)

# Extracting the first start and end positions
starts <- sapply(fw_matches, function(x) if(length(x) > 0) start(x)[1] else NA_integer_)
ends <- sapply(rev_matches, function(x) if(length(x) > 0) end(x)[1] else NA_integer_)
# Filter valid indices where both starts and ends are not NA
valid_indices <- !is.na(starts) & !is.na(ends)
V7_V9_ranges <- IRanges(start = starts[valid_indices], end = ends[valid_indices])
# V7_V9_ranges now contains the ranges for valid sequences

#I check how many sequences I've lost
length(ends) - length(V7_V9_ranges)
```

```r
#now I extract the sequences
all_sequences_filtered <- all_sequences[valid_indices]
V7_V9_sequences <- subseq(all_sequences_filtered, V7_V9_ranges)
V7_V9_sequences
```

```r
output_v7v9 <- '/Users/chiaraalbertini/Desktop/brilli/V7-V9.fasta'
writeXStringSet(V7_V9_sequences, filepath = output_v7v9, format = "fasta")
```

I exported the output and I aligned it with MAFFT. I uploaded the MAFFT output in R and calculated the pairwise identity values just as before.

```r
alignments_v7_v9 <- read.alignment("/Users/chiaraalbertini/Desktop/brilli/aligned_V7-V9.fasta")
pairwise_distances_v7_v9 <- dist.alignment(alignments_v7_v9)
distance_matrix_v7_v9 <- as.matrix(pairwise_distances_v7_v9)
distance_matrix_v7_v9_squared <- distance_matrix_v7_v9ˆ2
similarity_scores_v7_v9 <- (1 - distance_matrix_v7_v9_squared) * 100
diag(similarity_scores_v7_v9) <- NA
```

```r
preview79 <- similarity_scores_v7_v9[1:3, 1:3]
print(preview79)
```

I then built the regression model in the same way as before.

```r
ds1 <- similarity_scores[valid_indices,valid_indices]

dist_vector_full <- ds1[lower.tri(ds1)]
dist_vector_v7_v9 <- similarity_scores_v7_v9[lower.tri(similarity_scores_v7_v9)]

regression_model_v7_v9 <- lm(dist_vector_v7_v9 ~ dist_vector_full)
summary(regression_model_v7_v9)
```

As expected, and exactly as before, the R-squared value is very high, the p-value is significant and beta is close to 1; once again I have to check is beta is significantly different from 1.

```r
log2_ratios79 <- log2((ds1) / (similarity_scores_v7_v9))
log2_ratios_vector79 <- as.vector(log2_ratios79)[lower.tri(log2_ratios79)]
hist(log2_ratios_vector79,
  main = "Log2 ratios (16s identities / V7-V9 identities)",
  xlab = "Log2 ratio",
  ylab = "Frequency",
  col = "blue",
  breaks = 100)
```

Let’s now compute the Log2 full 16S pairwise distances / V7-V9 pairwise distances. Once again, the peak of the distribution does not fall on 0, meaning that the beta coefficient is significantly different from 1.

```r
linearHypothesis(regression_model_v7_v9, "dist_vector_full = 1")
```
Just like before, I get a very high F statistic value and a significant p-value, meaning that a correlation between the two groups have been detected but not enough to be able to safely use yet the hyper variable V7-V9 regions as proxy for the entire 16S gene.

### Establishing if the hyper variable regions can be used as a proxy
Now the first thing I have associate the sequence IDs to the species name because IDs are unique and don’t match with the species they come from.

```r
all_lines <- readLines("/Users/chiaraalbertini/Desktop/brilli_files/input_sequences.fasta")

headers_df <- grep("ˆ>", all_lines, value = TRUE) %>%
  data.frame(strings = .) %>%
  mutate(strings = str_remove(strings, "ˆ>"))

headers_df <- headers_df %>%
  separate(strings, into = c("ID", "species"), sep = " ", extra = "merge") %>%
  mutate(species = str_extract(species, "ˆ\\w+\\s+\\w+|ˆ[ˆ\\s]+\\s+[ˆ\\s]+\\s+[ˆ\\s]+"))
#I extract first 2 words -> genus and species
```

I want to create a three column table that has, for each row, the sequence IDs pair and the corresponding similarity value. I will then filter this table to only keep the rows where the sequences have a similarity value equal or higher than 97 (in general, 97% of sequence identity is the threshold used to say that two sequences belong to the same bacteria species). I then replace the IDs with the relative bacteria’s name. In a perfect world, the resulting table after all of these steps should contain only sequences pairs (whose IDs belong to the same bacteria) with value >= 97.

```r
long_format_similarity_scores_df <- melt(similarity_scores, varnames = c("ID1", "ID2"))

#filter sequences with value >= 97
long_format_similarity_scores_df_over97 <- long_format_similarity_scores_df %>%
  filter(value >=97)

#replace with real species name
sp1_sp2 <- long_format_similarity_scores_df_over97 %>%
  left_join(headers_df, by = c("ID1" = "ID"))%>%
  dplyr::rename(species1 = species) %>%
  left_join(headers_df, by = c("ID2" = "ID")) %>%
  dplyr::rename(species2 = species) %>%
  select(species1,species2)

exact_match_proportion <- mean(sp1_sp2$species1 == sp1_sp2$species2)
print(exact_match_proportion*100)
```

But since we don’t live in a perfect world, the table will contain also situations where two IDs have a similarity value over 97% while belonging to different bacteria. I calculate exact_match_proportion, which is the percentage of the sequences in the table that actually belong to the same species.
I the did the same thing for the V3-V4 sequences, to check if the exact_match_proportion value is different if I use only the sequences of the regions.

```r
long_format_similarity_scores_df_V3_V4 <- melt(similarity_scores_v3_v4, varnames = c("ID1", "ID2"))
long_format_similarity_scores_df_V3_V4_over97 <- long_format_similarity_scores_df_V3_V4 %>%
  filter(value >=97)

sp1_sp2_V3_V4 <- long_format_similarity_scores_df_V3_V4_over97 %>%
  left_join(headers_df, by = c("ID1" = "ID"))%>%
  dplyr::rename(species1 = species) %>%
  left_join(headers_df, by = c("ID2" = "ID")) %>%
  dplyr::rename(species2 = species) %>%
  select(species1,species2)

exact_match_proportion_V3_V4 <- mean(sp1_sp2_V3_V4$species1 == sp1_sp2_V3_V4$species2)
print(exact_match_proportion_V3_V4*100)
```

This time, exact_match_proportion returns a higher value than the one calculated using the full sequences. This means that, in this case, using the V3-V4 sequences allows me to classify as belonging to the same species a more correct number of sequences; therefore, I’d be inclined to use the hyper variable regions to tell different bacterial strains apart because it seems to be more precise.
A threshold correction here should not be necessary, since it would be bringing my value towards the one of the full sequences, which in this case is less ideal since I’ve already seen that using the V3-V4 region leads to better results than using the full length sequences.
A new threshold may be set using the linear regression formula with the corresponding intercept and beta parameters. After that, I check the corrected exact_match_proportion and I am expecting a lower value than before.

```r
new_tresh <- (97 - unname(regression$coefficients[1]))/unname(regression$coefficients[2])
long_format_similarity_scores_df_V3_V4_over97_corr <- long_format_similarity_scores_df_V3_V4 %>% filter(value >=97)

sp1_sp2_V3_V4_corr <- long_format_similarity_scores_df_V3_V4_over97_corr %>%
  left_join(headers_df, by = c("ID1" = "ID"))%>%
  dplyr::rename(species1 = species) %>%
  left_join(headers_df, by = c("ID2" = "ID")) %>%
  dplyr::rename(species2 = species) %>%
  select(species1,species2)

exact_match_proportion_V3_V4_corr <- mean(sp1_sp2_V3_V4_corr$species1 == sp1_sp2_V3_V4_corr$species2)
print(exact_match_proportion_V3_V4_corr*100)
```

Now let’s do the same for V7-V9.

```r
long_format_similarity_scores_df_V7_V9 <- melt(similarity_scores_v7_v9, varnames = c("ID1", "ID2"))
long_format_similarity_scores_df_V7_V9_over97 <- long_format_similarity_scores_df_V7_V9 %>% filter(value >=97)

sp1_sp2_V7_V9 <- long_format_similarity_scores_df_V7_V9_over97 %>%
  left_join(headers_df, by = c("ID1" = "ID"))%>%
  dplyr::rename(species1 = species) %>%
  left_join(headers_df, by = c("ID2" = "ID")) %>%
  dplyr::rename(species2 = species) %>%
  select(species1,species2) %>% filter(species1!="16S ribosomal",species2!="16S ribosomal")

exact_match_proportion_V7_V9 <- mean(sp1_sp2_V7_V9$species1 == sp1_sp2_V7_V9$species2)
print(exact_match_proportion_V7_V9*100)
```

```r
new_tresh <- (97 - unname(regression_model_v7_v9$coefficients[1]))/unname(regression_model_v7_v9$coeffic

long_format_similarity_scores_df_V7_V9_over97_corr <- long_format_similarity_scores_df_V7_V9 %>% filter(value >= new_tresh)

sp1_sp2_V7_V9_corr <- long_format_similarity_scores_df_V7_V9_over97_corr %>%
  left_join(headers_df, by = c("ID1" = "ID"))%>%
  dplyr::rename(species1 = species) %>%
  left_join(headers_df, by = c("ID2" = "ID")) %>%
  dplyr::rename(species2 = species) %>%
  select(species1,species2)

exact_match_proportion_V7_V9_corr <- mean(sp1_sp2_V7_V9_corr$species1 == sp1_sp2_V7_V9_corr$species2)
print(exact_match_proportion_V7_V9_corr*100)
```

The exact_match_proportion is very low (27%) when compared to the exact_match_proportion value calculated on the full length sequences, and after the correction it drops even lower than it originally was (24%). The reason may be found in the fact that lowering the threshold allows for more noise to be considered. In a case such as this one, I would not trust the V7-V9 sequences as a proxy and I would just use the full 16S gene sequences to do bacterial taxonomy.

### Conclusions

The pairwise distances calculated on all the 16S sequence pairs correlate with the pairwise distances calculated on both V3-V4 regions sequence pairs and V7-V9 sequence pairs; however, correlation by itself isn’t enough to affirm that just analyzing one of the hyper variable V regions found in the 16S gene may be enough differentiate bacteria species and strains, since I proved that the beta parameter in both the linear regression models I built is significantly different from 1, meaning that I can’t choose to just pick an hyper variable region as a proxy for the full 16S gene sequence without diving a little bit deeper.

In order to decide if it would be possible to use just the hyper variable sequences in place of the full 16S gene sequences to perform bacterial taxonomy, a preliminary step must be first performed. In this analysis I showed how to calculate a percentage value that I used as the parameter to decide if it was better, depending on the hyper variable region analyzed each time, to use the hyper variable region sequences or the full length ones. I suggest to calculate the percentage value for both the full length sequences and the sequences of the hyper variable regions of choice and to compare them:
- in the cases where the value is higher when calculated on the hyper variable sequences (like in the V3-V4 sequences case), I would use the hyper variable regions as a proxy for bacterial taxonomy. In a case such as this one we’d get the benefit of doing the taxonomy analysis on just a portion of the entire 16S gene, thus greatly cutting computational time and costs.
- in the cases where the value is higher when calculated on the full length sequences I’d suggest it would be best to either fix a new threshold or, like in the analyzed V7-V9 case, to just use the full length sequences.
