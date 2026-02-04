# 2024 nov
# dge in reproductive tissues
## libraries ----
library(tidyverse)
library(edgeR)
library(igraph)
library(VennDiagram)
library(gridExtra)
library(scales)

## import data& merge ----
# 1 - htseqcount data
F_Go_1029<- read.table("msc/firststep/data/htseqcount/Tps_F_Go_Ad_18-1029.txt", stringsAsFactors=F,header=F,sep="\t")
F_Go_1029 <- F_Go_1029 %>% rename(Gene = V1, Count = V2) %>% 
  mutate(sample = "F_Go_1029")

F_Go_1030 <- read.table("msc/firststep/data/htseqcount/Tps_F_Go_Ad_18-1030.txt", stringsAsFactors=F,header=F,sep="\t")
F_Go_1030 <- F_Go_1030 %>% rename(Gene = V1, Count = V2) %>% 
  mutate(sample = "F_Go_1030" )

F_Go_1031 <- read.table("msc/firststep/data/htseqcount/Tps_F_Go_Ad_18-1031.txt", stringsAsFactors=F,header=F,sep="\t")
F_Go_1031 <- F_Go_1031 %>% rename(Gene = V1, Count = V2) %>% 
  mutate(sample = "F_Go_1031" )

F_Go_1041 <- read.table("msc/firststep/data/htseqcount/Tps_F_Go_Ad_18-1041.txt", stringsAsFactors=F,header=F,sep="\t")
F_Go_1041 <- F_Go_1041 %>% rename(Gene = V1, Count = V2) %>% 
  mutate(sample = "F_Go_1041")

M_Te_1033 <- read.table("msc/firststep/data/htseqcount/Tps_M_Te_Ad_18-1033.txt", stringsAsFactors=F,header=F,sep="\t")
M_Te_1033 <- M_Te_1033 %>% rename(Gene = V1, Count = V2) %>% 
  mutate(sample ="M_Te_1033")

M_Te_1035 <- read.table("msc/firststep/data/htseqcount/Tps_M_Te_Ad_18-1035.txt", stringsAsFactors=F,header=F,sep="\t")
M_Te_1035 <- M_Te_1035 %>% rename(Gene = V1, Count = V2) %>% 
  mutate(sample = "M_Te_1035")

M_Te_1036 <- read.table("msc/firststep/data/htseqcount/Tps_M_Te_Ad_18-1036.txt", stringsAsFactors=F,header=F,sep="\t")
M_Te_1036 <- M_Te_1036 %>% rename(Gene = V1, Count = V2) %>% 
  mutate(sample = "M_Te_1036")

M_Te_1042 <- read.table("msc/firststep/data/htseqcount/Tps_M_Te_Ad_18-1042.txt", stringsAsFactors=F,header=F,sep="\t")
M_Te_1042 <- M_Te_1042 %>% rename(Gene = V1, Count = V2) %>% 
  mutate(sample = "M_Te_1042")

# combine df-s
list_of_dfs <- list(F_Go_1029, F_Go_1030,F_Go_1031,F_Go_1041,M_Te_1033,M_Te_1035,M_Te_1036, M_Te_1042)
all_counts <- list_of_dfs %>%
  bind_rows() %>%  
  select(Gene, Count, sample) %>%  #only relevant columns
  pivot_wider(names_from = sample, values_from = Count)  


# 2 - dupgenefinder
transposed <- read.table("msc/firststep/data/tps_dupgenefinder/Tps.transposed.pairs.txt", stringsAsFactors=F,header=T,sep="\t")
transposed <- transposed %>%  mutate( duplication_type = "transposed") %>% 
  rename(Duplicate.1 = Transposed, Duplicate.2 = Parental ) # need to remember this!

proximal <- read.table("msc/firststep/data/tps_dupgenefinder/Tps.proximal.pairs.txt", stringsAsFactors=F,header=T,sep="\t")
proximal <-  mutate(proximal, duplication_type = "proximal")

wgd <- read.table("msc/firststep/data/tps_dupgenefinder/Tps.wgd.pairs.txt", stringsAsFactors=F,header=T,sep="\t")
wgd <-  mutate(wgd, duplication_type = "wgd")
nrow(wgd)

tandem <- read.table("msc/firststep/data/tps_dupgenefinder/Tps.tandem.pairs.txt", stringsAsFactors=F,header=T,sep="\t")
tandem <-  mutate(tandem, duplication_type = "tandem")

dispersed <- read.table("msc/firststep/data/tps_dupgenefinder/Tps.dispersed.pairs.txt", stringsAsFactors=F,header=T,sep="\t")
dispersed <- mutate(dispersed, duplication_type = "other")

duplicate_pairs <-  bind_rows(transposed, proximal, wgd, tandem, dispersed)

duplicate_pairs$Duplicate.1 <- sub("\\..*$", "", duplicate_pairs$Duplicate.1)
duplicate_pairs$Duplicate.2 <- sub("\\..*$", "", duplicate_pairs$Duplicate.2)
# I hope these weren't important

str(duplicate_pairs)

# exclude wgd pairs 
duplicate_pairs <- filter(duplicate_pairs, duplication_type != "wgd")

write.table(duplicate_pairs, file="msc/firststep/data/duplicate_pairs.txt",quote=F, sep="\t")


## general data on paralogs in tps----
duplicate_proportions <- duplicate_pairs %>%
  group_by(duplication_type) %>%
  summarise(percentage = 100* n() / nrow(duplicate_pairs))

(duplicates_piechart <- ggplot(duplicate_proportions, aes(x = "", y = percentage, fill = duplication_type)) +
    geom_bar(stat = "identity", width = 1) +  # Create a bar chart
    coord_polar(theta = "y") +                # Convert to a pie chart
    labs(fill = "Duplication mechanism") +         # Add legend title
    theme_void()    +                        # Clean up the chart
    theme(plot.background = element_rect(fill = "white", color = NA),
          legend.title = element_text(size = 14),  # Adjust legend title size
          legend.text = element_text(size = 14),    # Adjust legend text size
          legend.key.size = unit(0.7, "cm"))+        # Adjust size of legend keys (color boxes)
    geom_text(aes(label = paste0(round(percentage, 1), "%")), 
              position = position_stack(vjust = 0.5),
              size = 5) + # Add percentage labels
    scale_fill_brewer(palette = "Paired"))    
ggsave("msc/firststep/plots/tps_duplication_types_proportions.png", plot = duplicates_piechart, bg = "white", width = 8, height = 6, dpi = 500, )


## filtering by gene length----

# Create edgeR object which will contain count data and normalization factors
expr <- DGEList(counts=all_counts)

# Read in gene length information and append to data
gene_length <- read.table("msc/firststep/data/Tps_gene_length.txt",stringsAsFactors=F)
# V1 - gene name, V2 - Length

expressed_genes <- all_counts$Gene

gene_length <- subset(gene_length, V1 %in% expressed_genes)
gene_length <- gene_length[match(expr$genes$Gene, gene_length$V1), ]# Is this ok?
gene_length <- na.omit(gene_length) # somehow there were na-s idk

gene_length_vector <- c(gene_length$V2)

all(gene_length$Gene == rownames(expr)) #should print TRUE - ok

valid_genes <- !is.na(match(expr$genes$Gene, gene_length$V1))
# Filter the expr object to include only valid genes i.e. those w a gene length
expr <- expr[valid_genes, ]

# verify
all(expr$genes$Gene == gene_length$V1)  # Should return TRUE - ok
anyNA(gene_length)                      # Should return FALSE - ok

# Calculate normalization factors that we need to account for
norm_expr <- calcNormFactors(expr)

# Convert counts to RPKM
rpkm_norm <- rpkm(norm_expr, log=FALSE,gene.length=gene_length_vector)

# convert the RPKM result into a data frame
rpkm_norm_df <- as.data.frame(rpkm_norm)
# add the gene names as a new column 
rpkm_norm_df$Gene <- expr$genes$Gene  # from expr$genes$Gene

write.table(rpkm_norm_df, file="msc/firststep/rpkm.txt",quote=F, sep="\t")

## filter lowly expressed genes. Keep genes that have min.2 RPKM expression in at least half of the samples-----

# Set the last column (Gene) as rownames
rownames(rpkm_norm_df) <- rpkm_norm_df$Gene
rpkm_norm_df <- rpkm_norm_df[, -ncol(rpkm_norm_df)]  

# set threshold for filtering (2 RPKM)
threshold <- 2

# no of samples
num_samples <- ncol(rpkm_norm_df)

# how many samples must have RPKM >= min_threshold
min_samples_required <- ceiling(num_samples / 2)

# filter genes based on the threshold condition
filtered_genes <- apply(rpkm_norm_df, 1, function(x) sum(x >= threshold) >= min_samples_required)

# subset the dataframe to include only the filtered genes
rpkm_filtered <- rpkm_norm_df[filtered_genes, ]

# add back the gene column
rpkm_filtered$Gene <- rownames(rpkm_filtered)

# Data for DGE: raw counts, but only keep the rows that appear in rpkm_filtered"Gene"
counts_filtered <- as.data.frame(all_counts[all_counts$Gene %in% rpkm_filtered$Gene, ])
write.table(counts_filtered, file="msc/firststep/data/filtered_counts.txt",quote=F, sep="\t")

## table on whether genes are DE or not ----
# set gene column as rownames again
rownames(counts_filtered) <- counts_filtered$Gene
counts_filtered <- counts_filtered %>% select(-Gene)
# information on which group each sample belongs to
conditions <- factor(c("F","F","F","F","M","M","M","M"))

# Create edgeR object which will contain count data and normalization factors
expr <- DGEList(counts=counts_filtered,group=conditions)
# Calculate normalization factors that we need to account for
norm_expr <- calcNormFactors(expr)

# Models how variable gene expression is across replicates

# estimate common dispersion 1st
norm_expr <- estimateCommonDisp(norm_expr)
# estimate gene-specific dispersion
norm_expr <- estimateTagwiseDisp(norm_expr)
et <- exactTest(norm_expr)
# does this do the test on the right variable (tagwise disp)?

# Obtain p-values for between group comparisons
p <- et$table$PValue
# Adjust p-values to correct for multiple testing (more info on this here: https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/05a_hypothesis_testing.html)
p_adj <- p.adjust(p, method = c("fdr"), n = length(p))
# Append adjusted p value to output table
table <- et$table
# table$Padj <- p_FDR original code but this wasnt defined before so i assume its p_adj?
table$Padj <- p_adj
table$Gene  <- rownames(table) # gene column back again
write.table(table, file="msc/firststep/data/differential_expression_testes_gonad.txt",quote=F, sep=",")

## DE stats on all genes----

# % of all genes that are DE
diff_expr_genes <- table %>% filter (Padj <= 0.05)
nrow(diff_expr_genes)/nrow(table) # 54 %, sounds ok

# add sex-biased/not
table<- table %>% 
  mutate(diff_exp_sexes = ifelse (Padj <= 0.05, "sex-biased", "unbiased")) %>% 
  # add paralog status
  mutate(duplicate_status = ifelse(
    ((Gene %in%duplicate_pairs$Duplicate.1)|(Gene %in% duplicate_pairs$Duplicate.2)),"paralog","singleton"))
table_duplicates <-
  table %>% 
  group_by(duplicate_status) %>% 
  summarise(count = length(duplicate_status))
table_duplicates

table_edited$duplicate_status <- as.factor(table_edited$duplicate_status)
str(table_edited$duplicate_status)
table_edited <- table %>% 
  filter(!((duplicate_status == "paralog") & !((Gene %in% duplicates_comparable$Duplicate.1) |(Gene %in% duplicates_comparable$Duplicate.1))))
nrow(table_edited)

# Bias type 
table$sex_bias <- "unbiased"
# if log2Foldchange > 1 and p < 0.05, male biased
table$sex_bias[table$logFC > 1 & table$Padj < 0.05] <- "M"
# if log2Foldchange < -1 and p < 0.05, female biased
table$sex_bias[table$logFC < 1 & table$Padj < 0.05] <- "F"


sum_for_volcano <- table %>% 
  group_by(sex_bias) %>% 
  summarise(count = n()) %>% 
  mutate(percentage = round((count / nrow(table)) * 100, 1))
sum_for_volcano
# Volcano plot
(volcanoplot <- ggplot(table, aes(x=logFC, y = -log10(Padj), col = sex_bias))+
    geom_point(size = 0.7)+
    theme_classic()+
    ylim(c(-8, 60))+
    ggtitle("a) Gonad (n = 8365) \n")+
    theme(legend.position = "none",
          plot.title = element_text(size = 13),
          axis.title.x = element_blank())+
    xlim(c(-10,10))+
    labs(x = "\nLog2FoldChange", y = "-log10 P(adj)\n")+
    scale_colour_manual(values = c("firebrick", "dodgerblue", "gray60"))+
    geom_hline(yintercept=-log10(0.05), col="red")+
    annotate("label", x = 7, y = 30, 
             label = "28.3%", 
             color = "black", fill = "white", label.size = 0.5) +
    annotate("label", x = -7, y = 30, 
             label = "25.3%", 
             color = "black", fill = "white", label.size = 0.5) +
    annotate("label", x = 0, y = -5, 
             label = "46.4%", 
             color = "black", fill = "white", label.size = 0.5))
ggsave("msc/firststep/plots/tps_volcano.png", width = 8, height = 6, dpi = 600)

## Data on duplicates and their expression-----

# filter duplicates that are expressed
duplicates_comparable <- duplicate_pairs %>% filter (Duplicate.1 %in% table$Gene & Duplicate.2 %in% table$Gene)

# what proportion of duplicates are expressed in reproductive tissues
nrow(duplicates_comparable)/nrow(duplicate_pairs) # 29.7%

# what proportion of dups are DE
diff_expr_duplicates <- diff_expr_genes %>%
  filter(Gene %in% duplicates_comparable$Duplicate.1 | Gene %in% duplicates_comparable$Duplicate.2)
nrow(diff_expr_duplicates)/nrow(duplicates_comparable)
# 65 % are sex-biased

# what proportion of singletons are DE
singletons <- table %>% filter (duplicate_status == "singleton")
(singletons_summary  <- singletons %>% 
    group_by(diff_exp_sexes) %>% 
    summarise(count = length(diff_exp_sexes)))
2044/(2044+1970)

# add new column to dups dataframe,: is duplicate 1 sex-biased? is duplicate 2 sex-biased?
duplicates_comparable<- duplicates_comparable %>% mutate( 
  Dup1_de = ifelse(Duplicate.1 %in% diff_expr_genes$Gene, "yes", "no"),
  Dup2_de = ifelse (Duplicate.2 %in%diff_expr_genes$Gene, "yes", "no"))

# and specific pattern
duplicates_comparable <- duplicates_comparable %>% 
  mutate(pattern = case_when(
    (Dup1_de == "yes"& Dup2_de == "yes") ~ 1,
    (Dup1_de == "no"& Dup2_de == "no") ~ 0,
    ((Dup1_de == "yes"& Dup2_de == "no") | (Dup1_de == "no"& Dup2_de == "yes")) ~2))

# Add logfc values
logfc_to_add <- diff_expr_genes %>% select(logFC, Gene) %>% rename(Duplicate.1 = Gene)
duplicate_pairs_logfc <- left_join(duplicates_comparable, logfc_to_add, by = "Duplicate.1")
duplicate_pairs_logfc <- rename(duplicate_pairs_logfc, logfc_dup1 = logFC)
logfc_to_add <- rename(logfc_to_add, Duplicate.2 = Duplicate.1)
duplicate_pairs_logfc_both <- left_join(duplicate_pairs_logfc, logfc_to_add, by = "Duplicate.2")
duplicate_pairs_logfc_both <- rename(duplicate_pairs_logfc_both, logfc_dup2 = logFC)
duplicates_comparable <- duplicate_pairs_logfc_both

# adding expression pattern
duplicates_comparable <- duplicates_comparable %>% 
  mutate(
    exp_bias_dup1 = case_when(
      logfc_dup1 < 1 ~ "F",
      logfc_dup1 >1 ~ "M", 
      is.na(logfc_dup1) ~ "unbiased"),
    exp_bias_dup2 = case_when(
      logfc_dup2 < 1 ~ "F",
      logfc_dup2 >1 ~ "M", 
      is.na(logfc_dup2) ~ "unbiased"),
    bias_pattern = case_when(
      (exp_bias_dup1 == "unbiased" & exp_bias_dup2 =="unbiased") ~ "UU",
      ((exp_bias_dup1 =="unbiased" & exp_bias_dup2 == "F") |
         (exp_bias_dup1 == "F" & exp_bias_dup2 == "unbiased")) ~ "UF",
      ((exp_bias_dup1 =="unbiased" & exp_bias_dup2 == "M") |
         (exp_bias_dup1 == "M" & exp_bias_dup2 == "unbiased")) ~ "UM",
      (exp_bias_dup1 == "M" & exp_bias_dup2 =="M") ~ "MM",
      (exp_bias_dup1 == "F" & exp_bias_dup2 =="F") ~ "FF",
      ((exp_bias_dup1 =="M" & exp_bias_dup2 == "F") |
         (exp_bias_dup1 == "F" & exp_bias_dup2 == "M")) ~ "FM"))

# discordant vs concordant 
duplicates_comparable<- duplicates_comparable%>% 
  mutate(discordant = ifelse(exp_bias_dup1 == exp_bias_dup2, "concordant", "discordant"))
(discordant_summary <- duplicates_comparable %>% 
    group_by(discordant) %>% 
    summarise(count = length(discordant)))
1684/(1684+1355)
# 55% show discordant expression

# pattern in transposed genes
transposed_exp <- duplicates_comparable %>% filter(duplication_type == "transposed") %>% 
  filter(pattern == 2)
summary_transposed <- transposed_exp %>% 
  group_by(Dup1_de) %>% 
  summarise(count = length(Dup1_de))
summary_transposed
# when only 1 transposed copy is DE,
# transposed is DE = 181,
# parental is DE = 111


## plots on duplicates----

# 1 - proportions of duplication types
duplicate_proportions_expressed <- duplicates_comparable %>% # this is from the filtered counts -
  # are those the ones that are expressed, or just the ones we can compare? Confused.
  group_by(duplication_type) %>%
  summarise(percentage = 100* n() / nrow(duplicate_pairs))
# Plot proportions
(expressed_duplicates_barplot <- ggplot(duplicate_proportions_expressed, aes(x = duplication_type, y = percentage, fill = duplication_type)) +
    geom_bar(stat = "identity") +  
    labs(x = "\nDuplication type", y = "% of all pairs\n") +
    scale_fill_brewer(palette = "Paired") +
    theme_classic())
ggsave("msc/firststep/plots/gonad_duplication_types_proportions.png", dpi = 500)

# 2 - basic patterns in duplicates that are expressed
duplicates_comparable$pattern<- as.factor( duplicates_comparable$pattern)
duplicates_comparable$duplication_type<- as.factor(duplicates_comparable$duplication_type)
(basic_plot <- ggplot (duplicates_comparable, aes(x = pattern, fill = duplication_type))+
    geom_bar(stat = "count")+
    labs(x = "\nPattern", y = "Count\n")+
    guides(fill = guide_legend(title = "Dupication type"))+
    scale_x_discrete(labels = c("Both no", "Both yes", "One yes, one not")) +
    theme_classic())

# 3 - Specific patterns
# Calculate proportions for each bias_pattern and duplication_type
duplicates_percent <- duplicates_comparable %>%
  group_by(bias_pattern, duplication_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(percent = 100 * count / sum(count))  

(expression_patterns_plot <- ggplot(duplicates_percent, aes(x = bias_pattern, y = percent, fill = duplication_type)) +
    geom_bar(stat = "identity") +  # Use calculated percentages
    labs(x = "\nPattern", y = "Percentage (%)\n") +
    ggtitle("a) Gonad\n")+
    scale_fill_brewer(palette = "Paired") +
    guides(fill = guide_legend(title = "Duplication type")) +
    theme_classic()+
    theme(legend.position = "none",
          plot.title = element_text(size = 13),
          axis.title.y = element_text(size = 13),
          axis.text = element_text(size = 12),
          axis.title.x = element_blank()))

# pie chart
duplicates_for_pie_gonad <- duplicates_percent %>%
  group_by(bias_pattern) %>%
  summarise(percent = sum(percent))

# reorder the bias_pattern factor
duplicates_for_pie_gonad$bias_pattern <- factor(
  duplicates_for_pie_gonad$bias_pattern, 
  levels = c("UU", "FF", "MM", "UM", "UF", "FM"))

(pie_gonad <- ggplot(duplicates_for_pie_gonad, aes(x = "", y = percent, fill = bias_pattern))+
    geom_bar(stat = "identity", width = 1)+
    coord_polar(theta = "y")+
    labs (x = NULL, y = NULL, title = "a) Gonad")+
    geom_text(aes(label = paste0(bias_pattern, "\n", round(percent, 1), "%")), 
              position = position_stack(vjust = 0.5), 
              size = 4) +
    guides(fill = guide_legend(title = "Expression pattern")) +
    theme_classic()+
    scale_fill_viridis_d(option = "plasma")+
    # scale_fill_brewer(palette = "Paired") +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          plot.title = element_text (size = 13)))

ggsave("msc/firststep/plots/tps_gonad_expression_patterns_dups.png", dpi=600)

# Create the pie chart
pie(percent, 
    labels = paste0(bias_pattern, "\n", round(percent, 1), "%"), 
    col = colors, 
    main = "a) Gonad")

(duplication_types_plot <- ggplot(duplicates_percent, aes(x = duplication_type, y = percent, fill = bias_pattern)) +
    geom_bar(stat = "identity", position = "fill") +  # Use calculated percentages
    labs(x = "\nPattern", y = "Percentage (%)") +
    scale_fill_brewer(palette = "Paired") +
    guides(fill = guide_legend(title = "Duplication type")) +
    theme_classic())
ggsave("msc/firststep/plots/tps_gonad_duplication_types_expression.png", dpi=600)

# concordant vs discordant
for_discordant_plot <- duplicates_comparable %>% 
  filter(! duplication_type == "other")

for_discordant_plot_summary <- duplicates_comparable %>% 
  group_by (duplication_type) %>% 
  summarise(count = length(duplication_type))
for_discordant_plot_summary

(discordant_plot <- ggplot(duplicates_comparable, aes(x = duplication_type, fill = discordant))+
    geom_bar(position = "fill", stat = "count")+
    labs (x = "\nDuplication type", y = "Percentage\n")+
    scale_y_continuous(labels = scales::percent) +  
    #  guides(fill = guide_legend(title = "Expression pattern")) +
    scale_x_discrete(labels = c( "other\n (2188)", "proximal\n (119)", "tandem\n (117)", "transposed\n (615)" ))+
    scale_fill_brewer(palette = "Paired")+
    theme_classic()+
    ggtitle ("a) Gonad\n")+
    theme(legend.position = "none",
          axis.title = element_blank(),
          plot.title = element_text(size= 18),
          #  axis.title.y = element_text(size = 14),
          axis.text = element_text(size = 14)))
ggsave("msc/firststep/plots/tps_gonad_discordant_duptypes.png", dpi = 600, width = 6, height =6)

transposed<- duplicates_comparable %>% filter(duplication_type == "transposed")
(transposed_sum <- transposed %>% 
    group_by(discordant) %>% 
    summarise(count = length(discordant)))
342/(342+273)


## Chi squared tests-----
# Paralog/singleton, sex-biased/unbiased
# step 1 - add duplicate/not duplicate status to counts_filtered
counts_filtered$Gene <- rownames(counts_filtered)
counts_filtered <- counts_filtered %>% mutate(
  duplicate = ifelse((Gene %in% duplicate_pairs$Duplicate.1 | Gene %in% duplicate_pairs$Duplicate.2),
                     "duplicated", "not_duplicated"),
  de = ifelse ((Gene %in% diff_expr_genes$Gene), "sex-biased", "unbiased"))

chisq_data <- counts_filtered %>% select(duplicate,de)
chisq_matrix <- table(chisq_data$duplicate, chisq_data$de)
chisq.test(chisq_matrix)
chisq_matrix
chisq.test(chisq_matrix)$expected
chisq.test(chisq_matrix)$residuals

# Within sex-biased genes: dup/singleton, M/F bias
chisq_data_2<- table %>% 
  select(duplicate_status, sex_bias)
chisq_matrix_2 <- table(chisq_data_2$duplicate_status, chisq_data_2$sex_bias)
chisq.test(chisq_matrix_2)
chisq.test(chisq_matrix_2)$residuals
chisq_matrix_2
# no association/ dif in whether sex-biased genes are female/male biased and whether they are duplicates/not.
# unlike wyman

chisq_data_3<- table %>%
  select(duplicate_status, sex_bias)
chisq_matrix_3 <- table(chisq_data_3$duplicate_status, chisq_data_3$sex_bias)
chisq.test(chisq_matrix_3)
chisq.test(chisq_matrix_3)$expected
chisq.test(chisq_matrix_3)$residuals
chisq_matrix_3

# discordant vs concordant

chisq_data_4 <- duplicates_comparable %>% 
  select(duplication_type, discordant)
chisq_matrix_4 <- table(chisq_data_4$duplication_type, chisq_data_4$discordant)
chisq.test(chisq_matrix_4)
chisq.test(chisq_matrix_4)$residuals


for_chisq_6 <- duplicates_comparable %>% 
  select(duplication_type, bias_pattern)
chisq_matrix_6 <- table(for_chisq_6$duplication_type, for_chisq_6$bias_pattern)
chisq.test(chisq_matrix_6)
chisq.test(chisq_matrix_6)$residuals

for_chisq_7 <- table %>% 
  select(duplicate_status, sex_bias)
chisq_matrix_7 <- table(for_chisq_7$duplicate_status, for_chisq_7$sex_bias)
chisq.test(chisq_matrix_7)
chisq.test(chisq_matrix_7)$residuals

length(which(duplicates_comparable$duplication_type == "transposed"))

# discordant vs concordant, based on tissue type
brain_pairs_chisq <- duplicates_comparable_brain %>% 
  select(Duplicate.1, Duplicate.2, duplication_type, discordant, bias_pattern) %>% 
  mutate(tissue = "brain")
gonad_pairs_chisq <- duplicates_comparable %>% 
  select(Duplicate.1, Duplicate.2, duplication_type, discordant, bias_pattern) %>% 
  mutate(tissue = "gonad")
antenna_pairs_chisq <-duplicates_comparable_antenna %>% 
  select(Duplicate.1, Duplicate.2, duplication_type, discordant, bias_pattern) %>% 
  mutate(tissue = "antenna")
gut_pairs_chisq <- duplicates_comparable_gut %>% 
  select(Duplicate.1, Duplicate.2, duplication_type, discordant, bias_pattern) %>% 
  mutate(tissue = "gut")

all_pairs_tissues <- bind_rows(gonad_pairs_chisq, antenna_pairs_chisq, brain_pairs_chisq, gut_pairs_chisq)
chisq_data_5 <- all_pairs_tissues %>% 
  select(tissue, discordant)
chisq_matrix_5 <- table(chisq_data_5$tissue, chisq_data_5$discordant)
chisq.test(chisq_matrix_5)
chisq.test(chisq_matrix_5)$expected
chisq.test(chisq_matrix_5)$residuals

length(which(duplicates_comparable$discordant == "discordant"))
length(which(duplicates_comparable$discordant == "concordant"))
1684/(1355+1684)

# 55% 
## 12-01: families attempt----
# strings of duplicate pairs
#  graph based approach
# create a list of gene pairs
edges <- unique(cbind(duplicate_pairs$Duplicate.1, duplicate_pairs$Duplicate.2))
# graph
gene_graph <- graph_from_edgelist(edges, directed = FALSE)
#  connected components
components <- components(gene_graph)
#  families
num_families <- components$no
cat("Number of gene families:", num_families, "\n")
# Get genes in each family
family_list <- split(names(components$membership), components$membership)
# Convert to a dataframe
max_family_size <- max(sapply(family_list, length))
family_df <- do.call(rbind, lapply(family_list, function(x) {
  c(x, rep(NA, max_family_size - length(x))) # Pad with NA for shorter families
}))
colnames(family_df) <- paste0("Gene", 1:ncol(family_df))
family_df <- as.data.frame(family_df)

mean(rowSums(!is.na(family_df)))
# on average, 6 genes in a family
max(rowSums(!is.na(family_df)))
# the largest family has 552 genes :')

## venn diagram----
duplicated_genes <- table$Gene[table$duplicate_status == "paralog"]
sex_biased_genes <- table$Gene[table$diff_exp_sexes == "sex-biased"]
all_genes <- nrow(table)


duplicated_genes <- table_edited$Gene[table_edited$duplicate_status == "paralog"]
sex_biased_genes <- table_edited$Gene[table_edited$diff_exp_sexes == "sex-biased"]
all_genes <- nrow(table_edited)
windows()

png("msc/firststep/plots/venn_diagram_gonad.png", width = 800, height = 600)
grid.newpage()
# Draw the Venn diagram without default labels
venn.plot <- draw.pairwise.venn(
  area1 = length(duplicated_genes),
  area2 = length(sex_biased_genes),
  cross.area = length(intersect(duplicated_genes, sex_biased_genes)),
  category = c("Paralogs (4485)", "Sex-biased genes (4351)"),
  fill = c("White", "Red"),
  alpha = 0.7,
  lwd = 0.5,
  cex = 2,
  cat.cex = 2,
  ind = TRUE,
  resolution = 2000,
  cat.pos = c(-20, 20),
  cat.dist = 0.05
)
dev.off()


# Significance of overlap
total <- all_genes
sex_biased <- 3493
paralogs <- 2604
overlap <- 1449

# p-value
pval <- phyper(overlap - 1, sex_biased, total - sex_biased, paralogs, lower.tail = FALSE)
print(pval)


## all plots -----
volcano_plots <- grid.arrange(volcanoplot, volcanoplot_antenna, volcanoplot_gut, volcanoplot_brain, nrow = 2, ncol = 2)
ggsave( "msc/firststep/plots/volcano_panel_samplesize.png", plot = volcano_plots, dpi = 600, width = 10, height = 7)

expression_plots <- grid.arrange(expression_patterns_plot, expression_patterns_plot_antenna, expression_patterns_plot_brain,
                                 expression_patterns_plot_gut, nrow = 2, ncol = 2)
ggsave( "msc/firststep/plots/expression_plots.png", plot = expression_plots, dpi = 600, width = 10, height = 7)

all_pairs_tissues$tissue <- factor(all_pairs_tissues$tissue, levels = c("gonad", "antenna", "gut", "brain"))
(disc_vs_conc_tissues  <- ggplot(all_pairs_tissues, aes(x = tissue, fill = discordant))+
    geom_bar(position= "fill", stat = "count" )+
    scale_fill_manual(values = c( "#7EC0EE", "#FFA500"))+
    labs(x = "\nTissue", y = "Percentage of duplicate pairs\n")+
    guides(fill = guide_legend(title = "Expression pattern")) +
    scale_y_continuous(labels = percent)+
    theme_classic()+
    theme(axis.text = element_text(size = 12),
          legend.title = element_text(size = 12),  # Increase legend title size
          legend.text = element_text(size = 12),
          axis.title = element_text(size = 12)))
ggsave("msc/firststep/plots/discordant_concordant_by_tissue.png", dpi = 600, width = 8, height = 5)


pies <- grid.arrange(pie_gonad, pie_antenna, pie_gut, pie_brain, nrow = 2, ncol = 2)
ggsave( "msc/firststep/plots/pies.png", plot = pies, dpi = 600, width = 10, height = 7)


dup_types_all_tissues <- grid.arrange(discordant_plot, discordant_plot_antenna, discordant_plot_gut, nrow = 1, ncol = 3)
ggsave( "msc/firststep/plots/dup_types_disc_v3.png", plot = dup_types_all_tissues, dpi = 600, width = 15, height = 5)


summary_fm <- table %>% 
  group_by(sex_bias) %>% 
  summarise(count = length (sex_bias))
summary_fm

table$Gene <- rownames(table)
genes_only_gonad <- table %>% 
  filter(! (Gene %in% table_antenna$Gene) ) %>% 
  filter(! (Gene %in% table_gut$Gene) ) %>% 
  filter(! (Gene %in% table_brain$Gene) )

summary_gonadspec <- genes_only_gonad %>% 
  group_by(sex_bias) %>% 
  summarise(count = length (sex_bias))
summary_gonadspec



# randomisation -----

# Function to generate a single dataset of l pairs without replacement
generate_random_pairs <- function(table, l) {
  sex_bias_pool <- table$sex_bias
  
  # Check if there are enough unique values to form l pairs
  if (length(sex_bias_pool) < 2 * l) {
    stop("not enough unique values in sex_bias to form the requested number of pairs without replacement.")
  }
  
  # Shuffle the sex_bias_pool to randomize the order
  sex_bias_pool <- sample(sex_bias_pool)
  
  # Select values in pairs without replacement
  random_pairs <- data.frame(
    Pair1 = sex_bias_pool[seq(1, 2 * l, by = 2)],
    Pair2 = sex_bias_pool[seq(2, 2 * l, by = 2)],
    stringsAsFactors = FALSE
  )
  
  return(random_pairs)
}

# create 1000 datasets
generate_multiple_datasets <- function(table, l, num_datasets) {
  datasets <- vector("list", num_datasets) # Initialize an empty list
  
  for (i in 1:num_datasets) {
    datasets[[i]] <- generate_random_pairs(table, l)
  }
  
  return(datasets)
}


set.seed(123) 
l <- nrow(duplicates_comparable_antenna) 
num_datasets <- 1000 

all_datasets <- generate_multiple_datasets(table_antenna, l, num_datasets)

View(all_datasets[[1]])

print(random_pairs_dataset)

calculate_discordance <- function(dataset) {
  discordant_pairs <- sum(dataset$Pair1 != dataset$Pair2)
  total_pairs <- nrow(dataset)
  return((discordant_pairs / total_pairs) * 100) 
}

observed_value <- 26.87

discordance_percentages <- sapply(all_datasets, calculate_discordance)

# Calculate the p-value: proportion of random datasets with discordance >= observed_value
p_value <- mean(discordance_percentages >= observed_value)
p_value

range(discordance_percentages)
View(duplicates_comparable_gut)
summary_gonad <- duplicates_comparable %>% 
  group_by(discordant) %>% 
  summarise(length(discordant))
summary_gonad
691/(691+1880)
# antenna 26.87%
1684/(1355+1684)
55.4%

