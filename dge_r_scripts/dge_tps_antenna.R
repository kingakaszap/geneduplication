## import data& merge ----
# 1 - htseqcount data
F_A_1029<- read.table("msc/firststep/data/htseq_tps_antenna/Tps_F_A_Ad_18-1029.txt", stringsAsFactors=F,header=F,sep="\t")
F_A_1029 <- F_A_1029 %>% rename(Gene = V1, Count = V2) %>% 
  mutate(sample = "F_A_1029")

F_A_1030 <- read.table("msc/firststep/data/htseq_tps_antenna/Tps_F_A_Ad_18-1030.txt", stringsAsFactors=F,header=F,sep="\t")
F_A_1030 <- F_A_1030 %>% rename(Gene = V1, Count = V2) %>% 
  mutate(sample = "F_A_1030" )

F_A_1031 <- read.table("msc/firststep/data/htseq_tps_antenna/Tps_F_A_Ad_18-1031.txt", stringsAsFactors=F,header=F,sep="\t")
F_A_1031 <- F_A_1031 %>% rename(Gene = V1, Count = V2) %>% 
  mutate(sample = "F_A_1031" )

F_A_1041 <- read.table("msc/firststep/data/htseq_tps_antenna/Tps_F_A_Ad_18-1041.txt", stringsAsFactors=F,header=F,sep="\t")
F_A_1041 <- F_A_1041 %>% rename(Gene = V1, Count = V2) %>% 
  mutate(sample = "F_A_1041")

M_A_1033 <- read.table("msc/firststep/data/htseq_tps_antenna/Tps_M_A_Ad_18-1033.txt", stringsAsFactors=F,header=F,sep="\t")
M_A_1033 <- M_A_1033 %>% rename(Gene = V1, Count = V2) %>% 
  mutate(sample ="M_A_1033")

M_A_1035 <- read.table("msc/firststep/data/htseq_tps_antenna/Tps_M_A_Ad_18-1035.txt", stringsAsFactors=F,header=F,sep="\t")
M_A_1035 <- M_A_1035 %>% rename(Gene = V1, Count = V2) %>% 
  mutate(sample = "M_A_1035")

M_A_1036 <- read.table("msc/firststep/data/htseq_tps_antenna/Tps_M_A_Ad_18-1036.txt", stringsAsFactors=F,header=F,sep="\t")
M_A_1036 <- M_A_1036 %>% rename(Gene = V1, Count = V2) %>% 
  mutate(sample = "M_A_1036")

M_A_1042 <- read.table("msc/firststep/data/htseq_tps_antenna/Tps_M_A_Ad_18-1042.txt", stringsAsFactors=F,header=F,sep="\t")
M_A_1042 <- M_A_1042 %>% rename(Gene = V1, Count = V2) %>% 
  mutate(sample = "M_A_1042")

# combine df-s
list_of_dfs <- list(F_A_1029, F_A_1030,F_A_1031,F_A_1041,M_A_1033,M_A_1035,M_A_1036, M_A_1042)
all_counts <- list_of_dfs %>%
  bind_rows() %>%  
  select(Gene, Count, sample) %>%  #only relevant columns
  pivot_wider(names_from = sample, values_from = Count)  


# 2 - dupgenefinder
duplicate_pairs <- read.table("msc/firststep/data/duplicate_pairs.txt", stringsAsFactors=F,header=T,sep="\t")

## filter by gene length ----
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

# trying to add back gene names idk where they went..
# convert the RPKM result into a data frame
rpkm_norm_df <- as.data.frame(rpkm_norm)
# add the gene names as a new column 
rpkm_norm_df$Gene <- expr$genes$Gene  # from expr$genes$Gene
# I hope this adds them back in the right order..

write.table(rpkm_norm_df, file="msc/firststep/rpkm_antenna.txt",quote=F, sep="\t")



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
counts_filtered_antenna <- as.data.frame(all_counts[all_counts$Gene %in% rpkm_filtered$Gene, ])
write.table(counts_filtered_antenna, file="msc/firststep/data/filtered_counts_antenna.txt",quote=F, sep="\t")

## table on whether genes are DE or not ----
# set gene column as rownames again
rownames(counts_filtered_antenna) <- counts_filtered_antenna$Gene
counts_filtered_antenna <- counts_filtered_antenna %>% select(-Gene)
# information on which group each sample belongs to
conditions <- factor(c("F","F","F","F","M","M","M","M"))
# does this understand that the first 4 samples are female the next 4 male?

# Create edgeR object which will contain count data and normalization factors
expr <- DGEList(counts=counts_filtered_antenna,group=conditions)
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
table_antenna <- et$table
# table$Padj <- p_FDR original code but this wasnt defined before so i assume its p_adj?
table_antenna$Padj <- p_adj
table_antenna$Gene  <- rownames(table_antenna) # gene column back again
write.table(table_antenna, file="msc/firststep/data/differential_expression_antenna.txt",quote=F, sep=",")

## DE stats on all genes----

# % of all genes that are DE
diff_expr_genes_antenna <- table_antenna %>% filter (Padj <= 0.05)
nrow(diff_expr_genes_antenna)/nrow(table_antenna) # 14 %

# add sex-biased/not
table_antenna<- table_antenna %>% 
  mutate(diff_exp_sexes = ifelse (Padj <= 0.05, "sex-biased", "unbiased")) %>% 
  # add paralog status
  mutate(duplicate_status = ifelse(
    ((Gene %in%duplicate_pairs$Duplicate.1)|(Gene %in% duplicate_pairs$Duplicate.2)),"paralog","singleton"))

# Bias type 
table_antenna$sex_bias <- "unbiased"
# if log2Foldchange > 1 and p < 0.05, male biased
table_antenna$sex_bias[table_antenna$logFC > 1 & table_antenna$Padj < 0.05] <- "M"
# if log2Foldchange < -1 and p < 0.05, female biased
table_antenna$sex_bias[table_antenna$logFC < 1 & table_antenna$Padj < 0.05] <- "F"

sum_for_volcano <- table_antenna %>% 
  group_by(sex_bias) %>% 
  summarise(count = n()) %>% 
  mutate(percentage = round((count / nrow(table_antenna)) * 100, 1))
sum_for_volcano
# Volcano plot
(volcanoplot_antenna <- ggplot(table_antenna, aes(x=logFC, y = -log10(Padj), col = sex_bias))+
    geom_point(size = 0.7)+
    theme_classic()+
    ylim(c(-8, 60))+
    xlim(c(-10,10))+
    theme(axis.title.y = element_blank(),
          legend.position = "none",
          axis.title.x = element_blank(),
          plot.title = element_text(size = 13))+
    ggtitle("b) Antenna (n = 7237)\n")+
    labs(x = "\nLog2FoldChange", y = "-log10 P(adj)\n")+
    scale_colour_manual(values = c("firebrick", "dodgerblue", "gray60"))+
    annotate("label", x = 5, y = 25, 
             label = "3.7%", 
             color = "black", fill = "white", label.size = 0.5) +
    annotate("label", x = -5, y = 25, 
             label = "10.5%", 
             color = "black", fill = "white", label.size = 0.5) +
    annotate("label", x = 0, y = -5, 
             label = "85.8%", 
             color = "black", fill = "white", label.size = 0.5)+
    geom_hline(yintercept=-log10(0.05), col="red"))
ggsave("msc/firststep/plots/tps_volcano_antenna.png", width = 8, height =6, dpi = 600)



## duplicate expression patterns ----
# filter duplicates that are expressed
duplicates_comparable_antenna <- duplicate_pairs %>% filter (Duplicate.1 %in% table_antenna$Gene & Duplicate.2 %in% table_antenna$Gene)
# 2571 

# what proportion of duplicates are expressed in antenna
nrow(duplicates_comparable_antenna)/nrow(duplicate_pairs) # 25%

# what proportion of dups are DE
diff_expr_duplicates_antenna <- diff_expr_genes_antenna %>%
  filter(Gene %in% duplicates_comparable_antenna$Duplicate.1 | Gene %in% duplicates_comparable_antenna$Duplicate.2)
nrow(diff_expr_duplicates_antenna)/nrow(duplicates_comparable_antenna)
# 20 % are sex-biased

# what proportion of expressed genes are dups
nrow(duplicates_comparable_antenna)/nrow(table_antenna) #35.5%

# what proportion of DE genes are duplicates
nrow(diff_expr_duplicates_antenna)/nrow(diff_expr_genes_antenna) # 50 %

# what proportion of unbiased genes are duplicates
(nrow(duplicates_comparable_antenna)- nrow(diff_expr_duplicates_antenna))/(nrow(table_antenna)- nrow(diff_expr_genes_antenna))
# 33%


# what proportion of singletons are DE
singletons <- table_antenna %>% filter (duplicate_status == "singleton")
(singletons_summary  <- singletons %>% 
    group_by(diff_exp_sexes) %>% 
    summarise(count = length(diff_exp_sexes)))
370/(370+3018) # 10.9%

# add new column to dups dataframe,: is duplicate 1 sex-biased? is duplicate 2 sex-biased?
duplicates_comparable_antenna<- duplicates_comparable_antenna %>% mutate( 
  Dup1_de = ifelse(Duplicate.1 %in% diff_expr_genes_antenna$Gene, "yes", "no"),
  Dup2_de = ifelse (Duplicate.2 %in%diff_expr_genes_antenna$Gene, "yes", "no"))

# and specific pattern
duplicates_comparable_antenna <- duplicates_comparable_antenna %>% 
  mutate(pattern = case_when(
    (Dup1_de == "yes"& Dup2_de == "yes") ~ 1,
    (Dup1_de == "no"& Dup2_de == "no") ~ 0,
    ((Dup1_de == "yes"& Dup2_de == "no") | (Dup1_de == "no"& Dup2_de == "yes")) ~2))

# Add logfc values
logfc_to_add <- diff_expr_genes_antenna %>% select(logFC, Gene) %>% rename(Duplicate.1 = Gene)
duplicate_pairs_logfc <- left_join(duplicates_comparable_antenna, logfc_to_add, by = "Duplicate.1")
duplicate_pairs_logfc <- rename(duplicate_pairs_logfc, logfc_dup1 = logFC)
logfc_to_add <- rename(logfc_to_add, Duplicate.2 = Duplicate.1)
duplicate_pairs_logfc_both <- left_join(duplicate_pairs_logfc, logfc_to_add, by = "Duplicate.2")
duplicate_pairs_logfc_both <- rename(duplicate_pairs_logfc_both, logfc_dup2 = logFC)
duplicates_comparable_antenna <- duplicate_pairs_logfc_both

# adding expression pattern
duplicates_comparable_antenna <- duplicates_comparable_antenna %>% 
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
duplicates_comparable_antenna<- duplicates_comparable_antenna %>% 
  mutate(discordant = ifelse(exp_bias_dup1 == exp_bias_dup2, "concordant", "discordant"))
(discordant_summary <- duplicates_comparable_antenna %>% 
    group_by(discordant) %>% 
    summarise(count = length(discordant)))
691/(691+1880)
# 26.87 % show discordant expression

# pattern in transposed genes
transposed_exp <- duplicates_comparable_antenna %>% filter(duplication_type == "transposed") %>% 
  filter(pattern == 2)
summary_transposed <- transposed_exp %>% 
  group_by(Dup1_de) %>% 
  summarise(count = length(Dup1_de))
summary_transposed
# when only 1 transposed copy is DE,
# transposed is DE = 60,
# parental is DE = 49

## plots about duplicates -----
### 1 - proportions of duplication types
duplicate_proportions_expressed_antenna <- duplicates_comparable_antenna %>% 
  group_by(duplication_type) %>%
  summarise(percentage = 100* n() / nrow(duplicate_pairs))
# Plot proportions
(expressed_duplicates_barplot_antenna <- ggplot(duplicate_proportions_expressed_antenna, aes(x = duplication_type, y = percentage, fill = duplication_type)) +
    geom_bar(stat = "identity") +  
    labs(x = "\nDuplication type", y = "% of all pairs\n") +
    scale_fill_brewer(palette = "Paired") +
    theme_classic())
ggsave("msc/firststep/plots/gonad_duplication_types_proportions_antenna.png", dpi = 500)

# 2 - basic patterns in duplicates that are expressed
duplicates_comparable_antenna$pattern<- as.factor( duplicates_comparable_antenna$pattern)
duplicates_comparable_antenna$duplication_type<- as.factor(duplicates_comparable_antenna$duplication_type)
(basic_plot <- ggplot (duplicates_comparable_antenna, aes(x = pattern, fill = duplication_type))+
    geom_bar(stat = "count")+
    labs(x = "\nPattern", y = "Count\n")+
    guides(fill = guide_legend(title = "Dupication type"))+
    scale_x_discrete(labels = c("Both no", "Both yes", "One yes, one not")) +
    theme_classic())

# 3 - Specific patterns
# get proportions for each bias pattern and duplication_type
duplicates_percent_antenna <- duplicates_comparable_antenna %>%
  group_by(bias_pattern, duplication_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(percent = 100 * count / sum(count))  

# Plot
(expression_patterns_plot_antenna <- ggplot(duplicates_percent_antenna, aes(x = bias_pattern, y = percent, fill = duplication_type)) +
    geom_bar(stat = "identity") +  # Use calculated percentages
    labs(x = "\nPattern", y = "Percentage (%)\n") +
    ggtitle("b) Antenna\n")+
    scale_fill_brewer(palette = "Paired") +
    guides(fill = guide_legend(title = "Duplication type")) +
    theme_classic()+
    theme(legend.position = "none",
          plot.title = element_text(size = 13),
          axis.text = element_text(size = 12),
          axis.title.y = element_blank(),
          axis.title.x = element_blank()))
ggsave("msc/firststep/plots/tps_antenna_expression_patterns_dups.png", dpi=600)

(duplication_types_plot_antenna <- ggplot(duplicates_percent_antenna, aes(x = duplication_type, y = percent, fill = bias_pattern)) +
    geom_bar(stat = "identity") +  # Use calculated percentages
    labs(x = "\nPattern", y = "Percentage (%)") +
    scale_fill_brewer(palette = "Paired") +
    guides(fill = guide_legend(title = "Duplication type")) +
    theme_classic())
ggsave("msc/firststep/plots/tps_gonad_duplication_types_expression_antenna.png", dpi=600)

for_discordant_plot_antenna <- duplicates_comparable_antenna %>% 
  filter( ! duplication_type == "other")

(discordant_plot_antenna <- ggplot(for_discordant_plot_antenna, aes(x = duplication_type, fill = discordant))+
    geom_bar(position = "fill", stat = "count")+
    labs (x = "\nDuplication type", y = "Percentage\n")+
    scale_y_continuous(labels = scales::percent) +  
    guides(fill = guide_legend(title = "Expression pattern")) +
    scale_fill_brewer(palette = "Paired")+
    theme_classic()+
    theme(legend.position = "bottom",
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14)))
ggsave("msc/firststep/plots/tps_antenna_discordant.png", dpi=600, width = 6, height = 6)


# pie
duplicates_for_pie_antenna <- duplicates_percent_antenna %>%
  group_by(bias_pattern) %>%
  summarise(percent = sum(percent))

# Reorder the bias_pattern factor
duplicates_for_pie_antenna$bias_pattern <- factor(
  duplicates_for_pie_antenna$bias_pattern, 
  levels = c("UU", "FF", "MM", "UM", "UF", "FM"))
duplicates_for_pie_antenna <- duplicates_for_pie_antenna %>%
  mutate(
    label = ifelse(bias_pattern %in% c("UU", "UF"),  # For "UU" and "UF"
                   paste0(bias_pattern, "\n", round(percent, 1), "%"),  # Show both label and percentage
                   as.character(bias_pattern)))  # For other categories, show only the label
View(duplicates_for_pie_antenna)
(pie_antenna <- ggplot(duplicates_for_pie_antenna, aes(x = "", y = percent, fill = bias_pattern))+
    geom_bar(stat = "identity", width = 1)+
    coord_polar(theta = "y")+
    labs (x = NULL, y = NULL, title = "b) Antenna")+
    # geom_text(aes(label = label),  # Use the new label column
 #  position = position_stack(vjust = 0.5),  # Place text in the middle of each slice
 #  size = 4) +   
    geom_text(aes(label = ifelse(bias_pattern %in% (c("UU", "UF", "UM")), 
                                 paste0(bias_pattern, "\n", round(percent, 1), "%"), NA)),
              position = position_stack(vjust = 0.5),
              size = 4) +
    guides(fill = guide_legend(title = "Expression pattern")) +
    theme_classic()+
    scale_fill_brewer(palette = "Paired") +
    theme(axis.line = element_blank(),
          legend.position = "none",
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text (size = 13)))
ggsave("msc/firststep/plots/tps_antenna_expression_patterns_dups.png", dpi=600)

# dup types
for_discordant_plot_antenna <- duplicates_comparable_antenna %>% 
  filter(! duplication_type == "other")

for_discordant_plot_antenna_summary <- duplicates_comparable_antenna %>% 
  group_by (duplication_type) %>% 
  summarise(count = length(duplication_type))
for_discordant_plot_antenna_summary

(discordant_plot_antenna <- ggplot(duplicates_comparable_antenna, aes(x = duplication_type, fill = discordant))+
    geom_bar(position = "fill", stat = "count")+
    labs (x = "\nDuplication type", y = "Percentage\n")+
    scale_y_continuous(labels = scales::percent) +  
    guides(fill = guide_legend(title = "Expression pattern")) +
     scale_fill_brewer(palette = "Paired")+
    scale_x_discrete(labels = c("other\n (1779)",  "proximal\n (123)", "tandem\n (129)", "transposed\n (540)" ))+
    theme_classic()+
    ggtitle ("b) Antenna\n")+
    theme(legend.position = "none",
          axis.title = element_blank(),
          plot.title = element_text(size = 18),
          # axis.title = element_text(size = 14),
          axis.text = element_text(size = 14)))
ggsave("msc/firststep/plots/tps_antenna_discordant_duptypes.png", dpi = 600, width = 6, height =6)


## chisq tests ----
# Paralog/singleton, sex-biased/unbiased
# step 1 - add duplicate/not duplicate status to counts_filtered
counts_filtered_antenna$Gene <- rownames(counts_filtered_antenna)
counts_filtered_antenna <- counts_filtered_antenna %>% mutate(
  duplicate = ifelse((Gene %in% duplicate_pairs$Duplicate.1 | Gene %in% duplicate_pairs$Duplicate.2),
                     "duplicated", "not_duplicated"),
  de = ifelse ((Gene %in% diff_expr_genes_antenna$Gene), "sex-biased", "unbiased"))
chisq_data <- counts_filtered_antenna %>% select(duplicate,de)
chisq_matrix <- table(chisq_data$duplicate, chisq_data$de)
chisq.test(chisq_matrix)
chisq_matrix
chisq.test(chisq_matrix)$expected
chisq.test(chisq_matrix)$residuals
# also very significant in antenna....

# Within sex-biased genes: dup/singleton, M/F bias
chisq_data_2<- table_antenna %>%
  # filter(sex_bias != "unbiased") %>% 
  select(duplicate_status, sex_bias)
chisq_matrix_2 <- table(chisq_data_2$duplicate_status, chisq_data_2$sex_bias)
chisq.test(chisq_matrix_2)
chisq.test(chisq_matrix_2)$residuals
chisq_matrix_2
# no association/ dif in whether sex-biased genes are female/male biased and whether they are duplicates/not.
# unlike wyman

chisq_data_3<- table_antenna %>%
  select(duplicate_status, sex_bias)
chisq_matrix_3 <- table(chisq_data_3$duplicate_status, chisq_data_3$sex_bias)
chisq.test(chisq_matrix_3)
chisq.test(chisq_matrix_3)$expected
chisq.test(chisq_matrix_3)$residuals
chisq_matrix_3
# so sexual dimorphism in antenna??

# discordant vs concordant
for_chisq_4 <- duplicates_comparable_antenna %>% 
  filter( (duplication_type != "other"))
chisq_data_4 <- for_chisq_4 %>% 
  select(duplication_type, discordant)

chisq_data_4 <- duplicates_comparable_antenna %>% 
  select(duplication_type, discordant)
chisq_matrix_4 <- table(chisq_data_4$duplication_type, chisq_data_4$discordant)
chisq.test(chisq_matrix_4)
chisq.test(chisq_matrix_4)$expected
chisq.test(chisq_matrix_4)$residuals
# but do I have a reason to exclude UU pairs?
# I dont think so

# venn diagram ----
library(VennDiagram)
library(grid)
View(table_antenna)
duplicated_genes <- table_antenna$Gene[table_antenna$duplicate_status == "paralog"]
sex_biased_genes <- table_antenna$Gene[table_antenna$diff_exp_sexes == "sex-biased"]
all_genes <- nrow(table_antenna)
length(duplicated_genes)
length(sex_biased_genes)

windows()

png("msc/firststep/plots/venn_diagram_antenna.png", width = 850, height = 600)
grid.newpage()
# Draw the Venn diagram without default labels
venn.plot <- draw.pairwise.venn(
  area1 = length(duplicated_genes),
  area2 = length(sex_biased_genes),
  cross.area = length(intersect(duplicated_genes, sex_biased_genes)),
  category = c("Paralogs (3849)", "Sex-biased genes\n(1031)"),
  fill = c("White", "Red"),
  alpha = c(0.5, 0.5),
  lwd = 0.5,
  cex = 2,
  cat.cex = 1.5,
  ind = TRUE,
  cat.pos = c(-20, 20),
  cat.dist = 0.05
)
dev.off()

total <- 7237
sex_biased <- 3493
paralogs <- 1031
overlap <- 661

# Calculate p-value
pval <- phyper(overlap - 1, sex_biased, total - sex_biased, paralogs, lower.tail = FALSE)
print(pval)

# sign
