# Characterising the link between gene duplication and sex-biased gene expression in Timema stick insects

## Project Summary
Intralocus sexual conflict arises when selection favours different optimal trait values for the same genetic locus in males and females, preventing both sexes from reaching their fitness optima. Gene duplication might resolve such conflicts, with one or both copies of a duplicate pair potentially evolving sex-biased expression and thereby contributing to the evolution of sexual dimorphism. In this project, I explored the link between gene duplication and sex-biased gene expression in Timema poppense, a non-model organism, by identifying paralog pairs in the genome and quantifying gene expression using RNAseq data. I found that duplicated genes were consistently more likely to be sex-biased than singletons in the gonad and three somatic tissues. Discordant patterns (where the two copies have different expression biases) were prevalent in all tissues, with unbiased-female biased pairs more common than expected based on prior studies. This might suggest stronger sexual, or sex-specific selection acting on females in this species, but further studies are recommended to confirm this. Contrary to expectations, pairs where one copy was relocated were not more likely to have discordant expression, but difficulties in categorising duplication mechanisms might limit this analysis. Overall, my findings support a role for gene duplication in resolving intralocus sexual conflict.

## Main Results
### Differential gene expression in four tissues
<img width="945" height="661" alt="image" src="https://github.com/user-attachments/assets/54c16bc7-440b-4157-a06a-b11b3941a5f0" />
Differential gene expression between the sexes across four tissues. Points represent expressed genes, with the x-axis representing the log2 fold change and the y-axis representing the -log10 adjusted p value. Genes significantly upregulated in males and females are shown in red and blue, respectively, based on a combined threshold of p <= 0.05 and log2fold change >1 for male-bias and <1 for female-bias. Genes with unbiased expression are shown in gray. 

### Overlap between duplicate-status and differential gene expression in four tissues
<img width="870" height="601" alt="image" src="https://github.com/user-attachments/assets/cada9afd-b4ca-491a-8c8e-7c8e670b0797" />
Overlap between paralogous genes and sex-biased genes in the four tissues. The total number of expressed genes in each tissue is indicated in parentheses. In all tissues, the observed overlap was significantly greater than expected under neutrality (hypergeometric test, p < 0.05 in all tissues). 

### Differential gene expression of duplicate pairs: concordant vs discordant expression
<img width="742" height="464" alt="image" src="https://github.com/user-attachments/assets/5c1e7f53-ae1c-4d69-89ce-082c4a5feae4" />
Concondant and discordant patterns of duplicate pairs across four tissues. The Y axis shows the percentage of pairs showing specific expression patterns relative to all expressed pairs in each tissue (i.e. duplicate pairs where both members are expressed in the given tissue). The total number of pairs (corresponding to 100%) is n = 3039 in the gonad, n = 1978 in the antenna, n = 2571 in the gut, and n = 2653 in the brain). 

### Specific expression patterns
<img width="775" height="593" alt="image" src="https://github.com/user-attachments/assets/37354260-0716-416d-95c7-1f2f480b9002" />
Specific expression patterns of duplicate pairs in the four tissues. „U” represents „unbiased”, „F” „female-biased” and „M” „male-biased” expression, with the combinations of these letters referring to expression patterns of paralogous pairs. 

## Code






