#Load libraries
library(biomformat)
library(readr)
library(phyloseq)
library(dplyr)
library(ggplot2)
library(forcats)
library(vegan)
library(ANCOMBC)

# Read the BIOM"phyloseq"# Read the BIOM file
biom_data <- read_biom("data/bracken.biom")
ps <- import_biom(biom_data)

# Read the SRA Metadata
metadata_df <- as.data.frame(read_tsv("data/metadata.tsv"))
rownames(metadata_df) <- metadata_df$Sample_ID
metadata_df$Sample_ID <- NULL

# Combine BIOM file with SRA Metadata
sample_names(ps)
rownames(metadata_df)

sample_data(ps) <- metadata_df
#sample_data(ps)[ , 2] <- sample_data(ps)[ ,1]
# Change Diet to factor
sample_data(ps)$Diet <- as.factor(sample_data(ps)$Diet)
class(sample_data(ps)$Diet)

# Rename taxonomy ranks
colnames(tax_table(ps)) <- c(
  "Kingdom",
  "Phylum",
  "Class",
  "Order",
  "Family",
  "Genus",
  "Species"
)

# Remove prefixes (eg. "p_"phylum)
ps@tax_table <- substring(ps@tax_table, 4)

# Rarefaction
otu_table <- as.data.frame(t(otu_table(ps)))
rare_curve <- rarecurve(otu_table, step = 1000, label = TRUE)
# all curves plateau

# Relative Abundance
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
ps_phy <- tax_glom(ps_rel, taxrank = "Phylum")
df <- psmelt(ps_phy)

# Retreieve top 10 Phyla by abundance
top_phyla <- df %>%
  group_by(Phylum) %>%
  summarise(total = sum(Abundance)) %>%
  arrange(desc(total)) %>%
  slice(1:10) %>%
  pull(Phylum)

# Group remaining Phyla into "Other"
df$Phylum <- ifelse(df$Phylum %in% top_phyla, df$Phylum, "Other")

# Plot Relative Abundance by Phylum
ggplot(df, aes(x = Sample, y = Abundance, fill = fct_reorder(Phylum, Abundance, .fun = sum, .desc = TRUE))) +
  geom_bar(stat = "identity", position = "stack") +
  labs(y = "Relative Abundance", x = "Sample", fill = "Phylum") +
  facet_wrap(~Diet, scales="free_x")

# Alpha Diversity
plot_richness(ps)

# Calculate Shannon and Simpson Diversity
alpha_div <- estimate_richness(ps, measures = c("Shannon", "Simpson"))

# Add Samples and Diet to alpha diversity df
alpha_div$Sample <- rownames(alpha_div)
alpha_div$Diet <- sample_data(ps)$Diet

# Calculate average Shannon diversity 
mean_Shannon <- alpha_div %>%
  group_by(Diet) %>%
  summarise(mean(Shannon))

# Calculate average Simpson diversity 
mean_Simpson <- alpha_div %>%
  group_by(Diet) %>%
  summarise(mean(Simpson))

# Test if mean alpha diversity differs between dietary groups
wilcox.test(Shannon ~ Diet, data = alpha_div)
wilcox.test(Simpson ~ Diet, data = alpha_div)

# Plot Shannon Diversity
ggplot(alpha_div, aes(x = Sample, y = Shannon)) +
  geom_point() +
  labs(y = "Shannon Diversity", x = "Sample") +
  facet_wrap(~Diet, scales="free_x")

# Plot Simpson Diversity
ggplot(alpha_div, aes(x = Sample, y = Simpson)) +
  geom_point() +
  labs(y = "Shannon Diversity", x = "Sample") +
  facet_wrap(~Diet, scales="free_x")

# Plot Alpha Diversity by 
plot_richness(ps, x = "Diet", measures = c("Shannon", "Simpson"))

# Beta Diversity
# PCoA with Bray-Curtis
ord.pcoa.bray <- ordinate(ps, method="PCoA", distance="bray")
plot_ordination(ps, ord.pcoa.bray, color = "Diet", title="Bray-Curtis PCoA") + geom_point(size = 3)

# PERMANOVA
metadata <- as(sample_data(ps), "data.frame")
adonis2(phyloseq::distance(ps, method = "bray") ~ Diet,
        data = metadata,
        permutations = 999)

# Differential Abundance
ancombc.out <- ancombc2(data = ps, tax_level = "Genus",
                        fix_formula = "Diet", rand_formula = NULL,
                        p_adj_method = "holm", pseudo_sens = TRUE,
                        prv_cut = 0, lib_cut = 1000, s0_perc = 0.05,
                        group = "Diet", struc_zero = TRUE, neg_lb = TRUE)

# Retrieve structural zeroes
ancombc.out$zero_ind

# Subset differences between diets
subset(ancombc.out$zero_ind, 
       `structural_zero (Diet = Omnivore)` != `structural_zero (Diet = Vegan)`)

# Review results
ancombc.out$res
ancombc.out$res$taxon

# Subset to only significant differences between diets (q < 0.05)
ancombc.sig <- subset(ancombc.out$res, q_DietVegan < 0.05)
ancombc.sig # zero results

 # Extract ANCOMBC results
res_df <- ancombc.out$res                               

# Create levels of higher taxa                          
levels <- c("Phylum", "Class", "Order", "Family")

# Clean taxa names                                  
res_df$taxon_clean <- sapply(res_df$taxon, function(x) {
  if (!grepl("_", x)) {
    return(paste0("Genus: ", x))
  }
  
  parts <- unlist(strsplit(x, "_"))
  parts <- parts[parts != ""]
  
  level <- levels[length(parts)]
  lowest <- tail(parts, 1)
  
  return(paste0(level, ": ", lowest))
})

# Plot all log-fold changes
ggplot(res_df, aes(x = lfc_DietVegan, y = reorder(taxon_clean, lfc_DietVegan))) +
  geom_point(aes(color = q_DietVegan < 0.05), size = 3) +
  geom_errorbar(aes(xmin = lfc_DietVegan - se_DietVegan, 
                    xmax = lfc_DietVegan + se_DietVegan)) +
  geom_vline(xintercept = 0, color = "red") +
  labs(x = "Log Fold Change (DietVegan)", 
       y = "Genus")

# Add column for absolute log-fold change 
res_df$abs_lfc <- abs(res_df$lfc_DietVegan)

# Extract top 20 taxa by decreasing absolute log-fold change
top_taxa <- res_df[order(res_df$abs_lfc, decreasing = TRUE), ][1:20, ]

# Plot top 20 taxa 
ggplot(top_taxa, aes(x = lfc_DietVegan, y = reorder(taxon_clean, lfc_DietVegan))) +
  geom_point(aes(color = q_DietVegan < 0.05), size = 3) +
  geom_errorbar(aes(xmin = lfc_DietVegan - se_DietVegan, 
                    xmax = lfc_DietVegan + se_DietVegan)) +
  geom_vline(xintercept = 0, color = "red") +
  labs(x = "Log Fold Change (DietVegan)", 
       y = "Genus")
