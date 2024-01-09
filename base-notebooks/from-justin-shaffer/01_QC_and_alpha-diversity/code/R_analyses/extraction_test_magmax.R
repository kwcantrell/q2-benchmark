#########################################################################################
# 2019-2020
# Extraction kit testing - round 3 - MagMAX 2vs20-min vs. PowerSoil
# Justin Shaffer
# justinparkshaffer@gmail.com
#########################################################################################
#########################################################################################

# Set working directory
#########################################################################################
getwd()
setwd("~/Google-Drive-UCSD/R/2020_extraction_test_magmax/")


# Install and load libraries needed for analysis
#########################################################################################

install.packages("remotes")
remotes::install_github("jbisanz/qiime2R")
install.packages("tidyverse")
install.packages("vegan")
install.packages("cluster")
install.packages("dendextend")
install.packages("indicspecies")
install.packages("gdata")
install.packages("reshape2")

library(qiime2R)
library(tidyverse)
library(vegan)
library(cluster)
library(dendextend)
library(indicspecies)

library(plyr)

# For hierarchical clustering
# NOTE: Make sure plyr is OFF
library(vegan)
library(cluster)
library(dendextend)
library(gdata)
library(reshape2)


# Append read count and alpha-diversity vectors to sample metadata
#######################################################################################################

# Read in  metadata
md <- read_tsv("metadata_12201_round3_qiitaIDs_2020.09.23.txt")


# Read in read count data
read_counts_16S <- read_tsv("extraction_test_magmax_read_counts_16S.txt")
read_counts_shotgun <- read_tsv("extraction_test_magmax_read_counts_shotgun.txt")


# Merge read count data with metadata
md_read_counts_16S <- left_join(md, read_counts_16S, by = c("sample_name_noQiitaIDs" = "sample_name"))
md_read_counts <- left_join(md_read_counts_16S, read_counts_shotgun, by = c("sample_name_noQiitaIDs" = "sample_name"))


# Read in alpha-diversity metrics for 16S and shotgun data
alpha_16S_faithspd <- read_qza("extraction_test_round_3_biom_lod_rar5k_alpha_faithspd.qza")
alpha_16S_shannon <- read_qza("extraction_test_round_3_biom_lod_rar5k_alpha_shannon.qza")
alpha_16S_richness <- read_qza("extraction_test_round_3_biom_lod_rar5k_alpha_richness.qza")

alpha_shotgun_high_faithspd <- read_qza("gotu_profile_updated_sampleIDs_highBiomass_rar35K_alpha_faithspd.qza")
alpha_shotgun_high_shannon <- read_qza("gotu_profile_updated_sampleIDs_highBiomass_rar35K_alpha_shannon.qza")
alpha_shotgun_high_richness <- read_qza("gotu_profile_updated_sampleIDs_highBiomass_rar35K_alpha_richness.qza")

alpha_shotgun_low_faithspd <- read_qza("gotu_profile_updated_sampleIDs_lowBiomass_rar20K_alpha_faithspd.qza")
alpha_shotgun_low_shannon <- read_qza("gotu_profile_updated_sampleIDs_lowBiomass_rar20K_alpha_shannon.qza")
alpha_shotgun_low_richness <- read_qza("gotu_profile_updated_sampleIDs_lowBiomass_rar20K_alpha_richness.qza")


# Concatenate alpha-diversity metrics from low- and high-biomass samples for metagenomics data
alpha_shotgun_faithspd <- rbind((rownames_to_column(as.data.frame(alpha_shotgun_high_faithspd$data), "sample_name_shotgun")), (rownames_to_column(as.data.frame(alpha_shotgun_low_faithspd$data), "sample_name_shotgun")))

alpha_shotgun_shannon <- rbind((rownames_to_column(as.data.frame(alpha_shotgun_high_shannon$data), "sample_name_shotgun")), (rownames_to_column(as.data.frame(alpha_shotgun_low_shannon$data), "sample_name_shotgun")))

alpha_shotgun_richness <- rbind((rownames_to_column(as.data.frame(alpha_shotgun_high_richness$data), "sample_name_shotgun")), (rownames_to_column(as.data.frame(alpha_shotgun_low_richness$data), "sample_name_shotgun")))


# Combine alpha-diversity metrics with metadata
alpha_shotgun_faithspd %>%
  right_join(md_read_counts) -> md_read_counts

alpha_shotgun_shannon %>%
  right_join(md_read_counts) -> md_read_counts

alpha_shotgun_richness %>%
  right_join(md_read_counts) -> md_read_counts

colnames(md_read_counts)[2] <- "alpha_shotgun_rar35k_20k_richness"
colnames(md_read_counts)[3] <- "alpha_shotgun_rar35k_20k_shannon"
colnames(md_read_counts)[4] <- "alpha_shotgun_rar35k_20k_faithspd"

alpha_16S_faithspd$data %>%
  as.data.frame() %>%
  rownames_to_column("sample_name") %>%
  right_join(md_read_counts) -> md_read_counts

alpha_16S_shannon$data %>%
  as.data.frame() %>%
  rownames_to_column("sample_name") %>%
  right_join(md_read_counts) -> md_read_counts

alpha_16S_richness$data %>%
  as.data.frame() %>%
  rownames_to_column("sample_name") %>%
  right_join(md_read_counts) -> md_read_counts

colnames(md_read_counts)[2] <- "alpha_16S_rar5k_richness"
colnames(md_read_counts)[3] <- "alpha_16S_rar5k_shannon"
colnames(md_read_counts)[4] <- "alpha_16S_rar5k_faithspd"


# Export new metadata file
write_tsv(md_read_counts, path = "metadata_12201_round3_read_counts_alpha_diversity_2020.09.23.txt")


# Filter metadata for plotting
#######################################################################################################

# Read in metadata with read counts and alpha-diversity metrics
md <- read_tsv("metadata_12201_round3_read_counts_alpha_diversity_2020.09.23.txt")

# Re-order levels for sample biomass
md$biomass_sample <- factor(md$biomass_sample, levels = c("high", "low"))

# Re-order levels for sample type 2
md$sample_type_2 <- factor(md$sample_type_2, levels = c("doorknob", "floor", "sourkraut", "yogurt", "seawater", "freshwater", "bare soil", "rhizosphere soil", "mouse jejunum tissue", "mouse feces", "cat feces", "human feces", "human female urine", "human male urine", "human foot", "human armpit", "human lower arm", "human forehead", "human nares", "human throat", "human saliva", "human oral saline rinse", "human oral saline rinse in VTM", "single strain", "two strains", "PCR extraction control"))

# filter blanks
md_noBlanks <- subset(md, md$sample_type != "control blank")

# filter mock
md_noBlanks_noMock <- subset(md_noBlanks, md_noBlanks$sample_type != "cultured cells")

# filter failures
md_noBlanks_noZeros <- subset(md_noBlanks, md_noBlanks$dna_ng_ul != "0")
md_noBlanks_noMock_noZeros <- subset(md_noBlanks_noMock, md_noBlanks_noMock$dna_ng_ul != "0")


# Re-shape metadata for pairwise comparisons of continuous variables between protocols
#######################################################################################################
md_noBlanks_unique_sample_id <- melt(md_noBlanks, id.vars = c("unique_sample_id", "extraction_protocol_no_spaces", "biomass_sample_long", "sample_type", "sample_type_2", "sample_type_3"), measure.vars = c("dna_ng_ul", "rna_ng_ul", "read_count_16s", "read_count_shotgun", "alpha_16S_rar5k_richness", "alpha_16S_rar5k_shannon", "alpha_16S_rar5k_faithspd", "alpha_shotgun_rar35k_20k_richness", "alpha_shotgun_rar35k_20k_shannon", "alpha_shotgun_rar35k_20k_faithspd"), values_drop_na = FALSE)

md_noBlanks_unique_sample_id_noCovid <- subset(md_noBlanks_unique_sample_id, 
                                        md_noBlanks_unique_sample_id$sample_type_2 != "human forehead" &
                                        md_noBlanks_unique_sample_id$sample_type_2 != "human nares" &
                                        md_noBlanks_unique_sample_id$sample_type_2 != "human throat")

md_noBlanks_unique_sample_id_noCovid_wide <- pivot_wider(md_noBlanks_unique_sample_id_noCovid, c("unique_sample_id", "extraction_protocol_no_spaces", "biomass_sample_long", "sample_type", "sample_type_2", "sample_type_3"), names_from = c(extraction_protocol_no_spaces, variable), values_from = value, values_fn = mean)


# Plot DNA yield
#######################################################################################################

# Scatterplot - 1-to-1:
## x = dna_ng_ul_PS
## y = dna_ng_ul_MM_2
## facets = sample_type_2
ggplot(md_noBlanks_unique_sample_id_noCovid_wide, aes(x = PS_dna_ng_ul, y = MM_2_dna_ng_ul)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  geom_point(size = 0.75) +
  facet_wrap(~sample_type_2) +
  xlab("PowerSoil DNA yield (ng/µL) ") +
  ylab("MagMAX (2-min.) DNA yield (ng/µL) ") +
  #xlim(c(0,80)) +
  #ylim(c(0,80)) +
  theme_bw() +
  theme(legend.position = "none", legend.title = element_blank())


# Scatterplot - 1-to-1:
## x = log10(dna_ng_ul_PS + 1)
## y = log10(dna_ng_ul_MM_2 + 1)
## facets = sample_type_2
ggplot(md_noBlanks_unique_sample_id_noCovid_wide, aes(x = log10(PS_dna_ng_ul + 1), y = log10(MM_2_dna_ng_ul + 1))) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  geom_point(size = 0.75) +
  facet_wrap(~sample_type_2) +
  xlab("log(PowerSoil DNA yield [ng/µL]) ") +
  ylab("log(MagMAX [2-min.] DNA yield [ng/µL]) ") +
  #xlim(c(0,80)) +
  #ylim(c(0,80)) +
  theme_bw() +
  theme(legend.position = "none", legend.title = element_blank())


# Scatterplot:
## x = unique_sample_id
## y = dna_ng_ul
## color = extraction_protocol_long
## facets = sample_type_2
ggplot(md_noBlanks, aes(reorder(x = unique_sample_id, dna_ng_ul), y = dna_ng_ul, color = extraction_protocol)) +
  geom_point(size = 0.75) +
  facet_wrap(~sample_type_2, scales = "free") +
  xlab("sample") +
  ylab("DNA yield (ng/µL)") +
  scale_color_manual(values = c("orange", "green3", "blue"), aes(alpha = 0.5)) +
  theme_bw() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position = "top", legend.title = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=3)))


# Scatterplot:
## x = unique_sample_id
## y = dna_ng_ul
## color = extraction_protocol_long
## facets = biomass_sample_long
ggplot(md_noBlanks, aes(reorder(x = unique_sample_id, dna_ng_ul), y = dna_ng_ul, color = extraction_protocol)) +
  geom_point(size = 0.75) +
  facet_wrap(~biomass_sample_long, scales = "free", ncol = 1) +
  xlab("sample") +
  ylab("DNA yield (ng/µL)") +
  scale_color_manual(values = c("orange", "green3", "blue"), aes(alpha = 0.5)) +
  theme_bw() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position = "top", legend.title = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=3)))


# Boxplot:
## x = extraction_protocol_short
## y = dna_ng_ul
## facets = biomass_sample_long
ggplot(md_noBlanks_noMock_noZeros, aes(x = extraction_protocol_short, y = as.numeric(dna_ng_ul))) +
  geom_boxplot(outlier.alpha = 0.25) +
  facet_wrap(~biomass_sample_long, ncol = 2, scales = "free_y") +
  stat_summary(fun = mean, geom = "point", color = "black", fill = "red", shape = 21) +
  xlab("extraction protocol") +
  ylab("DNA yield (ng/µL)") +
  theme_bw() +
  theme(text = element_text(size = 14))


# Boxplot:
## x = extraction_protocol_short
## y = log10(dna_ng_ul + 1)
## facets = biomass_sample_long
ggplot(md_noBlanks_noMock_noZeros, aes(x = extraction_protocol_short, y = log10(as.numeric(dna_ng_ul) + 1))) +
  geom_boxplot(outlier.alpha = 0.25) +
  facet_wrap(~biomass_sample_long, ncol = 2, scales = "free_y") +
  xlab("extraction protocol") +
  ylab("log(DNA yield [ng/µL])") +
  theme_bw() +
  theme(text = element_text(size = 14))


# Boxplot:
## x = extraction_protocol_short
## y = dna_ng_ul
## facets = sample_type_2
md_noBlanks_noMock_noZeros_noArm <- subset(md_noBlanks_noMock_noZeros,
                                           md_noBlanks_noMock_noZeros$sample_type_2 != "human lower arm")
ggplot(md_noBlanks_noMock_noZeros_noArm, aes(x = extraction_protocol_short, y = as.numeric(dna_ng_ul))) +
  geom_boxplot(outlier.alpha = 0, aes(fill = extraction_protocol_short)) +
  geom_jitter(, width = 0.2, size = 0.5) +
  facet_wrap(~sample_type_2, scales = "free_y") +
  stat_summary(fun = mean, geom = "point", color = "black", fill = "red", shape = 21) +
  xlab("extraction protocol") +
  ylab("DNA yield (ng/µL)") +
  scale_fill_manual(values = c("orange", "green4", "blue")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 14),
        strip.text.x = element_text(size = 6),
        strip.background =element_rect(fill="white"))


# Scatterplot:
## x = sample_type_2
## y = dna_ng_ul
## facets = extraction_protocol~biomass_sample_long
ggplot(md_noBlanks_noMock_noZeros, aes(x = sample_type_2, y = dna_ng_ul)) +
  geom_point() +
  facet_wrap(extraction_protocol~biomass_sample_long, scales = "free", ncol = 2) +
  theme(axis.text.x = element_text(angle=35, hjust=1)) +
  xlab("sample type") +
  ylab("DNA yield (ng/µL)")


# Calculate means and SDs for DNA yield
#######################################################################################################
aggregate(as.numeric(dna_ng_ul) ~ biomass_sample_long:extraction_protocol_short, FUN = mean, data = md_noBlanks_noMock_noZeros)
aggregate(as.numeric(dna_ng_ul) ~ biomass_sample_long:extraction_protocol_short, FUN = sd, data = md_noBlanks_noMock_noZeros)


# Wilcoxon  test for differences in gDNA across extraction protocols
######################################################################################################

## Paired tests across all data except blanks and non-paired samples ##
wilcox.test(md_noBlanks_unique_sample_id_noCovid_wide$MM_2_dna_ng_ul, md_noBlanks_unique_sample_id_noCovid_wide$MM_20_dna_ng_ul, paired = TRUE)

wilcox.test(md_noBlanks_unique_sample_id_noCovid_wide$MM_2_dna_ng_ul, md_noBlanks_unique_sample_id_noCovid_wide$PS_dna_ng_ul, paired = TRUE)

wilcox.test(md_noBlanks_unique_sample_id_noCovid_wide$MM_20_dna_ng_ul, md_noBlanks_unique_sample_id_noCovid_wide$PS_dna_ng_ul, paired = TRUE)


# Subset data by sample biomass
md_noBlanks_noMock_noZeros_highBiomass <- subset(md_noBlanks_noMock_noZeros, md_noBlanks_noMock_noZeros$biomass_sample == "high")

md_noBlanks_noMock_noZeros_lowBiomass <- subset(md_noBlanks_noMock_noZeros, md_noBlanks_noMock_noZeros$biomass_sample == "low")


# Subset data further by extraction kit (high biomass)
md_noBlanks_noMock_noZeros_highBiomass_2v20 <- subset(md_noBlanks_noMock_noZeros_highBiomass, md_noBlanks_noMock_noZeros_highBiomass$extraction_protocol_short != "PS")

md_noBlanks_noMock_noZeros_highBiomass_2vPS <- subset(md_noBlanks_noMock_noZeros_highBiomass, md_noBlanks_noMock_noZeros_highBiomass$extraction_protocol_short != "MM (20-min.)")

md_noBlanks_noMock_noZeros_highBiomass_20vPS <- subset(md_noBlanks_noMock_noZeros_highBiomass, md_noBlanks_noMock_noZeros_highBiomass$extraction_protocol_short != "MM (2-min.)")


# Subset data further by extraction kit (low biomass)
md_noBlanks_noMock_noZeros_lowBiomass_2v20 <- subset(md_noBlanks_noMock_noZeros_lowBiomass, md_noBlanks_noMock_noZeros_lowBiomass$extraction_protocol_short != "PS")

md_noBlanks_noMock_noZeros_lowBiomass_2vPS <- subset(md_noBlanks_noMock_noZeros_lowBiomass, md_noBlanks_noMock_noZeros_lowBiomass$extraction_protocol_short != "MM (20-min.)")

md_noBlanks_noMock_noZeros_lowBiomass_20vPS <- subset(md_noBlanks_noMock_noZeros_lowBiomass, md_noBlanks_noMock_noZeros_lowBiomass$extraction_protocol_short != "MM (2-min.)")


# Perform tests (high biomass)
wilcox.test(as.numeric(md_noBlanks_noMock_noZeros_highBiomass_2v20$dna_ng_ul) ~ md_noBlanks_noMock_noZeros_highBiomass_2v20$extraction_protocol_short)

wilcox.test(as.numeric(md_noBlanks_noMock_noZeros_highBiomass_2vPS$dna_ng_ul) ~ md_noBlanks_noMock_noZeros_highBiomass_2vPS$extraction_protocol_short)

wilcox.test(as.numeric(md_noBlanks_noMock_noZeros_highBiomass_20vPS$dna_ng_ul) ~ md_noBlanks_noMock_noZeros_highBiomass_20vPS$extraction_protocol_short)


# Perform tests (low biomass)
wilcox.test(as.numeric(md_noBlanks_noMock_noZeros_lowBiomass_2v20$dna_ng_ul) ~ md_noBlanks_noMock_noZeros_lowBiomass_2v20$extraction_protocol_short)

wilcox.test(as.numeric(md_noBlanks_noMock_noZeros_lowBiomass_2vPS$dna_ng_ul) ~ md_noBlanks_noMock_noZeros_lowBiomass_2vPS$extraction_protocol_short)

wilcox.test(as.numeric(md_noBlanks_noMock_noZeros_lowBiomass_20vPS$dna_ng_ul) ~ md_noBlanks_noMock_noZeros_lowBiomass_20vPS$extraction_protocol_short)


# Plot RNA yield
#######################################################################################################

# boxplot - RNA yield by sample biomass (Figure 5)
ggplot(md_noBlanks_noMock, aes(x = extraction_protocol_short, y = rna_ng_ul)) +
  geom_boxplot(outlier.alpha = 0.25, aes(fill = extraction_protocol_short)) +
  facet_wrap(~biomass_sample_long, ncol = 2, scales = "free_y") +
  stat_summary(fun = mean, geom = "point", color = "black", fill = "red", shape = 21) +
  scale_fill_manual(values = c("orange", "green4", "blue")) +
  xlab("extraction protocol") +
  ylab("RNA yield (ng/µL)") +
  ylim(c(0,60)) +
  theme_bw() +
  theme(text = element_text(size = 14),
        legend.position = "none", 
        legend.title = element_blank())


# scatterplot - RNA across protocols by sample
ggplot(md_noBlanks, aes(reorder(x = unique_sample_id, rna_ng_ul), y = rna_ng_ul, color = extraction_protocol)) +
  geom_point(size = 0.75) +
  facet_wrap(~sample_type_2, scales = "free") +
  xlab("sample") +
  ylab("RNA yield (ng/µL)") +
  scale_color_manual(values = c("orange", "green3", "blue"), aes(alpha = 0.5)) +
  theme_bw() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position = "top", legend.title = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=3)))


# Mann-Whitney U test for differences in RNA yield across extraction protocols
######################################################################################################

## Paired tests across all data except blanks and non-paired samples ##
wilcox.test(md_noBlanks_unique_sample_id_noCovid_wide$MM_2_rna_ng_ul, md_noBlanks_unique_sample_id_noCovid_wide$PS_rna_ng_ul, paired = TRUE)

wilcox.test(md_noBlanks_unique_sample_id_noCovid_wide$MM_20_rna_ng_ul, md_noBlanks_unique_sample_id_noCovid_wide$PS_rna_ng_ul, paired = TRUE)

wilcox.test(md_noBlanks_unique_sample_id_noCovid_wide$MM_2_rna_ng_ul, md_noBlanks_unique_sample_id_noCovid_wide$MM_20_rna_ng_ul, paired = TRUE)


# Perform tests (high biomass)
wilcox.test(as.numeric(md_noBlanks_noMock_noZeros_highBiomass_2v20$rna_ng_ul) ~ md_noBlanks_noMock_noZeros_highBiomass_2v20$extraction_protocol_short)

wilcox.test(as.numeric(md_noBlanks_noMock_noZeros_highBiomass_2vPS$rna_ng_ul) ~ md_noBlanks_noMock_noZeros_highBiomass_2vPS$extraction_protocol_short)

wilcox.test(as.numeric(md_noBlanks_noMock_noZeros_highBiomass_20vPS$rna_ng_ul) ~ md_noBlanks_noMock_noZeros_highBiomass_20vPS$extraction_protocol_short)


# Perform tests (low biomass)
wilcox.test(as.numeric(md_noBlanks_noMock_noZeros_lowBiomass_2v20$rna_ng_ul) ~ md_noBlanks_noMock_noZeros_lowBiomass_2v20$extraction_protocol_short)

wilcox.test(as.numeric(md_noBlanks_noMock_noZeros_lowBiomass_2vPS$rna_ng_ul) ~ md_noBlanks_noMock_noZeros_lowBiomass_2vPS$extraction_protocol_short)

wilcox.test(as.numeric(md_noBlanks_noMock_noZeros_lowBiomass_20vPS$rna_ng_ul) ~ md_noBlanks_noMock_noZeros_lowBiomass_20vPS$extraction_protocol_short)


# Plot RNA quality
#######################################################################################################

rin <- read_tsv("rna_quality.txt")
levels(as.factor(rin$sample_type))
rin$sample_type <- factor(rin$sample_type, levels = c("single isolate", "yogurt", "sourkraut", "mouse tissue", "mouse feces", "cat feces", "human feces", "human saliva", "bare soil", "rhizosphere soil"))


ggplot(rin, aes(x = sample_type, y = rin)) +
  geom_bar(stat = "identity", aes(fill = extraction_protocol), position = "dodge") +
  xlab("sample type") +
  ylab("\nRIN") +
  theme_bw() +
  theme(text = element_text(size = 14)) +
  scale_fill_manual(values = c("orange", "green4")) +
  theme(axis.text.x = element_text(angle=35, hjust=1),
        legend.position = "top", 
        legend.title = element_blank())



# Plot read counts
######################################################################################################

# Make vector for colors for extraction protocols
extraction_protocol_colors <- c("orange", "green4", "blue")

# boxplot - read counts by extraction protocol - 16S (Figure 1A)
ggplot(md_noBlanks_noMock_noZeros_noArm, aes(x = extraction_protocol_short, y = read_count_16s)) +
  geom_boxplot(outlier.alpha = 0, aes(fill = extraction_protocol_short)) +
  geom_jitter(width = 0.2, size = 0.5) +
  facet_wrap(~sample_type_2, scales = "free_y") +
  stat_summary(fun = mean, geom = "point", color = "black", fill = "red", shape = 21) +
  xlab("extraction protocol") +
  ylab("quality-filtered 16S reads") +
  geom_hline(yintercept = 1E4, linetype = "dashed", color = "gray") +
  scale_fill_manual(values = c("orange", "green4", "blue")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 14),
        strip.text.x = element_text(size = 6),
        strip.background =element_rect(fill="white"))


# boxplot - read counts by extraction protocol - shotgun (Figure 1B)
ggplot(md_noBlanks_noMock_noZeros_noArm, aes(x = extraction_protocol_short, y = read_count_shotgun)) +
  geom_boxplot(outlier.alpha = 0, aes(fill = extraction_protocol_short)) +
  geom_jitter(width = 0.2, size = 0.5) +
  facet_wrap(~sample_type_2, scales = "free_y") +
  stat_summary(fun = mean, geom = "point", color = "black", fill = "red", shape = 21) +
  xlab("extraction protocol") +
  ylab("host- and quality-filtered\nmetagenomic reads") +
  geom_hline(yintercept = 1E6, linetype = "dashed", color = "gray") +
  scale_fill_manual(values = c("orange", "green4", "blue")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 14),
        strip.text.x = element_text(size = 6),
        strip.background =element_rect(fill="white"))


# boxplot - read counts by extraction protocol (log10 reads) - shotgun
ggplot(md_noBlanks_noMock_noZeros, aes(x = extraction_protocol_short, y = log10(read_count_shotgun))) +
  geom_boxplot(outlier.alpha = 0.25) +
  xlab("extraction protocol") +
  ylab("log(host- and quality-filtered\nmetagenomic reads)") +
  theme_bw() +
  theme(text = element_text(size = 14))


# Scatter plot - Shotgun reads per sample - by biomass (Figure S1A)
ggplot(md_noBlanks_noMock_noZeros_noArm, aes(reorder(x = unique_sample_id, read_count_16s), y = read_count_16s, color = extraction_protocol_short)) +
  geom_point(size = 0.75) +
  facet_wrap(~sample_type_2, scales = "free") +
  xlab("sample") +
  ylab("quality-filtered 16S reads") +
  geom_hline(yintercept = 1E4, linetype = "dashed", color = "gray") +
  theme_bw() +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        legend.position = "top", 
        legend.title = element_blank(), 
        text = element_text(size = 14),
        strip.text.x = element_text(size = 6),
        strip.background =element_rect(fill="white")) +
  scale_color_manual(values = c("orange", "green3", "blue"))  +
  guides(colour = guide_legend(override.aes = list(size=3)))


# Scatter plot - Shotgun reads per sample - by biomass (Figure S1B)
ggplot(md_noBlanks_noMock_noZeros_noArm, aes(reorder(x = unique_sample_id, read_count_shotgun), y = read_count_shotgun, color = extraction_protocol_short)) +
  geom_point(size = 0.75) +
  facet_wrap(~sample_type_2, scales = "free") +
  xlab("sample") +
  ylab("host- and quality-filtered\nmetagenomic reads") +
  geom_hline(yintercept = 1E6, linetype = "dashed", color = "gray") +
  theme_bw() +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        legend.position = "top", 
        legend.title = element_blank(), 
        text = element_text(size = 14),
        strip.text.x = element_text(size = 6),
        strip.background =element_rect(fill="white")) +
  scale_color_manual(values = c("orange", "green3", "blue")) +
  guides(colour = guide_legend(override.aes = list(size=3)))


# Mann-Whitney U test for differences in read counts across extraction protocols
######################################################################################################

## 16S data - paired tests across all data except blanks and non-paired samples ##
wilcox.test(md_noBlanks_unique_sample_id_noCovid_wide$MM_2_read_count_16s, md_noBlanks_unique_sample_id_noCovid_wide$PS_read_count_16s, paired = TRUE)

wilcox.test(md_noBlanks_unique_sample_id_noCovid_wide$MM_20_read_count_16s, md_noBlanks_unique_sample_id_noCovid_wide$PS_read_count_16s, paired = TRUE)

wilcox.test(md_noBlanks_unique_sample_id_noCovid_wide$MM_2_read_count_16s, md_noBlanks_unique_sample_id_noCovid_wide$MM_20_read_count_16s, paired = TRUE)


# 16S data, high biomass
wilcox.test(as.numeric(md_noBlanks_noMock_noZeros_highBiomass_2v20$read_count_16s) ~ md_noBlanks_noMock_noZeros_highBiomass_2v20$extraction_protocol_short)

wilcox.test(as.numeric(md_noBlanks_noMock_noZeros_highBiomass_2vPS$read_count_16s) ~ md_noBlanks_noMock_noZeros_highBiomass_2vPS$extraction_protocol_short)

wilcox.test(as.numeric(md_noBlanks_noMock_noZeros_highBiomass_20vPS$read_count_16s) ~ md_noBlanks_noMock_noZeros_highBiomass_20vPS$extraction_protocol_short)


# 16S data, low biomass
wilcox.test(as.numeric(md_noBlanks_noMock_noZeros_lowBiomass_2v20$read_count_16s) ~ md_noBlanks_noMock_noZeros_lowBiomass_2v20$extraction_protocol_short)

wilcox.test(as.numeric(md_noBlanks_noMock_noZeros_lowBiomass_2vPS$read_count_16s) ~ md_noBlanks_noMock_noZeros_lowBiomass_2vPS$extraction_protocol_short)

wilcox.test(as.numeric(md_noBlanks_noMock_noZeros_lowBiomass_20vPS$read_count_16s) ~ md_noBlanks_noMock_noZeros_lowBiomass_20vPS$extraction_protocol_short)


## Shotgun data - paired tests across all data except blanks and non-paired samples ##
wilcox.test(md_noBlanks_unique_sample_id_noCovid_wide$MM_2_read_count_shotgun, md_noBlanks_unique_sample_id_noCovid_wide$PS_read_count_shotgun, paired = TRUE)

wilcox.test(md_noBlanks_unique_sample_id_noCovid_wide$MM_20_read_count_shotgun, md_noBlanks_unique_sample_id_noCovid_wide$PS_read_count_shotgun, paired = TRUE)

wilcox.test(md_noBlanks_unique_sample_id_noCovid_wide$MM_2_read_count_shotgun, md_noBlanks_unique_sample_id_noCovid_wide$MM_20_read_count_shotgun, paired = TRUE)


# Shotgun data, high biomass
wilcox.test(as.numeric(md_noBlanks_noMock_noZeros_highBiomass_2v20$read_count_shotgun) ~ md_noBlanks_noMock_noZeros_highBiomass_2v20$extraction_protocol_short)

wilcox.test(as.numeric(md_noBlanks_noMock_noZeros_highBiomass_2vPS$read_count_shotgun) ~ md_noBlanks_noMock_noZeros_highBiomass_2vPS$extraction_protocol_short)

wilcox.test(as.numeric(md_noBlanks_noMock_noZeros_highBiomass_20vPS$read_count_shotgun) ~ md_noBlanks_noMock_noZeros_highBiomass_20vPS$extraction_protocol_short)


# Shotgun data, low biomass
wilcox.test(as.numeric(md_noBlanks_noMock_noZeros_lowBiomass_2v20$read_count_shotgun) ~ md_noBlanks_noMock_noZeros_lowBiomass_2v20$extraction_protocol_short)

wilcox.test(as.numeric(md_noBlanks_noMock_noZeros_lowBiomass_2vPS$read_count_shotgun) ~ md_noBlanks_noMock_noZeros_lowBiomass_2vPS$extraction_protocol_short)

wilcox.test(as.numeric(md_noBlanks_noMock_noZeros_lowBiomass_20vPS$read_count_shotgun) ~ md_noBlanks_noMock_noZeros_lowBiomass_20vPS$extraction_protocol_short)


# Plot alpha-diversity
######################################################################################################

# boxplot - alpha-diversity by extraction protocol - 16S (Figure 4)
md_noBlanks_noMock_noZeros_alpha_16S_all_groups <- subset (md_noBlanks_noMock_noZeros,
                                    md_noBlanks_noMock_noZeros$sample_type_2 != "doorknob" &
                                    md_noBlanks_noMock_noZeros$sample_type_2 != "mouse jejunum tissue" &
                                    md_noBlanks_noMock_noZeros$sample_type_2 != "mouse feces" &
                                    md_noBlanks_noMock_noZeros$sample_type_2 != "human male urine" &
                                    md_noBlanks_noMock_noZeros$sample_type_2 != "human lower arm")

ggplot(md_noBlanks_noMock_noZeros_alpha_16S_all_groups, aes(x = extraction_protocol_short, y = alpha_16S_rar5k_faithspd)) +
  geom_boxplot(outlier.alpha = 0, aes(fill = extraction_protocol_short)) +
  geom_jitter(width = 0.2, size = 0.5) +
  facet_wrap(~sample_type_2, scales = "free_y") +
  stat_summary(fun = mean, geom = "point", color = "black", fill = "red", shape = 21) +
  xlab("extraction protocol") +
  ylab("Faith's Phylogenetic Diversity") +
  scale_fill_manual(values = c("orange", "green4", "blue")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 14),
        strip.text.x = element_text(size = 6),
        strip.background =element_rect(fill="white"))


# boxplot - alpha-diversity by extraction protocol - shotgun (Figure 4)
md_noBlanks_noMock_noZeros_alpha_shotgun_all_groups <- subset (md_noBlanks_noMock_noZeros,
                            md_noBlanks_noMock_noZeros$sample_type_2 != "human male urine" &
                            md_noBlanks_noMock_noZeros$sample_type_2 != "human lower arm")

ggplot(md_noBlanks_noMock_noZeros_alpha_shotgun_all_groups, aes(x = extraction_protocol_short, y = alpha_shotgun_rar35k_20k_faithspd)) +
  geom_boxplot(outlier.alpha = 0, aes(fill = extraction_protocol_short)) +
  geom_jitter(width = 0.2, size = 0.5) +
  facet_wrap(~sample_type_2, scales = "free_y") +
  stat_summary(fun = mean, geom = "point", color = "black", fill = "red", shape = 21) +
  xlab("extraction protocol") +
  ylab("Faith's Phylogenetic Diversity") +
  scale_fill_manual(values = c("orange", "green4", "blue")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 14),
        strip.text.x = element_text(size = 6),
        strip.background =element_rect(fill="white"))


# Scatter plot - alpha-diversity by sample facets sample_type_2 - 16S (Figure S7A)
md_noBlanks_noMock_noZeros_alpha_16S_all_groups_noMissing <- subset(md_noBlanks_noMock_noZeros_alpha_16S_all_groups, md_noBlanks_noMock_noZeros_alpha_16S_all_groups$alpha_16S_rar5k_faithspd >= 0)

ggplot(md_noBlanks_noMock_noZeros_alpha_16S_all_groups_noMissing, aes(reorder(x = unique_sample_id, alpha_16S_rar5k_faithspd), y = alpha_16S_rar5k_faithspd, color = extraction_protocol_short)) +
  geom_point(size = 0.75) +
  facet_wrap(~sample_type_2, scales = "free") +
  xlab("sample") +
  ylab("Faith's Phylogenetic Diversity") +
  theme_bw() +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        legend.position = "top", 
        legend.title = element_blank(), 
        text = element_text(size = 14),
        strip.text.x = element_text(size = 6),
        strip.background =element_rect(fill="white")) +
  scale_color_manual(values = c("orange", "green3", "blue"))  +
  guides(colour = guide_legend(override.aes = list(size=3)))


# Scatter plot - alpha-diversity by sample facets sample_type_2 - shotgun (Figure S7B)
md_noBlanks_noMock_noZeros_alpha_shotgun_all_groups_noMissing <- subset(md_noBlanks_noMock_noZeros_alpha_shotgun_all_groups, md_noBlanks_noMock_noZeros_alpha_shotgun_all_groups$alpha_shotgun_rar35k_20k_faithspd >= 0)

ggplot(md_noBlanks_noMock_noZeros_alpha_shotgun_all_groups_noMissing, aes(reorder(x = unique_sample_id, alpha_shotgun_rar35k_20k_faithspd), y = alpha_shotgun_rar35k_20k_faithspd, color = extraction_protocol_short), drop = TRUE) +
  geom_point(size = 0.75) +
  facet_wrap(~sample_type_2, scales = "free", drop = TRUE) +
  xlab("sample") +
  ylab("Faith's Phylogenetic Diversity") +
  theme_bw() +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        legend.position = "top", 
        legend.title = element_blank(), 
        text = element_text(size = 14),
        strip.text.x = element_text(size = 6),
        strip.background =element_rect(fill="white")) +
  scale_color_manual(values = c("orange", "green3", "blue"))  +
  guides(colour = guide_legend(override.aes = list(size=3)))


# Mann-Whitney U test for differences in alpha-diversity across extraction protocols
######################################################################################################

## 16S data - paired tests across all data except blanks and non-paired samples ##
wilcox.test(md_noBlanks_unique_sample_id_noCovid_wide$MM_2_alpha_16S_rar5k_faithspd, md_noBlanks_unique_sample_id_noCovid_wide$PS_alpha_16S_rar5k_faithspd, paired = TRUE)

wilcox.test(md_noBlanks_unique_sample_id_noCovid_wide$MM_20_alpha_16S_rar5k_faithspd, md_noBlanks_unique_sample_id_noCovid_wide$PS_alpha_16S_rar5k_faithspd, paired = TRUE)

wilcox.test(md_noBlanks_unique_sample_id_noCovid_wide$MM_2_alpha_16S_rar5k_faithspd, md_noBlanks_unique_sample_id_noCovid_wide$MM_20_alpha_16S_rar5k_faithspd, paired = TRUE)


# 16S data, high biomass
wilcox.test(as.numeric(md_noBlanks_noMock_noZeros_highBiomass_2v20$alpha_16S_rar5k_faithspd) ~ md_noBlanks_noMock_noZeros_highBiomass_2v20$extraction_protocol_short)

wilcox.test(as.numeric(md_noBlanks_noMock_noZeros_highBiomass_2vPS$alpha_16S_rar5k_faithspd) ~ md_noBlanks_noMock_noZeros_highBiomass_2vPS$extraction_protocol_short)

wilcox.test(as.numeric(md_noBlanks_noMock_noZeros_highBiomass_20vPS$alpha_16S_rar5k_faithspd) ~ md_noBlanks_noMock_noZeros_highBiomass_20vPS$extraction_protocol_short)


# 16S data, low biomass
wilcox.test(as.numeric(md_noBlanks_noMock_noZeros_lowBiomass_2v20$alpha_16S_rar5k_faithspd) ~ md_noBlanks_noMock_noZeros_lowBiomass_2v20$extraction_protocol_short)

wilcox.test(as.numeric(md_noBlanks_noMock_noZeros_lowBiomass_2vPS$alpha_16S_rar5k_faithspd) ~ md_noBlanks_noMock_noZeros_lowBiomass_2vPS$extraction_protocol_short)

wilcox.test(as.numeric(md_noBlanks_noMock_noZeros_lowBiomass_20vPS$alpha_16S_rar5k_faithspd) ~ md_noBlanks_noMock_noZeros_lowBiomass_20vPS$extraction_protocol_short)


# Shotgun data, high biomass
wilcox.test(as.numeric(md_noBlanks_noMock_noZeros_highBiomass_2v20$alpha_shotgun_rar35k_20k_faithspd) ~ md_noBlanks_noMock_noZeros_highBiomass_2v20$extraction_protocol_short)

wilcox.test(as.numeric(md_noBlanks_noMock_noZeros_highBiomass_2vPS$alpha_shotgun_rar35k_20k_faithspd) ~ md_noBlanks_noMock_noZeros_highBiomass_2vPS$extraction_protocol_short)

wilcox.test(as.numeric(md_noBlanks_noMock_noZeros_highBiomass_20vPS$alpha_shotgun_rar35k_20k_faithspd) ~ md_noBlanks_noMock_noZeros_highBiomass_20vPS$extraction_protocol_short)


# Shotgun data, low biomass
wilcox.test(as.numeric(md_noBlanks_noMock_noZeros_lowBiomass_2v20$alpha_shotgun_rar35k_20k_faithspd) ~ md_noBlanks_noMock_noZeros_lowBiomass_2v20$extraction_protocol_short)

wilcox.test(as.numeric(md_noBlanks_noMock_noZeros_lowBiomass_2vPS$alpha_shotgun_rar35k_20k_faithspd) ~ md_noBlanks_noMock_noZeros_lowBiomass_2vPS$extraction_protocol_short)

wilcox.test(as.numeric(md_noBlanks_noMock_noZeros_lowBiomass_20vPS$alpha_shotgun_rar35k_20k_faithspd) ~ md_noBlanks_noMock_noZeros_lowBiomass_20vPS$extraction_protocol_short)


# Survival analysis (Kaplan-Meier estimator) - 16S
######################################################################################################
install.packages('survminer')
library('survival')
library('survminer')

fit_16S <- survfit(Surv(read_count_16s) ~ extraction_protocol_short, data = md)
ggsurvplot(fit_16S, data = md,
           conf.int = TRUE,
           ggtheme = theme_bw(),
           xlim = c(0,20000),
           xlab = "quality-filtered 16S reads",
           ylab = "proportion of samples",
           break.time.by = 5000)


# Binary biomass - High
fit_16S_high <- survfit(Surv(read_count_16s) ~ extraction_protocol_short, data = md_2_high)
ggsurvplot(fit_16S_high, data = md_2_high,
           conf.int = TRUE,
           ggtheme = theme_bw(),
           xlim = c(0,20000),
           xlab = "quality-filtered 16S reads",
           ylab = "proportion of samples",
           break.time.by = 5000)


# Binary biomass - low
fit_16S_low <- survfit(Surv(read_count_16s) ~ extraction_protocol_short, data = md_2_low)
ggsurvplot(fit_16S_low, data = md_2_low,
           conf.int = TRUE,
           ggtheme = theme_bw(),
           xlim = c(0,20000),
           xlab = "quality-filtered 16S reads",
           ylab = "proportion of samples",
           break.time.by = 5000)


# Three biomass groups - High
fit_16S_high2 <- survfit(Surv(read_count_16s) ~ extraction_protocol_short, data = md_2_high2)
ggsurvplot(fit_16S_high2, data = md_2_high2,
           conf.int = TRUE,
           ggtheme = theme_bw(),
           xlim = c(0,20000),
           xlab = "quality-filtered 16S reads",
           ylab = "proportion of samples",
           break.time.by = 5000)


# Three biomass groups - Medium
fit_16S_med <- survfit(Surv(read_count_16s) ~ extraction_protocol_short, data = md_2_med)
ggsurvplot(fit_16S_med, data = md_2_med,
           conf.int = TRUE,
           ggtheme = theme_bw(),
           xlim = c(0,20000),
           xlab = "quality-filtered 16S reads",
           ylab = "proportion of samples",
           break.time.by = 5000)


# Three biomass groups - Low
fit_16S_low2 <- survfit(Surv(read_count_16s) ~ extraction_protocol_short, data = md_2_low2)
ggsurvplot(fit_16S_low2, data = md_2_low2,
           conf.int = TRUE,
           ggtheme = theme_bw(),
           xlim = c(0,20000),
           xlab = "quality-filtered 16S reads",
           ylab = "proportion of samples",
           break.time.by = 5000)


# Survival analysis (Kaplan-Meier estimator) - shotgun
#########################################################################################
fit_shotgun <- survfit(Surv(read_count_shotgun) ~ extraction_protocol_short, data = md_2)
ggsurvplot(fit_shotgun, data = md_2,
           conf.int = TRUE,
           ggtheme = theme_bw(),
           xlim = c(0,2000000),
           xlab = "human- and quality-filtered\nmetagenomic reads",
           ylab = "proportion of samples",
           break.time.by = 500000)


# Binary biomass groups - High
fit_shotgun_high <- survfit(Surv(read_count_shotgun) ~ extraction_protocol_short, data = md_2_high)
ggsurvplot(fit_shotgun_high, data = md_2_high,
           conf.int = TRUE,
           ggtheme = theme_bw(),
           xlim = c(0,2000000),
           xlab = "human- and quality-filtered\nmetagenomic reads",
           ylab = "proportion of samples",
           break.time.by = 500000)


# Binary biomass groups - Low
fit_shotgun_low <- survfit(Surv(read_count_shotgun) ~ extraction_protocol_short, data = md_2_low)
ggsurvplot(fit_shotgun_low, data = md_2_low,
           conf.int = TRUE,
           ggtheme = theme_bw(),
           xlim = c(0,2000000),
           xlab = "human- and quality-filtered\nmetagenomic reads",
           ylab = "proportion of samples",
           break.time.by = 500000)


# Three biomass groups - High
fit_shotgun_high2 <- survfit(Surv(read_count_shotgun) ~ extraction_protocol_short, data = md_2_high2)
ggsurvplot(fit_shotgun_high2, data = md_2_high2,
           conf.int = TRUE,
           ggtheme = theme_bw(),
           xlim = c(0,2000000),
           xlab = "human- and quality-filtered\nmetagenomic reads",
           ylab = "proportion of samples",
           break.time.by = 500000)


# Three biomass groups - Medium
fit_shotgun_med <- survfit(Surv(read_count_shotgun) ~ extraction_protocol_short, data = md_2_med)
ggsurvplot(fit_shotgun_med, data = md_2_med,
           conf.int = TRUE,
           ggtheme = theme_bw(),
           xlim = c(0,2000000),
           xlab = "human- and quality-filtered\nmetagenomic reads",
           ylab = "proportion of samples",
           break.time.by = 500000)


# Three biomass groups - Low
fit_shotgun_low2 <- survfit(Surv(read_count_shotgun) ~ extraction_protocol_short, data = md_2_low2)
ggsurvplot(fit_shotgun_low2, data = md_2_low2,
           conf.int = TRUE,
           ggtheme = theme_bw(),
           xlim = c(0,2000000),
           xlab = "human- and quality-filtered\nmetagenomic reads",
           ylab = "proportion of samples",
           break.time.by = 500000)


# Well-to-well contamination - heatmap visualization
#########################################################################################
library(viridis)

# Read in data
magmax_2min_w2w <- read.csv("well_to_well_ggplot2_magmax2min.txt", sep = "\t")
magmax_20min_w2w <- read.csv("well_to_well_ggplot2_magmax20min.txt", sep = "\t")
powersoil_w2w <- read.csv("well_to_well_ggplot2_powersoil.txt", sep = "\t")
all_w2w <- read.csv("well_to_well_ggplot2_allProtocols.txt", sep = "\t")


# Plot data for each protocol - log10 transformed
ggplot(magmax_2min_w2w, aes(x = column, y = reorder(row, desc(row)), fill = plasmid_reads_log10)) +
  geom_tile() +
  scale_fill_viridis(name = "log10(plasmid reads)", discrete = F) +
  scale_x_discrete(limits = c(1:12)) +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5)) +
  ggtitle("MagMAX (2-min.)")


ggplot(magmax_20min_w2w, aes(x = column, y = reorder(row, desc(row)), fill = plasmid_reads_log10)) +
  geom_tile() +
  scale_fill_viridis(name = "log10(plasmid reads)", discrete = F) +
  scale_x_discrete(limits = c(1:12)) +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5)) +
  ggtitle("MagMAX (20-min.)")


ggplot(powersoil_w2w, aes(x = column, y = reorder(row, desc(row)), fill = plasmid_reads_log10)) +
  geom_tile() +
  scale_fill_viridis(name = "log10(plasmid reads)", discrete = F) +
  scale_x_discrete(limits = c(1:12)) +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5)) +
  ggtitle("PowerSoil")



# Plot data for each protocol - raw read counts
ggplot(magmax_2min_w2w, aes(x = column, y = reorder(row, desc(row)), fill = plasmid_reads)) +
  geom_tile() +
  scale_fill_viridis(name = "plasmid reads", discrete = F) +
  scale_x_discrete(limits = c(1:12)) +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5)) +
  ggtitle("MagMAX (2-min.)")


ggplot(magmax_20min_w2w, aes(x = column, y = reorder(row, desc(row)), fill = plasmid_reads)) +
  geom_tile() +
  scale_fill_viridis(name = "plasmid reads", discrete = F) +
  scale_x_discrete(limits = c(1:12)) +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5)) +
  ggtitle("MagMAX (20-min.)")


ggplot(powersoil_w2w, aes(x = column, y = reorder(row, desc(row)), fill = plasmid_reads)) +
  geom_tile() +
  scale_fill_viridis(name = "plasmid reads", discrete = F) +
  scale_x_discrete(limits = c(1:12)) +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5)) +
  ggtitle("PowerSoil")


# Plot data for each protocol - percent plasmid reads
ggplot(magmax_2min_w2w, aes(x = column, y = reorder(row, desc(row)), fill = plasmid_reads_percent)) +
  geom_tile() +
  scale_fill_viridis(name = "% plasmid reads", discrete = F) +
  scale_x_discrete(limits = c(1:12)) +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5)) +
  ggtitle("MagMAX (2-min.)")


ggplot(magmax_20min_w2w, aes(x = column, y = reorder(row, desc(row)), fill = plasmid_reads_percent)) +
  geom_tile() +
  scale_fill_viridis(name = "% plasmid reads", discrete = F) +
  scale_x_discrete(limits = c(1:12)) +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5)) +
  ggtitle("MagMAX (20-min.)")


ggplot(powersoil_w2w, aes(x = column, y = reorder(row, desc(row)), fill = plasmid_reads_percent)) +
  geom_tile() +
  scale_fill_viridis(name = "% plasmid reads", discrete = F) +
  scale_x_discrete(limits = c(1:12)) +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5)) +
  ggtitle("PowerSoil")


ggplot(all_w2w, aes(x = column, y = reorder(row, desc(row)), fill = plasmid_reads_percent)) +
  geom_tile() +
  facet_wrap(~extraction_protocol) +
  scale_fill_viridis(name = "% plasmid reads", discrete = F) +
  scale_x_discrete(limits = c(1:12)) +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(hjust = 0.5))


# Well-to-well contamination - boxplot
#########################################################################################
library(tidyverse)
install.packages("PMCMRplus")
library(PMCMRplus)

all_w2w <- read.csv("well_to_well_ggplot2_allProtocols.txt", sep = "\t")
all_w2w_sinkOnly <- subset(all_w2w, all_w2w$well_type == "sink")
  
# boxplot
ggplot(all_w2w_sinkOnly, aes(x = extraction_protocol, y = plasmid_reads_log10)) +
  geom_boxplot(outlier.alpha = 0, aes(fill = extraction_protocol)) +
  geom_jitter(width = 0.2, size = 0.5) +
  #facet_wrap(~sample_type_2, scales = "free_y") +
  stat_summary(fun = mean, geom = "point", color = "black", fill = "red", shape = 21) +
  xlab("extraction protocol") +
  ylab("log(spike-in read count\nin sink well)") +
  ylim(c(-0.01,5)) +
  scale_fill_manual(values = c("orange", "green4", "blue")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_blank(),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 14),
        strip.text.x = element_text(size = 6),
        strip.background =element_rect(fill="white"))

kruskal.test(all_w2w_sinkOnly$plasmid_reads ~ all_w2w_sinkOnly$extraction_protocol)
posthoc.kruskal.dunn.test(all_w2w_sinkOnly$plasmid_reads, 
                          as.factor(all_w2w_sinkOnly$extraction_protocol), 
                          p.adjust.method = "hochberg")


# Technical replicate distance analysis
#########################################################################################

colors = c("#d6604d","#0571b0","#92c5de")

# Jaccard
tech_rep_dist_jaccard <- read_tsv("technical_replicate_distances_jaccard.tsv")

levels(as.factor(tech_rep_dist_jaccard$sample_type_2))
tech_rep_dist_jaccard$sample_type_2 <- factor(tech_rep_dist_jaccard$sample_type_2, levels = c("floor", "sourkraut", "yogurt", "seawater", "freshwater", "bare soil", "rhizosphere soil", "cat feces", "human feces", "human female urine", "human foot", "human armpit", "human forehead", "human nares", "human throat", "human saliva", "human oral saline rinse", "human oral saline rinse in VTM", "single strain"))

ggplot(tech_rep_dist_jaccard, aes(x = sample1_extraction_protocol_short, y = value)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(aes(color = biomass_sample_long), width = 0.2, size = 0.3) +
  facet_wrap(~sample_type_2) +
  ylab("Within sample\nJaccard distance") +
  scale_color_manual(values = c("orange2", "blue", "green4")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "top", 
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))
  
  
# RPCA
tech_rep_dist_rpca <- read_tsv("technical_replicate_distances_rpca.tsv")

levels(as.factor(tech_rep_dist_rpca$sample_type_2))
tech_rep_dist_rpca$sample_type_2 <- factor(tech_rep_dist_rpca$sample_type_2, levels = c("floor", "sourkraut", "yogurt", "seawater", "freshwater", "bare soil", "rhizosphere soil", "cat feces", "human feces", "human female urine", "human foot", "human armpit", "human forehead", "human nares", "human throat", "human saliva", "human oral saline rinse", "human oral saline rinse in VTM", "single strain", "two strains"))

ggplot(tech_rep_dist_rpca, aes(x = sample1_extraction_protocol_short, y = value)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(aes(color = biomass_sample_long), width = 0.2, size = 0.3) +
  facet_wrap(~sample_type_2) +
  ylab("Within sample\nRPCA distance") +
  scale_color_manual(values = c("orange2", "blue", "green4")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "top", 
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))


# unweighted UniFrac
tech_rep_dist_unifrac <- read_tsv("technical_replicate_distances_unweighted_unifrac.tsv")

levels(as.factor(tech_rep_dist_unifrac$sample_type_2))
tech_rep_dist_unifrac$sample_type_2 <- factor(tech_rep_dist_unifrac$sample_type_2, levels = c("floor", "sourkraut", "yogurt", "seawater", "freshwater", "bare soil", "rhizosphere soil", "cat feces", "human feces", "human female urine", "human foot", "human armpit", "human forehead", "human nares", "human throat", "human saliva", "human oral saline rinse", "human oral saline rinse in VTM", "single strain", "two strains"))

ggplot(tech_rep_dist_unifrac, aes(x = sample1_extraction_protocol_short, y = value)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(aes(color = biomass_sample_long), width = 0.2, size = 0.3) +
  facet_wrap(~sample_type_2) +
  ylab("Within sample\nunweighted UniFrac distance") +
  scale_color_manual(values = c("orange2", "blue", "green4", "black")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "top", 
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))


# weighted UniFrac
tech_rep_dist_wunifrac <- read_tsv("technical_replicate_distances_weighted_unifrac.tsv")

levels(as.factor(tech_rep_dist_wunifrac$sample_type_2))
tech_rep_dist_wunifrac$sample_type_2 <- factor(tech_rep_dist_wunifrac$sample_type_2, levels = c("floor", "sourkraut", "yogurt", "seawater", "freshwater", "bare soil", "rhizosphere soil", "cat feces", "human feces", "human female urine", "human foot", "human armpit", "human forehead", "human nares", "human throat", "human saliva", "human oral saline rinse", "human oral saline rinse in VTM", "single strain", "two strains"))

ggplot(tech_rep_dist_wunifrac, aes(x = sample1_extraction_protocol_short, y = value)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(aes(color = biomass_sample_long), width = 0.2, size = 0.3) +
  facet_wrap(~sample_type_2) +
  ylab("Within sample\nweighted UniFrac distance") +
  scale_color_manual(values = c("orange2", "blue", "green4", "black")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "top", 
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))


# Shotgun

# Jaccard
tech_rep_dist_shotgun_high_jaccard <- read_tsv("technical_replicate_distances_shotgun_high_biomass_jaccard.tsv")
tech_rep_dist_shotgun_low_jaccard <- read_tsv("technical_replicate_distances_shotgun_low_biomass_jaccard.tsv")
tech_rep_dist_shotgun_jaccard <- rbind(tech_rep_dist_shotgun_high_jaccard, tech_rep_dist_shotgun_low_jaccard)

levels(as.factor(tech_rep_dist_shotgun_jaccard$sample_type_2))
tech_rep_dist_shotgun_jaccard$sample_type_2 <- factor(tech_rep_dist_shotgun_jaccard$sample_type_2, levels = c("doorknob", "floor", "sourkraut", "yogurt", "seawater", "freshwater", "rhizosphere soil", "mouse jejunum tissue", "mouse feces", "cat feces", "human feces", "human female urine", "human foot", "human armpit", "human forehead", "human nares", "human throat", "human saliva", "human oral saline rinse", "human oral saline rinse in VTM"))

ggplot(tech_rep_dist_shotgun_jaccard, aes(x = sample1_extraction_protocol_short, y = value)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(aes(color = biomass_sample_long), width = 0.2, size = 0.3) +
  facet_wrap(~sample_type_2) +
  ylab("Within sample\nJaccard distance") +
  scale_color_manual(values = c("orange2", "blue", "green4")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "top", 
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))


# RPCA
tech_rep_dist_shotgun_high_rpca <- read_tsv("technical_replicate_distances_shotgun_high_biomass_rpca.tsv")
tech_rep_dist_shotgun_low_rpca <- read_tsv("technical_replicate_distances_shotgun_low_biomass_rpca.tsv")
tech_rep_dist_shotgun_rpca <- rbind(tech_rep_dist_shotgun_high_rpca, tech_rep_dist_shotgun_low_rpca)

levels(as.factor(tech_rep_dist_shotgun_rpca$sample_type_2))
tech_rep_dist_shotgun_rpca$sample_type_2 <- factor(tech_rep_dist_shotgun_rpca$sample_type_2, levels = c("doorknob", "floor", "sourkraut", "yogurt", "seawater", "freshwater", "rhizosphere soil", "mouse jejunum tissue", "mouse feces", "cat feces", "human feces", "human female urine", "human foot", "human armpit", "human forehead", "human nares", "human throat", "human saliva", "human oral saline rinse", "human oral saline rinse in VTM"))

ggplot(tech_rep_dist_shotgun_rpca, aes(x = sample1_extraction_protocol_short, y = value)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(aes(color = biomass_sample_long), width = 0.2, size = 0.3) +
  facet_wrap(~sample_type_2) +
  ylab("Within sample\nRPCA distance") +
  scale_color_manual(values = c("orange2", "blue", "green4")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "top", 
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))


# unweighted UniFrac
tech_rep_dist_shotgun_high_unifrac <- read_tsv("technical_replicate_distances_shotgun_high_biomass_unweighted_unifrac.tsv")
tech_rep_dist_shotgun_low_unifrac <- read_tsv("technical_replicate_distances_shotgun_low_biomass_unweighted_unifrac.tsv")
tech_rep_dist_shotgun_unifrac <- rbind(tech_rep_dist_shotgun_high_unifrac, tech_rep_dist_shotgun_low_unifrac)

levels(as.factor(tech_rep_dist_shotgun_unifrac$sample_type_2))
tech_rep_dist_shotgun_unifrac$sample_type_2 <- factor(tech_rep_dist_shotgun_unifrac$sample_type_2, levels = c("doorknob", "floor", "sourkraut", "yogurt", "seawater", "freshwater", "rhizosphere soil", "mouse jejunum tissue", "mouse feces", "cat feces", "human feces", "human female urine", "human foot", "human armpit", "human forehead", "human nares", "human throat", "human saliva", "human oral saline rinse", "human oral saline rinse in VTM"))

ggplot(tech_rep_dist_shotgun_unifrac, aes(x = sample1_extraction_protocol_short, y = value)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(aes(color = biomass_sample_long), width = 0.2, size = 0.3) +
  facet_wrap(~sample_type_2) +
  ylab("Within sample\nunweighted UniFrac distance") +
  scale_color_manual(values = c("orange2", "blue", "green4")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "top", 
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))


# weighted UniFrac
tech_rep_dist_shotgun_high_wunifrac <- read_tsv("technical_replicate_distances_shotgun_high_biomass_weighted_unifrac.tsv")
tech_rep_dist_shotgun_low_wunifrac <- read_tsv("technical_replicate_distances_shotgun_low_biomass_weighted_unifrac.tsv")
tech_rep_dist_shotgun_wunifrac <- rbind(tech_rep_dist_shotgun_high_wunifrac, tech_rep_dist_shotgun_low_wunifrac)

levels(as.factor(tech_rep_dist_shotgun_wunifrac$sample_type_2))
tech_rep_dist_shotgun_wunifrac$sample_type_2 <- factor(tech_rep_dist_shotgun_wunifrac$sample_type_2, levels = c("doorknob", "floor", "sourkraut", "yogurt", "seawater", "freshwater", "rhizosphere soil", "mouse jejunum tissue", "mouse feces", "cat feces", "human feces", "human female urine", "human foot", "human armpit", "human forehead", "human nares", "human throat", "human saliva", "human oral saline rinse", "human oral saline rinse in VTM"))

ggplot(tech_rep_dist_shotgun_wunifrac, aes(x = sample1_extraction_protocol_short, y = value)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(aes(color = biomass_sample_long), width = 0.2, size = 0.3) +
  facet_wrap(~sample_type_2) +
  ylab("Within sample\nweighted UniFrac distance") +
  scale_color_manual(values = c("orange2", "blue", "green4")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "top", 
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 14)) +
  guides(colour = guide_legend(override.aes = list(size=3)))
