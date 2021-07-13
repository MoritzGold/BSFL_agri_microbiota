# Analayses of DNA-based amplicon sequencing results for the manuscript 
# "Effect of rearing system and microbial inoculants on black soldier fly 
# growth and microbiota with two agro-food byproducts".

# The data set also includes samples previously published in Gold et al. (2020): 
# Identification of bacteria in two food waste black soldier fly larvae rearing residues
# published in Frontiers in Microbiology

# Libraries and functions -------------------------------------------------

# clean/reset environment 

rm(list=ls()) 

# R and Bioconductor libraries 

require("ggplot2")
require("vegan")
#require("knitr")
require("phyloseq")
require("microbiome")
#require("ampvis2")
#require("betapart")
#require("dplyr")
#require("magrittr")
#require("doBy")
#require("plotly")
#require("cooccur")
#require("scales")
#require("picante")
library(tidyverse)
#library(metacoder)
library(ampvis2)
library(openxlsx)
library(UpSetR)
library(ComplexHeatmap)
library(ggplot2)
library(ggord)
library(corrplot)
library(usdm)
library(conflicted)
library(gridExtra)
library(metagMisc)
library(cowplot)
library(gridExtra)
library(gridExtra)
library(ggpubr)
library(factoextra)
library(cluster)
library(fastcluster)
library(usedist)
library(clValid)
library(factoextra)
library(immunarch)
library(cluster)
library(pvclust)
library(fpc)
library(ggalt)
library(ungeviz)


conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("intersect", "dplyr")
conflict_prefer("lag", "dplyr")
conflict_prefer("mutate","dplyr")
conflict_prefer("summarise","dplyr")

# Functions

source("code/Rfunctions.R")


# Data import in R --------------------------------------------------------

zotufile    <- "data/p529_run190225_run200606_16S_ZOTU_Count_Sintax.txt"

mapfile     <- "data/MapFile.txt"

treefile    <- "data/p529_run190225_run200606_16S_ZOTU_MSA.tre"


d <- import_qiime(otufilename = zotufile, mapfilename = mapfile, treefilename = treefile)


# Remove samples from zotufile 

  # some samples in run had clogging filters and are therefore not representative, these are samples
  # from Gold et al. (2020): Identification of bacteria in two food waste black soldier fly larvae rearing residues
  # published in Frontiers in Microbiology
  
  d <- subset_samples(d, !SampleID %in% c("12SL2","12SR2","9SL1","9SL2","9SR1","9SR2"))
  d <- prune_taxa(taxa_sums(d) > 0, d)

  
# Remove rare taxa < 10 reads ---------------------------------------------------------
  
  d <- prune_taxa(taxa_sums(d) >= 10, d)      
  
# save runs as individual sets ------------------------------------
  
  # run 1 (Gold et al., 2020)
  
  d_1 <- subset_samples(d, run == "1")
  d_1 <- prune_taxa(taxa_sums(d_1) > 0, d_1)
  
  # run 2 (this manuscript)
  
  d_2 <- subset_samples(d, run == "2a"| run == "2b")
  d_2 <- prune_taxa(taxa_sums(d_2) > 0, d_2)
                    
# Look at controls -------------------------------------------------  
  
  # run 2
  
  
    # positives (mock) - Zymo6311
  
    d_2_mocks <- subset_samples(d_2, type == "Mock341")
    d_2_mocks <- prune_taxa(taxa_sums(d_2_mocks) > 0, d_2_mocks)
    
    sort(taxa_sums(d_2_mocks), TRUE) 
    
    tax_table(d_2_mocks)[c("OTU3","OTU90","OTU72","OTU204","OTU46","OTU11","OTU505","OTU164","OTU65","OTU524","OTU522")]
    
    #1 Listeria monocytogenes: 95.9 = #1 OTU3
    #2 Pseudomonas aeruginosa: 2.8 = #2 OTU 90
    #3 Bacillus subtilis: 1.2 = #3 OTU72
    #4 Saccharomyces cerevisiae: NA
    #5 Escherichia coli: 0.069 = # 4 OTU46
    #6 Salmonella enterica: 0.07 = #5 OTU11
    #7 Lactobacillus fermentum: 0.012 = #6 OTU 164
    #8 Enterococcus faecalis: 0.00067
    #9 Cryptococcus neoformans: NA
    #10 Staphylococcus aureus: 0.0001
  
    d_2_mocks_norm <- transform_sample_counts(d_2_mocks, function(x) x / sum(x))
  
    d_2_mocks_norm <- filter_taxa(d_2_mocks_norm, function(x) mean(x) > 0.001, TRUE)
    
    otu_table(d_2_mocks_norm)
  
    
    # negatives (water)
  
    d_2_controls <- subset_samples(d_2, SampleID =="Blank1a" | SampleID =="Blank2a" | SampleID =="KitBlank")
    d_2_controls <- prune_taxa(taxa_sums(d_2_controls) > 0, d_2_controls)
    
    summary(sample_sums(d_2_controls)) 
    
    sample_sums(d_2_controls)
    
    # very low number of reads in the controls
 
    # Samples with read counts below reads in controls
    
    d_2_low <- prune_samples(sample_sums(d_2)<=150, d_2)
    
    sample_sums(d_2_low)
  
    #MG004: wine residue, experiment 4, closed, replicate 2
    #MG056: wine waste, replicate 1
    #MG057: wine waste, replicate 2
    
    # no results in both runs, propably due to low number of DNA
  
  
# Remove samples ------------------- 
  
  # Larval samples - possible cross-contamination
  
  d <- subset_samples(d, !(run == "1" & type == "larvae"))  
  
  # chicken feed substrate - not relevant
  
  d <- subset_samples(d, !waste == "CF")
  
  # Residue, canteen waste, day 12, replicate 3 - deviates very much from the other three samples
  
  d <- subset_samples(d, !(type=="residue"&type_detail=="residue"&waste=="C"&Day=="12"&replicate=="3"))
  
  # Residue, household waste, control without BSFL - only replictae which had a different feeding rate in comparison to all other containers, not representative
  
  d <- subset_samples(d, !(type_detail=="No_larvae"&waste=="H"))
  
  # Samples with less than 2000 reads and update reads table, only including taxa with > 0 reads
  
  d = prune_samples(sample_sums(d)>=2000, d)
  
  # phylum cynaobacterium (=plant chloroplasts)
  
  d = subset_taxa(d, !Phylum=="Cyanobacteria")
  
  # family Mitochondria (=plant mitochondria)
  
  d = subset_taxa(d, !Family=="Mitochondria")

  # Update reads table, only including taxa with > 0 reads
  
  d <- prune_taxa(taxa_sums(d) > 0, d)  
  
  # Remove controls
  
  d <-  subset_samples(d, !(type =="lib_prep_blank"))
  
  d <-  subset_samples(d, !(type =="extraction_blank"))

  # Remove Mocks
  
  d <-  subset_samples(d, !(type =="Mock319"))
  
  d <-  subset_samples(d, !(type =="Mock341"))
  
  
  # update separate runs
  
  
  # run 2
  
  d_2 <- subset_samples(d, run == "2a"| run == "2b")
  d_2 <- prune_taxa(taxa_sums(d_2) > 0, d_2)

  
     
# Sample summary -------------------   
  
  
  # run 2

  # Number of read counts per taxon
  
  taxa_sums(d_2)
  
  # Number of read counts per sample
  
  sample_sums(d_2)
  
  # Descriptive statistics on read counts per sample
  
  summary(sample_sums(d_2))
  
  # Summary of phyloseq object
  
  summarize_phyloseq(d_2)
 
  d_2
  
  
# Rarefaction curves ------------------------------------------------------------
  
  
  # run 2
   
  rare_d_2 <- rarecurve(t(otu_table(d_2)), step=50, cex=0.5,label=FALSE)
  

# alpha diversity ------------
  
  
  # subset samples of run 2 
  
  alpha_d_2 <- plot_richness(d_2, measures=c("Shannon","Chao1"))
  
  
  alpha_d_2 <- alpha_d_2$data 
  

  # calculate mean and sd (n>2)
  
  alpha_d_2_mean <-
    
    alpha_d_2 %>% 
    
    group_by(waste,type,experiment,system,variable,others,type_detail) %>% 
    
    summarise(n=n(),
              mean=mean(value),
              sd=sd(value))
  
  write.xlsx(alpha_d_2_mean , file="output/alpha_d_2_mean.xlsx")
  
  
  # introduce function to ensure that adequate number of decimals are displayed
  
  f <- function(x){
    format(round(x, 0), nsmall=0)
  }
  

  # plot alpha diversity metrics of experiment 4 (WWP, TP, DFW without inoculant) -----
  
  # subset residue samples
  
  alpha_d_2_exp4_residue <-
    
    alpha_d_2 %>% 
    
    filter(type== "residue" & experiment== "4") %>% 
    
    mutate(waste=case_when(waste=="tomato"~"TP",
                           waste=="foodwaste"~"DFW",
                           waste=="wine"~"WP",
                           TRUE~""))
  
  # subset larvae samples
  
  alpha_d_2_exp4_larvae <- 
    
    alpha_d_2 %>% 
    
    filter(type== "larvae" & experiment== "4") %>% 
    
    mutate(waste=case_when(waste=="tomato"~"TP",
                           waste=="foodwaste"~"DFW",
                           waste=="wine"~"WP",
                           TRUE~""))
  
  # plot residue alpha diversity of residue samples
  
  alpha_d_2_mean %>% 
    
    filter(type== "residue" & experiment== "4") %>% 
    
    mutate(waste=case_when(waste=="tomato"~"TP",
                           waste=="foodwaste"~"DFW",
                           waste=="wine"~"WP",
                           TRUE~"")) %>% 
    
    ggplot(aes(waste,mean)) +
    
    geom_hpline(data=alpha_d_2_exp4_residue,
                aes(waste, value, shape=system),
                position = position_dodge(0.6),
                width = 0.2, size = 0.7,
                stat = "summary") +
    
    geom_jitter(data=alpha_d_2_exp4_residue, 
                aes(waste, value, shape = system), 
                position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.6),
                size = 2.5) +
    
    stat_summary(data=alpha_d_2_exp4_residue,
                 aes(waste,value, shape = system),
                 fun.data="mean_sdl",  fun.args = list(mult=1), 
                 geom = "errorbar",  size = 0.4, width=.2,
                 position = position_dodge(0.6)) +
    
    theme_bw(base_size = 15) +
    
    ylab("") +
    
    xlab("") +
    
    facet_wrap(~variable,scales = "free_y",nrow = 2,ncol = 1) +
    
    scale_shape_manual(values=c(16, 1),name="Rearing\nsystem")
  
  ggsave("output/exp4_alpha_residue.jpeg", height = 10, width = 9, units = 'cm',dpi = "print")
  
  
  # larvae
  
  alpha_d_2_mean %>% 
    
    filter(type== "larvae" & experiment== "4") %>% 
    
    mutate(waste=case_when(waste=="tomato"~"TP",
                           waste=="foodwaste"~"DFW",
                           waste=="wine"~"WP",
                           TRUE~"")) %>% 
    
    ggplot(aes(waste,mean)) +
    
    geom_hpline(data=alpha_d_2_exp4_larvae,
                aes(waste, value, shape=system),
                position = position_dodge(0.6),
                width = 0.2, size = 0.7,
                stat = "summary") +
    
    geom_jitter(data=alpha_d_2_exp4_larvae, 
                aes(waste, value, shape = system), 
                position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
                size = 2.5) +
    
    stat_summary(data=alpha_d_2_exp4_larvae,
                 aes(waste,value, shape = system),
                 fun.data="mean_sdl",  fun.args = list(mult=1), 
                 geom = "errorbar",  size = 0.4, width=.2,
                 position = position_dodge(0.6)) +
    
    theme_bw(base_size = 15) +
    
    ylab("") +
    
    xlab("") +
    
    scale_shape_discrete(name="Rearing\nsystem") +
    
    facet_wrap(~variable,scales = "free_y",nrow = 2,ncol = 1) +
    
    scale_shape_manual(values=c(16, 1),name="Rearing\nsystem")
  
  ggsave("output/exp4_alpha_larvae.jpeg", height = 10, width = 9, units = 'cm',dpi = "print")
  
  
  #plot alpha diversity metrics of experiment 5 (WWP, TP with inoculant) -----
    
    # residue
    
    alpha_d_2_exp5_residue <-
      
      alpha_d_2 %>% 
      
      filter(type== "residue" & experiment== "5") %>% 
      
      filter(!type_detail == "treatment_foodwaste") %>% 
      
      unite("type_detail",c("waste","type_detail")) %>% 
      
      mutate(type_detail=case_when(type_detail=="wine_sterile inoculant"~"WP (+)",
                                   type_detail=="wine_inoculant"~"WP (o)",
                                   type_detail=="tomato_sterile inoculant"~"TP (+)",
                                   type_detail=="tomato_inoculant"~"TP (o)",
                                   TRUE~"")) 
    
    alpha_d_2_mean %>% 
      
      filter(type== "residue" & experiment== "5") %>% 
      
      filter(!type_detail == "treatment_foodwaste") %>% 
      
      unite("type_detail",c("waste","type_detail")) %>% 
      
      mutate(type_detail=case_when(type_detail=="wine_sterile inoculant"~"WP (+)",
                                   type_detail=="wine_inoculant"~"WP (o)",
                                   type_detail=="tomato_sterile inoculant"~"TP (+)",
                                   type_detail=="tomato_inoculant"~"TP (o)",
                                   TRUE~"")) %>% 
    
    
    ggplot(aes(type_detail,mean)) +
      
      geom_hpline(data=alpha_d_2_exp5_residue,
                  aes(type_detail, value, shape=system),
                  position = position_dodge(0.6),
                  width = 0.2, size = 0.7,
                  stat = "summary") +
      
      geom_jitter(data=alpha_d_2_exp5_residue, 
                  aes(type_detail, value, shape = system), 
                  position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.6),
                  size = 2.5
      ) +
      
      stat_summary(data=alpha_d_2_exp5_residue,
                   aes(type_detail, value, shape = system),
                   fun.data="mean_sdl",  fun.args = list(mult=1), 
                   geom = "errorbar",  size = 0.4, width=.2,
                   position = position_dodge(0.6)
      ) +
      
      theme_bw(base_size = 15) +
      
      ylab("") +
      
      xlab("") +
      
      facet_wrap(~variable,scales = "free_y",nrow = 2,ncol = 1) +
      
      geom_vline(xintercept = c(2.5), linetype="dashed") +

      scale_shape_manual(values=c(16, 21),name="Rearing\nsystem") 
    
    ggsave("output/exp5_alpha_residue.jpeg", height = 10, width = 12, units = 'cm',dpi = "print")
    
    
    # larvae
   
    alpha_d_2_exp5_larvae <-
      
      alpha_d_2 %>% 
      
      filter(type== "larvae" & experiment== "5") %>% 
      
      filter(!type_detail == "treatment_foodwaste") %>% 
      
      unite("type_detail",c("waste","type_detail")) %>% 
      
      mutate(type_detail=case_when(type_detail=="wine_sterile inoculant"~"WP (+)",
                                   type_detail=="wine_inoculant"~"WP (o)",
                                   type_detail=="tomato_sterile inoculant"~"TP (+)",
                                   type_detail=="tomato_inoculant"~"TP (o)",
                                   TRUE~"")) 
    
    alpha_d_2_mean %>% 
      
      filter(type== "larvae" & experiment== "5") %>% 
      
      filter(!type_detail == "treatment_foodwaste") %>% 
      
      unite("type_detail",c("waste","type_detail")) %>% 
      
      mutate(type_detail=case_when(type_detail=="wine_sterile inoculant"~"WP (+)",
                                   type_detail=="wine_inoculant"~"WP (o)",
                                   type_detail=="tomato_sterile inoculant"~"TP (+)",
                                   type_detail=="tomato_inoculant"~"TP (o)",
                                   TRUE~"")) %>% 
      
      
      ggplot(aes(type_detail,mean)) +
      
      geom_hpline(data=alpha_d_2_exp5_larvae,
                  aes(type_detail, value, shape=system),
                  position = position_dodge(0.6),
                  width = 0.2, size = 0.7,
                  stat = "summary") +
      
      geom_jitter(data=alpha_d_2_exp5_larvae, 
                  aes(type_detail, value, shape = system), 
                  position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.6),
                  size = 2.5
      ) +
      
      stat_summary(data=alpha_d_2_exp5_larvae,
                   aes(type_detail, value, shape = system),
                   fun.data="mean_sdl",  fun.args = list(mult=1), 
                   geom = "errorbar",  size = 0.4, width=.2,
                   position = position_dodge(0.6)
      ) +
      
      theme_bw(base_size = 15) +
      
      ylab("") +
      
      xlab("") +

      facet_wrap(~variable,scales = "free_y",nrow = 2,ncol = 1) +
      
      geom_vline(xintercept = c(2.5), linetype="dashed") +
    
      scale_shape_manual(values=c(16, 21),name="Rearing\nsystem") 
      
    
    ggsave("output/exp5_alpha_larvae.jpeg", height = 10, width = 12, units = 'cm',dpi = "print")
    
    
    
    
    
# Convert reads in relative abundance per sample ----------------------------------------------
  
  # overall
  
  d_norm <- transform_sample_counts(d, function(x) x / sum(x))
  
  d_norm <- prune_taxa(taxa_sums(d_norm) > 0, d_norm)

  
  experiments <- c("1","3","4","5")
  
  clusters <- c("1","2","3","4","5")
  

  # clustering -----
  
  mapfile <- read_tsv(mapfile)
  
  # experiment 4 ----
  
  # get sample names
  
  # larvae
  
  d_norm_2_exp4_mapfile_larvae <-   as_tibble(sample_data(d_norm_2_exp4_larvae)) %>%   
    
    unite("sample",c("waste","system"))
  
  # residue
  
  d_norm_2_exp4_mapfile_residue <-   as_tibble(sample_data(d_norm_2_exp4_residue)) %>%   
    
    unite("sample",c("waste","system"))
  
  
  # caclulate distance matrix (weighed unifrac distance considering abundance and richness)
  
  # larvae
  
  d_norm_2_exp4_larvae_unifrac = phyloseq::distance(d_norm_2_exp4_larvae, method="wunifrac")
  
  d_norm_2_exp4_larvae_matrix <- as.data.frame(otu_table(d_norm_2_exp4_larvae)) %>% t()
  
  # residue
  
  d_norm_2_exp4_residue_unifrac = phyloseq::distance(d_norm_2_exp4_residue, method="wunifrac")
  
  d_norm_2_exp4_residue_matrix <- as.data.frame(otu_table(d_norm_2_exp4_residue)) %>% t()
  
  
  # add sample names to distance matrix
  
  # larvae  
  
  dist_matrix_exp4_larvae <- dist_setNames(d_norm_2_exp4_larvae_unifrac, d_norm_2_exp4_mapfile_larvae$sample)
  
  d_norm_2_exp4_larvae_unifrac <- as.dist(dist_matrix_exp4_larvae)
  
  # residue
  
  dist_matrix_exp4_residue <- dist_setNames(d_norm_2_exp4_residue_unifrac, d_norm_2_exp4_mapfile_residue$sample)
  
  d_norm_2_exp4_residue_unifrac <- as.dist(dist_matrix_exp4_residue)
  
  # calculate clustering
  
  # larvae
  
  exp4_larvae_hclust     <- hclust(d_norm_2_exp4_larvae_unifrac, method="average")
  
  # residue
  
  exp4_residue_hclust     <- hclust( dist_matrix_exp4_residue, method="average")
  
  
  # determine clusters 
  
  # look at plot
  
  plot(exp4_larvae_hclust) # 4?
  
  plot(exp4_residue_hclust) # 4?
  
  # look at silhoutte plot
  
  # larvae 
  
  cut_exp4_larvae <- cutree(exp4_larvae_hclust,k=3)
  
  fviz_nbclust(d_norm_2_exp4_larvae_matrix, FUN = hcut, method = "silhouette")
  
  
  # residue
  
  cut_exp4_residue <- cutree(exp4_residue_hclust,k=2)
  plot(silhouette(cut_exp4_residue ,exp4_residue_hclust))
  
  fviz_nbclust(d_norm_2_exp4_residue_matrix, FUN = hcut, method = "silhouette")
  
  
  # calcilate prediction.strength 
  
  # larvae
  
  prediction.strength(d_norm_2_exp4_larvae_unifrac,
                      classification = "averagedist",
                      clustermethod = hclustCBI,method="average",
                      distances=TRUE,Gmin=3,Gmax=6,M=50)
  
  # residue
  
  prediction.strength(d_norm_2_exp4_residue_unifrac,
                      classification = "averagedist",
                      clustermethod = hclustCBI,method="average",
                      distances=TRUE,Gmin=3,Gmax=8,M=50)
  
  
  # calculate Jaccard similarity (>0.75)
  
  # larvae
  
  JS_exp4_larvae <- clusterboot(d_norm_2_exp4_larvae_unifrac,
                                bootmethod = "boot",
                                distances=TRUE,clustermethod = hclustCBI,
                                method="average",k=5)
  
  JS_exp4_larvae
  
  JS_exp4_larvae$result
  
  # residue
  
  JS_exp4_residue <- clusterboot(d_norm_2_exp4_residue_unifrac,
                                 bootmethod = "boot",
                                 distances=TRUE,
                                 clustermethod = hclustCBI,
                                 method="average",k=8)
  
  JS_exp4_residue
  
  JS_exp4_residue$result 
  
  
  # plot final dendogram 
  
  # larvae
  
  fviz_dend(exp4_larvae_hclust, cex = 0.5, k = 4, color_labels_by_k = TRUE,
            k_colors = c("#999999", "#E69F00", "#56B4E9", "#009E73"),main="")
  
  ggsave("output/dendo_exp4_larvae.jpeg", height = 10, width = 14, units = 'cm',dpi = 1200)
  
  
  
  exp4_larvae_hclust_cut <- cutree(exp4_larvae_hclust, k = 4)
  
  
  
  # residue
  
  fviz_dend(exp4_residue_hclust, cex = 0.5, k = 5, color_labels_by_k = TRUE,
            k_colors= c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442"), main="")
  
  
  ggsave("output/dendo_exp4_residue.jpeg", height = 10, width = 14, units = 'cm',dpi = 1200)
  
  
  # plot PCoA  
  
  fviz_cluster(hc.cut, ellipse.type = "convex")
  
  ?hcut
  
  # Use hcut() which compute hclust and cut the tree
  larvae_exp4_cut <- hcut(d_norm_2_exp4_larvae_matrix, k = 4, 
                          hc_method = "average",
                          hc_metric = "euclidian")
  
  # Visualize dendrogram
  fviz_dend(larvae_exp4_cut, show_labels = FALSE, rect = TRUE)
  # Visualize cluster
  
  fviz_cluster(list(data = d_norm_2_exp4_larvae_matrix, cluster = exp4_larvae_hclust_cut),
               palette = "jco", 
               ellipse.type = "convex", # Concentration ellipse
               repel = TRUE, # Avoid label overplotting (slow)
               show.clust.cent = FALSE, ggtheme = theme_minimal())
  
  ?fviz_cluster
  
  
  # experiment 5 ----
  
  
  # get sample names
  
  # larvae
  
  d_norm_2_exp5_mapfile_larvae <-   as_tibble(sample_data(d_norm_2_exp5_larvae)) %>%   
    
    unite("sample",c("waste","system","type_detail"))
  
  # residue
  
  d_norm_2_exp5_mapfile_residue <-   as_tibble(sample_data(d_norm_2_exp5_residue)) %>%   
    
    unite("sample",c("waste","system","type_detail"))
  
  
  # caclulate distance matrix (weighed unifrac distance considering abundance and richness)
  
  # larvae
  
  d_norm_2_exp5_larvae_unifrac = phyloseq::distance(d_norm_2_exp5_larvae, method="wunifrac")
  
  d_norm_2_exp5_larvae_matrix <- as.data.frame(otu_table(d_norm_2_exp5_larvae)) %>% t()
  
  # residue
  
  d_norm_2_exp5_residue_unifrac = phyloseq::distance(d_norm_2_exp5_residue, method="wunifrac")
  
  d_norm_2_exp5_residue_matrix <- as.data.frame(otu_table(d_norm_2_exp5_residue)) %>% t()
  
  
  # add sample names to distance matrix
  
  # larvae  
  
  dist_matrix_exp5_larvae <- dist_setNames(d_norm_2_exp5_larvae_unifrac, d_norm_2_exp5_mapfile_larvae$sample)
  
  d_norm_2_exp5_larvae_unifrac <- as.dist(dist_matrix_exp5_larvae)
  
  # residue
  
  dist_matrix_exp5_residue <- dist_setNames(d_norm_2_exp5_residue_unifrac, d_norm_2_exp5_mapfile_residue$sample)
  
  d_norm_2_exp5_residue_unifrac <- as.dist(dist_matrix_exp5_residue)
  
  # calculate clustering
  
  # larvae
  
  exp5_larvae_hclust     <- hclust(d_norm_2_exp5_larvae_unifrac, method="average")
  
  # residue
  
  exp5_residue_hclust     <- hclust(d_norm_2_exp5_residue_unifrac, method="average")
  
  
  # determine clusters 
  
  # look at plot
  
  plot(exp5_larvae_hclust) # 2?
  
  plot(exp5_residue_hclust) # 2?
  
  # look at silhoutte plot
  
  # larvae 
  
  cut_exp5_larvae <- cutree(exp5_larvae_hclust,k=2)
  
  fviz_nbclust(d_norm_2_exp5_larvae_matrix, FUN = hcut, method = "silhouette")
  
  
  ?fviz_nbclust
  
  # residue
  
  cut_exp5_residue <- cutree(exp5_residue_hclust,k=2)
  plot(silhouette(cut_exp5_residue ,exp5_residue_hclust))
  
  fviz_nbclust(d_norm_2_exp5_residue_matrix, FUN = hcut, method = "silhouette")
  
  
  # calcilate prediction.strength 
  
  # larvae
  
  prediction.strength(d_norm_2_exp5_larvae_unifrac,
                      classification = "averagedist",
                      clustermethod = hclustCBI,method="average",
                      distances=TRUE,Gmin=2,Gmax=3,M=50)
  
  # residue
  
  prediction.strength(d_norm_2_exp5_residue_unifrac,
                      classification = "averagedist",
                      clustermethod = hclustCBI,method="average",
                      distances=TRUE,Gmin=2,Gmax=3,M=50)
  
  
  # calculate Jaccard similarity (>0.75)
  
  # larvae
  
  JS_exp5_larvae <- clusterboot(d_norm_2_exp5_larvae_unifrac,
                                bootmethod = "boot",
                                distances=TRUE,clustermethod = hclustCBI,
                                method="average",k=2)
  
  JS_exp5_larvae
  
  JS_exp5_larvae$result
  
  # residue
  
  JS_exp5_residue <- clusterboot(d_norm_2_exp5_residue_unifrac,
                                 bootmethod = "boot",
                                 distances=TRUE,
                                 clustermethod = hclustCBI,
                                 method="average",k=3)
  
  JS_exp5_residue
  
  JS_exp5_residue$result 
  
  
  # plot final dendogram 
  
  # larvae
  
  fviz_dend(exp5_larvae_hclust, cex = 0.5, k = 3, color_labels_by_k = TRUE,
            k_colors = c("#999999", "#E69F00", "#56B4E9", "#009E73"),main="")
  
  ggsave("output/dendo_exp5_larvae.jpeg", height = 10, width = 14, units = 'cm',dpi = 1200)
  
  
  
  exp5_larvae_hclust_cut <- cutree(exp5_larvae_hclust, k = 4)
  
  
  
  # residue
  
  fviz_dend(exp5_residue_hclust, cex = 0.5, k = 3, color_labels_by_k = TRUE,
            k_colors= c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442"), main="")
  
  
  ggsave("output/dendo_exp5_residue.jpeg", height = 10, width = 14, units = 'cm',dpi = 1200)
  
  ?fviz_dend
  
  # plot PCoA  
  
  fviz_cluster(hc.cut, ellipse.type = "convex")
  
  ?hcut
  
  # Use hcut() which compute hclust and cut the tree
  larvae_exp4_cut <- hcut(d_norm_2_exp4_larvae_matrix, k = 4, 
                          hc_method = "average",
                          hc_metric = "euclidian")
  
  # Visualize dendrogram
  fviz_dend(larvae_exp4_cut, show_labels = FALSE, rect = TRUE)
  # Visualize cluster
  
  fviz_cluster(list(data = d_norm_2_exp4_larvae_matrix, cluster = exp4_larvae_hclust_cut),
               palette = "jco", 
               ellipse.type = "convex", # Concentration ellipse
               repel = TRUE, # Avoid label overplotting (slow)
               show.clust.cent = FALSE, ggtheme = theme_minimal())
  
  ?fviz_cluster
  
  
  # 
  # result <- pvclust(d_norm_2_exp4_larvae, 
  #                   method.dist = function(x){
  #                     vegdist(d_norm_2_exp4_larvae_matrix, "bray")
  #                   },
  #                   method.hclust = "ward.D",
  #                   nboot=1000)
  # 
  # result
  # 
  # 
  # # Create the dendrogram with p values
  # plot(result, hang=-1, main="Alberta Natural Subregions")
  # # add rectangles around groups highly supported by the data
  # pvrect(result, alpha=.95)
  
  # permanova
  
  # waste = tomato
  
  # https://www.nicholas-ollberding.com/post/introduction-to-the-statistical-analysis-of-microbiome-data-in-r/
  # https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/adonis
  # https://mattsigal.github.io/eqcov_supp/betadisp-ex.html
  # researchgate.net/post/How_can_I_do_PerMANOVA_pairwise_contrasts_in_R
  
  d_norm_2_exp4_tomato <-  subset_samples(d_norm_2_exp4, waste=="tomato" & type == "residue")
  
  d_norm_2_exp4_tomato <- prune_taxa(taxa_sums(d_norm_2_exp4_tomato) > 0, d_norm_2_exp4_tomato)                        
  
  # calculate distance matrix
  
  d_norm_2_exp4_tomato_unifrac = phyloseq::distance(d_norm_2_exp4_tomato, method="unifrac", weighted=F)
  
  # homogeniety of variance, Dispersion test and plot, adonis is different within-group variation (dispersion) 
  
  # What to do when there are significantly different dispersion between groups
  # https://www.researchgate.net/post/PERMANOVA_with_non-homogeneous_dispersion_among_test_groups_betadisper
  
  
  # d_norm_2_exp4_tomato_unifrac_dispr <- vegan::betadisper(d_norm_2_exp4_tomato_unifrac, phyloseq::sample_data(d_norm_2_exp4_tomato)$system)
  # 
  # d_norm_2_exp4_tomato_unifrac_dispr
  # 
  # plot(d_norm_2_exp4_tomato_unifrac_dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
  # 
  # boxplot(d_norm_2_exp4_tomato_unifrac_dispr, main = "", xlab = "")
  # 
  # permutest(d_norm_2_exp4_tomato_unifrac_dispr)
  
  # adonis 
  
  vegan::adonis(d_norm_2_exp4_tomato_unifrac ~ phyloseq::sample_data(d_norm_2_exp4_tomato)$system,permutations = 100)
  
  
  
# beta diversity -----    
  

# experiment 4: WWP, TP and DFW without inoculants in open/closed rearing system ----
    
    # subset samples of experiment 4
    
    d_norm_2_exp4 <-  subset_samples(d_norm, 
                                     
                                     (run=="2a" | run == "2b") & 
                                       
                                     experiment =="4" & 
                                       
                                     !waste=="5DOL")
    
    # subset samples - separate for larvae and residue samples
    
    d_norm_2_exp4_larvae <-  subset_samples(d_norm_2_exp4, type == "larvae")
    d_norm_2_exp4_residue <-  subset_samples(d_norm_2_exp4, type == "residue")
    
    # convert clusters into factors
    
    sample_data(d_norm_2_exp4_larvae)$Cluster <- factor(x=sample_data(d_norm_2_exp4_larvae)$Cluster, levels=clusters)
    sample_data(d_norm_2_exp4_residue)$Cluster <- factor(x=sample_data(d_norm_2_exp4_residue)$Cluster, levels=clusters)
    
    # remove taxa with no reads
    
    d_norm_2_exp4_larvae <- prune_taxa(taxa_sums(d_norm_2_exp4_larvae) > 0, d_norm_2_exp4_larvae)                        
    d_norm_2_exp4_residue <- prune_taxa(taxa_sums(d_norm_2_exp4_residue) > 0, d_norm_2_exp4_residue)   
    
    # PCoA
    
    d_norm_2_exp4_larvae_ordinate_pcoa <- ordinate(d_norm_2_exp4_larvae, "PCoA", "unifrac", weighted=TRUE)
    d_norm_2_exp4_residue_ordinate_pcoa <- ordinate(d_norm_2_exp4_residue, "PCoA", "unifrac", weighted=TRUE)

    
    # plot PCoA

      # larvae samples
      
      exp4_pcoa_larvae <- plot_ordination(d_norm_2_exp4_larvae, d_norm_2_exp4_larvae_ordinate_pcoa,
                                          shape="system", title="") + 
        
        geom_point(size=5) +
        
        geom_text(mapping = aes(label = Day), size = 3, vjust = 1.5) +
        
        theme_bw(base_size = 15) +
        
        geom_encircle(aes(group=Cluster,fill=Cluster),alpha=.3) +
        
        annotate(geom="text", x=-.12, y=.18, label="DFW",
                 color="black") +
      
        annotate(geom="text", x=-.12, y=-.19, label="TP",
                 color="black") +
        
        annotate(geom="text", x=.3, y=.09, label="WWP",
                 color="black") +
        
        scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73")) +
        
        scale_shape_manual(values=c(16, 1),name="Rearing\nsystem")
      
        # remove first layer of points (otherwise double open symbol)
      
        exp4_pcoa_larvae$layers <- exp4_pcoa_larvae$layers[-1]
  
        exp4_pcoa_larvae
      
        # save
      
        ggsave("output/PCoA_exp4_larvae.jpeg", height = 9, width = 11, units = 'cm',dpi = "print")
      
     
    # residue samples   
    
      exp4_pcoa_residue <-   
        
      plot_ordination(d_norm_2_exp4_residue, d_norm_2_exp4_residue_ordinate_pcoa,shape="system", 
                      title="") + 
        
        geom_point(size=5) +
        
        geom_text(mapping = aes(label = Day), size = 3, vjust = 1.5) +
        
        theme_bw(base_size = 15) +
        
        geom_encircle(aes(group=Cluster,fill=Cluster),alpha=.3) +
        
        annotate(geom="text", x=0.03, y=.25, label="TP",
                 color="black") +
        
        annotate(geom="text", x=-.20, y=-.16, label="WWP",
                 color="black") +
        
        annotate(geom="text", x=.3, y=-.09, label="DFW",
                 color="black") +
      
        scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442")) +
        
        scale_shape_manual(values=c(16, 1),name="Rearing\nsystem")
      
      
      # remove first layer of points (otherwise double open symbol)
      
      exp4_pcoa_residue$layers <- exp4_pcoa_residue$layers[-1]
      
      exp4_pcoa_residue
      
      ggsave("output/PCoA_exp4_residue.jpeg", height = 9, width = 11, units = 'cm',dpi = 1200)
   
      
# experiment 5: Influence of inoculant on community --------
  
    # define symbols indicating type of treatment (o = autoclaved/sterile inoculant, + = inoculant)  
    
    treatment_detail = c("o","+")
      
    # subset samples of experiment 5 
      
    d_norm_2_exp5 <-  subset_samples(d_norm,  (run=="2a" | run == "2b") & 
                                     
                                            experiment =="5" & 
                                     
                                            (type=="residue"|type=="larvae") &
                                     
                                            !type_detail=="treatment_foodwaste")
    
    
    # subset samples - separate for larvae and residue samples
      
    d_norm_2_exp5_larvae <-  subset_samples(d_norm_2_exp5, type == "larvae")
    d_norm_2_exp5_residue <-  subset_samples(d_norm_2_exp5, type == "residue")
      
    # convert clusters into factors   
      
    sample_data(d_norm_2_exp5_larvae)$Cluster <- factor(x=sample_data(d_norm_2_exp5_larvae)$Cluster, levels=clusters)
    sample_data(d_norm_2_exp5_residue)$Cluster <- factor(x=sample_data(d_norm_2_exp5_residue)$Cluster, levels=clusters)
    
    # remove taxa with no reads
      
    d_norm_2_exp5_larvae <- prune_taxa(taxa_sums(d_norm_2_exp5_larvae) > 0, d_norm_2_exp5_larvae)                        
    d_norm_2_exp5_residue <- prune_taxa(taxa_sums(d_norm_2_exp5_residue) > 0, d_norm_2_exp5_residue)
    
    # PCoA
    
    d_norm_2_exp5_larvae_ordinate_pcoa <- ordinate(d_norm_2_exp5_larvae, "PCoA", "unifrac", weighted=TRUE) 
    d_norm_2_exp5_residue_ordinate_pcoa <- ordinate(d_norm_2_exp5_residue, "PCoA", "unifrac", weighted=TRUE) 
    
    
    # plot PCoA
    
    # larvae samples
    
      exp5_pcoa_larvae <- 
    
      plot_ordination(d_norm_2_exp5_larvae, d_norm_2_exp5_larvae_ordinate_pcoa,shape="system",
                    title="") + 
        
      geom_point(size=5) +

      theme_bw(base_size = 15) +
        
      geom_text(aes(label=treatment_detail),nudge_y=-0.025,nudge_x = .025,check_overlap = TRUE) +

      geom_encircle(aes(group=Cluster,fill=Cluster),alpha=.3) +

      annotate(geom="text", x=-.25, y=-0.14, label="TP",
               color="black") +

      annotate(geom="text", x=-.06, y=.12, label="WWP",
               color="black") +

      scale_fill_manual(values=c("#999999", "#E69F00")) +
        
      scale_shape_manual(values=c(16, 1),name="Rearing\nsystem")
      
      # remove first layer of points (otherwise double open symbol)
      
      exp5_pcoa_larvae$layers <- exp5_pcoa_larvae$layers[-1]
      
      exp5_pcoa_larvae
    
    ggsave("output/PCoA_exp5_larvae.jpeg", height = 9, width = 11, units = 'cm',dpi = "print")
    
  
    # residue samples
    
    exp5_pcoa_residue <- 
    
    plot_ordination(d_norm_2_exp5_residue, d_norm_2_exp5_residue_ordinate_pcoa,shape="system",
                    title="") + 
      
      geom_point(size=5) +
      
      theme_bw(base_size = 15) +
      
      geom_encircle(aes(group=Cluster,fill=Cluster),alpha=.3) +
      
      geom_text(aes(label=treatment_detail),nudge_y=-0.025,nudge_x = .025,check_overlap = TRUE) +
      
      annotate(geom="text", x=.23, y=.01, label="TP",
               color="black") +

      annotate(geom="text", x=.02, y=-.14, label="WWP",
               color="black") +
      
      scale_fill_manual(values=c("#999999","#E69F00")) +
      
      scale_shape_manual(values=c(16,1),name="Rearing\nsystem")
    
    # remove first layer of points (otherwise double open symbol)
    
    exp5_pcoa_residue$layers <- exp5_pcoa_residue$layers[-1]
    
    exp5_pcoa_residue
    
    ggsave("output/PCoA_exp5_residue.jpeg", height = 9, width = 11, units = 'cm',dpi = "print")
    