# Libraries and functions -------------------------------------------------

## clean/reset environment 

# rm(list=ls()) 



## R and Bioconductor libraries 


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
library(metacoder)
library(ampvis2)
library(openxlsx)
library(UpSetR)
library(ComplexHeatmap)
library(ggplot2)
library(grid)
library(cowplot)
library(plyr)

## Functions
source("code/Rfunctions.R")




# Read phyloseq object from pyloseq.R

d

# after removal of samples


# Create ampvis file ------------------------------------------------------


#Help: (https://madsalbertsen.github.io/ampvis2/articles/faq.html)


# run 2 

otu_table_d_df <- 
  
  data.frame(ZOTU = rownames(phyloseq::otu_table(d)@.Data),
             phyloseq::otu_table(d)@.Data,
             phyloseq::tax_table(d)@.Data,
             check.names = FALSE)

metadata_d_df <- data.frame(phyloseq::sample_data(d),
                              check.names = FALSE)


amp_d <- amp_load(otu_table_d_df, metadata_d_df)



# Heat maps ----------------------------------------------------------------


# overall

  # waste

  amp_wastes <-amp_subset_samples(amp_d, type == "waste")

  #genus level
  
  amp_heatmap(amp_wastes,
              group_by = c("replicate"),
              facet_by = "waste",
              tax_aggregate = "Genus",
              tax_show = 20,
              plot_na = TRUE,
              normalise= TRUE) +
    
    theme(axis.text.x = element_text(angle = 45, size=12, vjust = 1),
          axis.text.y = element_text(size=12,face = "italic"),
          legend.position="right") +
    
    labs(fill = "Relative\nabundance")

  #phylum level
  
  amp_heatmap(amp_wastes,
              group_by = c("replicate"),
              facet_by = "waste",
              tax_aggregate = "Phylum",
              tax_show = 4,
              plot_na = TRUE,
              normalise= TRUE) +
    
    theme(axis.text.x = element_text(angle = 45, size=12, vjust = 1),
          axis.text.y = element_text(size=12,face = "italic"),
          legend.position="right") +
    
    labs(fill = "Relative\nabundance")


# run 2 ------------


# wastes ----


  amp_wastes_d_2 <-amp_subset_samples(amp_d, 
                                    
                                    ( run== "2a" | run == "2b") & 
                                      
                                    ( type == "waste") &
                                      
                                    (  waste == "foodwaste"| waste == "tomato")
                                    
                                    )
  #phylum level
  
  amp_heatmap(amp_wastes_d_2,
              group_by = c("replicate"),
              facet_by = "waste",
              tax_aggregate = "Phylum",
              tax_show = 4,
              plot_na = TRUE,
              normalise= TRUE) +
    
    theme(axis.text.x = element_text(angle = 45, size=12, vjust = 1),
          axis.text.y = element_text(size=12,face = "italic"),
          legend.position="right") +
    
    labs(fill = "Relative\nabundance")
  
  ggsave("output/waste_heatmap_phylum.jpeg", height = 10, width = 16, units = 'cm',dpi = "print")
  
  #genus level

  amp_wastes_d_2$metadata$waste   <- mapvalues(amp_wastes_d_2$metadata$waste, from = c("wine", "tomato","foodwaste"), to = c("WWP", "TP","DFW"))
  
  
    amp_heatmap(amp_wastes_d_2,
              group_by = c("replicate"),
              facet_by = "waste",
              tax_aggregate = "Genus",
              tax_show = 10,
              plot_na = TRUE,
              normalise= TRUE) +
    
    theme(axis.text.x = element_text(angle = 45, size=12, vjust = 1),
          axis.text.y = element_text(size=12,face = "italic"),
          legend.position="right") +
    
    labs(fill = "Relative\nabundance")

     ggsave("output/waste_heatmap_genus.jpeg", height = 10, width = 12, units = 'cm',dpi = "print")
    

     #Class level
     
     amp_heatmap(amp_wastes_d_2,
                 group_by = c("replicate"),
                 facet_by = "waste",
                 tax_aggregate = "Class",
                 tax_show = 10,
                 plot_na = TRUE,
                 normalise= TRUE) +
       
       theme(axis.text.x = element_text(angle = 45, size=12, vjust = 1),
             axis.text.y = element_text(size=12,face = "italic"),
             legend.position="right") +
       
       labs(fill = "Relative\nabundance")
     
     ggsave("output/waste_heatmap_class.jpeg", height = 10, width = 12, units = 'cm',dpi = "print")
     
       
     
# experiment 4 ------     
     
     
 amp_exp4 <-amp_subset_samples(amp_d, experiment=="4" & !type=="5DOL")
 amp_exp4_residue <-amp_subset_samples(amp_d, experiment=="4" & type=="residue")
 amp_exp4_larvae <-amp_subset_samples(amp_d, experiment=="4" &  type=="larvae")
    
 
 amp_exp4$metadata$waste   <- mapvalues(amp_exp4$metadata$waste, from = c("wine", "tomato","foodwaste"), to = c("White wine\npomace", "Tomato\npomace","Digested\nfood waste"))

 
  # differences among systems for digested food waste
 
 amp_exp4_larvae <-amp_subset_samples(amp_d, experiment=="4" &  type=="larvae" & waste =="foodwaste")
 
 
 amp_exp4_larvae_system <- 
 
 amp_heatmap(amp_exp4_larvae,
             group_by = c("replicate"),
             facet_by = "system",
             tax_aggregate = "Family",
             tax_show = 10,
             plot_na = TRUE,
             normalise= TRUE) +
   
   theme(axis.text.x = element_text(angle = 45, size=12, vjust = 1),
         axis.text.y = element_text(size=12,face = "italic"),
         legend.position="right") +
   
   labs(fill = "Relative\nabundance")
 
 amp_exp4_larvae_system$data 
 
  #phylum level
     

 
     amp_heatmap(amp_exp4,
                 group_by = c("name"),
                 facet_by = "type",
                 tax_aggregate = "Phylum",
                 tax_show = 4,
                 plot_na = TRUE,
                 normalise= TRUE) +
       
       theme(axis.text.x = element_text(angle = 45, size=12, vjust = 1),
             axis.text.y = element_text(size=12,face = "italic"),
             legend.position="right") +
       
       labs(fill = "Relative\nabundance")
     
     ggsave("output/exp4_heatmap_phylum.jpeg", height = 10, width = 14, units = 'cm',dpi = "print")
     
    
     
     amp_exp4 <-amp_subset_samples(amp_d, experiment=="4" & !type=="5DOL")
     amp_exp4_residue <-amp_subset_samples(amp_d, experiment=="4" & type=="residue")
     amp_exp4_larvae <-amp_subset_samples(amp_d, experiment=="4" &  type=="larvae")
     
    #Family level
     
     
 
     amp_heatmap(amp_exp4,
                 group_by = c("name"),
                 facet_by = "type",
                 tax_aggregate = "Family",
                 tax_show = 10,
                 plot_na = TRUE,
                 normalise= TRUE) +
       
       theme(axis.text.x = element_text(angle = 45, size=12, vjust = 1),
             axis.text.y = element_text(size=12,face = "italic"),
             legend.position="right") +
       
       labs(fill = "Relative\nabundance") 
      
     ggsave("output/exp4_heatmap_family.jpeg", height = 10, width = 16, units = 'cm',dpi = "print")
     
     
     
     #genus level
     
      # overall
     
     amp_heatmap(amp_exp4,
                 group_by = c("name"),
                 facet_by = "type",
                 tax_aggregate = "Genus",
                 tax_show = 15,
                 plot_na = TRUE,
                 normalise= TRUE) +
       
       theme(axis.text.x = element_text(angle = 45, size=12, vjust = 1),
             axis.text.y = element_text(size=12,face = "italic"),
             legend.position="right") +
       
       labs(fill = "Relative\nabundance")
     
     ggsave("output/exp4_heatmap_genus.jpeg", height = 16, width = 16, units = 'cm',dpi = "print")
    
      # residue
     
     
     amp_exp4_residue$metadata$waste   <- mapvalues(amp_exp4_residue$metadata$waste, from = c("wine", "tomato","foodwaste"), to = c("White wine\npomace", "Tomato\npomace","Digested\nfood waste"))
     amp_exp4_larvae$metadata$waste   <- mapvalues(amp_exp4_larvae$metadata$waste, from = c("wine", "tomato","foodwaste"), to = c("White wine\npomace", "Tomato\npomace","Digested\nfood waste"))
     

      
      amp_heatmap(amp_exp4_residue,
                 group_by = c("name"),
                 tax_aggregate = "Genus",
                 tax_show = 15,
                 plot_na = TRUE,
                 normalise= TRUE) +
       
       theme(axis.text.x = element_text(angle = 45, size=12, vjust = 1),
             axis.text.y = element_text(size=12,face = "italic"),
             legend.position="right") +
       
       labs(fill = "Relative\nabundance")
     
     ggsave("output/exp4_heatmap_genus_residue_defense.jpeg", height = 14, width = 12, units = 'cm',dpi = "print")
     

     
      # larvae
     
     amp_heatmap(amp_exp4_larvae,
                 group_by = c("name"),
                 tax_aggregate = "Genus",
                 tax_show = 15,
                 plot_na = TRUE,
                 normalise= TRUE) +
       
       theme(axis.text.x = element_text(angle = 45, size=12, vjust = 1),
             axis.text.y = element_text(size=12,face = "italic"),
             legend.position="right") +
       
       labs(fill = "Relative\nabundance")
     
     ggsave("output/exp4_heatmap_genus_larvae.jpeg", height = 14, width = 14, units = 'cm',dpi = "print")
     
     
      #class level
     
     amp_heatmap(amp_exp4,
                 group_by = c("waste"),
                 facet_by = "type",
                 tax_aggregate = "Family",
                 tax_show = 10,
                 plot_na = TRUE,
                 normalise= TRUE) +
       
       theme(axis.text.x = element_text(angle = 45, size=12, vjust = 1),
             axis.text.y = element_text(size=12,face = "italic"),
             legend.position="right") +
       
       labs(fill = "Relative\nabundance")
     
     ggsave("output/exp4_heatmap_family.jpeg", height = 10, width = 16, units = 'cm',dpi = "print")
     
     
       
# residue vs. inoculant (exp 1 and 3) ----

  amp_inoc_1 <-amp_subset_samples(amp_d, 
                                
                                (run == "2a" | run== "2b") &
                                
                                (others == "inoculant" | others == "residue") & 
                                  
                                (experiment=="1" | experiment=="3") & !waste=="fly")


  #genus level

  amp_heatmap(amp_inoc_1,
            group_by = c("others"),
            facet_by = "waste",
            tax_aggregate = "Genus",
            tax_show = 15,
            plot_na = TRUE,
            normalise= TRUE) +
  
  theme(axis.text.x = element_text(angle = 45, size=12, vjust = 1),
        axis.text.y = element_text(size=12,face = "italic"),
        legend.position="right") +
  
  labs(fill = "Relative\nabundance")
  
  ggsave("output/inoc1_heatmap_tomato_genus.jpeg", height = 15, width = 23, units = 'cm',dpi = "print")



 # residue vs. inoculant (exp 4 and 5) ----

  
  #tomato (closed) ----
  
  amp_inoc_2_tomato <-amp_subset_samples(amp_d,
                                  
                                  (run== "2a" | run== "2b") &
                                    
                                    waste == "tomato" &
                                    
                                    ((others == "inoculant" | others == "residue") & system == "closed") & 
                                    
                                    (experiment == "4" | experiment == "5") 
                                    )


  # genus level    
  
  amp_inoc_2_tomato$metadata$type   <- factor(amp_inoc_2_tomato$metadata$type  , levels = c("residue","inoculant"))
                 
                                    
  amp_heatmap_inoc_tomato <-
                                                     
  amp_heatmap(amp_inoc_2_tomato,
              group_by = c("type"),
              facet_by = "system",
              tax_aggregate = "Genus",
              tax_show = 50,
              plot_na = TRUE,
              normalise= TRUE) +
    
    theme(axis.text.x = element_text(angle = 45, size=12, vjust = 1),
          axis.text.y = element_text(size=12,face = "italic"),
          legend.position="right") +
    
    labs(fill = "Relative\nabundance")
  
    ggsave("output/inoc2_heatmap_tomato_genus.jpeg", height = 10, width = 12, units = 'cm',dpi = "print")
  
  
  # wine (open and closed) ----
  
  amp_inoc_2_wine <-amp_subset_samples(amp_d,
                                  
                                  (run== "2a" | run== "2b") &
                                    
                                    waste == "wine" &
                                    
                                    ((others == "inoculant" | others == "residue")) & 
                                    
                                    (experiment == "4" | experiment == "5") 
                                  )
    amp_inoc_2_wine$metadata$type   <- factor(amp_inoc_2_wine$metadata$type  , levels = c("residue","inoculant"))
    
    
  # genus level    
  
  amp_heatmap_inoc_wine <- 
    
  amp_heatmap(amp_inoc_2_wine,
              group_by = c("type"),
              facet_by = "system",
              tax_aggregate = "Genus",
              tax_show = 15,
              plot_na = TRUE,
              normalise= TRUE) +
    
    theme(axis.text.x = element_text(angle = 45, size=12, vjust = 1),
          axis.text.y = element_text(size=12,face = "italic"),
          legend.position="right") +
    
    labs(fill = "Relative\nabundance")
  
  ggsave("output/inoc2_heatmap_wine_genus.jpeg", height = 10, width = 13, units = 'cm',dpi = "print")
  
  # Do we find typical gut microbes
  
  amp_heatmap_inoc_tomato$data %>% 
    
    filter(Display=="Providencia")
  
  
  
  # food waste (open)
  
  amp_inoc_2 <-amp_subset_samples(amp_d,
                                  
                                  (run== "2a" | run== "2b") &
                                    
                                    waste == "foodwaste" &
                                    
                                    ((others == "inoculant" | others == "residue") & system == "open") & 
                                    
                                    (experiment == "4" | experiment == "5") 
                                    )
  
  # genus level    
  
  amp_heatmap(amp_inoc_2,
              group_by = c("others"),
              tax_aggregate = "Genus",
              facet_by = "waste",
              tax_show = 15,
              plot_na = TRUE,
              normalise= TRUE) +
    
    theme(axis.text.x = element_text(angle = 45, size=12, vjust = 1),
          axis.text.y = element_text(size=12,face = "italic"),
          legend.position="right") +
    
    labs(fill = "Relative\nabundance")
  
  ggsave("output/inoc2_heatmap_foodwaste_genus.jpeg", height = 15, width = 23, units = 'cm',dpi = "print")
  
  

# Venn diagramms -----------------------------------------------------------



#Experiment 1 & 3 inoculant vs. residue -----

  # tomato ----

  amp_inoc_1 <-amp_subset_samples(amp_d,
                                  
                                  (run == "2a" | run== "2b") &
                                  
                                  waste=="tomato" & 
                                    
                                  (others == "inoculant" | others == "residue") & 
                                    
                                  (experiment=="1" | experiment=="3"))
  
  amp_inoc_1_venn <- amp_venn(amp_inoc_1,group_by = "others",cut_a = 0.01,cut_f = 80,detailed_output = TRUE)
  
  
  
  amp_inoc_1_venn$Otutable %>% filter(Shared=="Core")

  ggsave("output/venn_residue_inoc_1_tomato.jpeg", height = 16, width = 17, units = 'cm',dpi = 1200)

  
  # wine ----
  
  amp_inoc_1 <-amp_subset_samples(amp_d,
                                  
                                  (experiment == "1a" | experiment == "1b") |
                                    
                                    waste=="wine" & 
                                    
                                    (others == "inoculant" | others == "residue") & 
                                    
                                    (experiment=="1" | experiment=="3"))
  
  amp_inoc_1_venn <- amp_venn(amp_inoc_1,group_by = "others",cut_a = 0.01,cut_f = 80,detailed_output = TRUE)
  
  
  
  amp_inoc_1_venn$Otutable %>% filter(Shared=="Core") 
  
  ggsave("output/venn_residue_inoc1_wine.jpeg", height = 16, width = 17, units = 'cm',dpi = 1200)
  

  
#Experiment 4 & 5 inoculant vs. residue -----

  # tomato ----
  
  amp_inoc_2 <-amp_subset_samples(amp_d,
                                  
                                  (run== "2a" | run== "2b") &
                                    
                                  waste == "tomato" &
                                    
                                  ((others == "inoculant" | others == "residue") & system == "closed") & 
                                    
                                  (experiment == "4" | experiment == "5") 
                                    
                                  )
  
  amp_inoc_2_venn <- amp_venn(amp_inoc_2,group_by = "others",cut_a = 0.01,cut_f = 80,detailed_output = TRUE)
  
  
  
  amp_inoc_2_venn$Otutable %>% filter(Shared=="Core") %>% View()
  
  ggsave("output/venn_residue_inoc2_tomato.jpeg", height = 16, width = 17, units = 'cm',dpi = 1200)
  
  
  # wine
  
  amp_inoc_2 <-amp_subset_samples(amp_d,
                                  
                                  (run== "2a" | run== "2b") &
                                    
                                    waste == "wine" &
                                    
                                    ((others == "inoculant" | others == "residue") & system == "open") & 
                                    
                                    (experiment == "4" | experiment == "5") 
                                  
                                    )
  
  amp_inoc_2_venn <- amp_venn(amp_inoc_2,group_by = "others",cut_a = 0.1,cut_f = 80,detailed_output = TRUE)
  
  
  
  amp_inoc_2_venn$Otutable %>% filter(Shared=="Core") %>% View()
  
  ggsave("output/venn_residue_inoc2_wine_open.jpeg", height = 16, width = 17, units = 'cm',dpi = 1200)
  
  
  # food waste
  
  amp_inoc_2 <-amp_subset_samples(amp_d,
                                  
                                  (run== "2a" | run== "2b") &
                                    
                                    waste == "foodwaste" &
                                    
                                    ((others == "inoculant" | others == "residue") & system == "open") & 
                                    
                                    (experiment == "4" | experiment == "5") 
                                  
                                    )
  
  amp_inoc_2_venn <- amp_venn(amp_inoc_2,group_by = "others",cut_a = 0.01,cut_f = 80,detailed_output = TRUE)
  
  
  
  amp_inoc_2_venn$Otutable %>% filter(Shared=="Core") %>% View()
  
  ggsave("output/venn_residue_inoc2_foodwaste_open.jpeg", height = 16, width = 17, units = 'cm',dpi = 1200)
  
  
  
#Experiment 4: waste vs. residue ---- 

  #  tomato (change system)

  amp_waste_residue <- amp_subset_samples(amp_d, 
                                          
                                          (experiment == "1a" | experiment == "1b") |
                                          
                                          (type == "waste" & waste == "tomato") | 
                                            
                                          (experiment =="4" & system == "open" & type == "residue" & 
                                               
                                           waste == "tomato")
      
                                          )
  
  amp_waste_residue <- amp_venn(amp_waste_residue,group_by = "type",cut_a = 0.01,cut_f = 80,detailed_output = TRUE)
  
  
  amp_waste_residue$Otutable %>% filter(Shared=="Core") %>% View()
  
  
  ggsave("output/venn_tomato_residue_waste.jpeg", height = 16, width = 17, units = 'cm',dpi = 1200)
  

  
  #  food waste (change system)
  
  amp_waste_residue <- amp_subset_samples(amp_d, 
                                          
                                          (experiment == "1a" | experiment == "1b") |
                                            
                                            (type == "waste" & waste == "foodwaste") | 
                                            
                                            (experiment =="4" & system == "open" & type == "residue" & 
                                               
                                               waste == "foodwaste")
                                          
  )
  
  amp_waste_residue <- amp_venn(amp_waste_residue,group_by = "type",cut_a = 0.01,cut_f = 80,detailed_output = TRUE)
  
  
  amp_waste_residue$Otutable %>% filter(Shared=="Core") %>% View()
  
  
  ggsave("output/venn_foodwaste_residue.jpeg", height = 16, width = 17, units = 'cm',dpi = 1200)
