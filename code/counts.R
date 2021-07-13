##Analysis CFU

#load libraries ----

library(openxlsx)
library(tidyverse)
library(stringr)
library(hepre)
library(ggpubr)
library(gridExtra)

##cleans R brains, empties environment

rm(list=ls()) 

#load data from excel sheets -----


counts_raw <- read.xlsx(xlsxFile =here::here("data/plate_counts.xlsx"),sheet = 1)


#define columns as factors

experiment_fac <- c("0","3","4","5")

waste_fac <- c("tomato","wine","foodwaste")

treatment_fac <- c("control","waste_fly","wine","tomato","initital_waste",
                   "waste_inoculant","waste_inoculant_digestate","foodwaste")

type_fac <- c("waste","residue","inoculant")

systems <- c("closed","open")


counts <-

counts_raw %>% 
  
  select(experiment, waste, treatment, replicate, type,system,result_per_plate_CFUg,result_per_sample_CFUg) %>% 

  mutate(treatment=factor(treatment,levels=treatment_fac)) %>% 
  

  mutate(experiment=factor(experiment,levels=experiment_fac)) %>% 
  
  mutate(waste=factor(waste,levels=waste_fac)) %>% 
  
  mutate(system=factor(system,levels=systems)) %>% 
  
  mutate(type=factor(type,levels=type_fac)) %>% 
  
  gather(7:8,key = parameter,value = value,na.rm = TRUE) %>% 
  
  mutate(log_value=log10(value))




# CFU inoculants ----------------------------------------------------------

counts_inoc <-

counts %>% 
  
  filter(type=="inoculant") %>% 
  
  filter(experiment=="5") %>% 
  
  group_by(experiment,treatment,system,waste) %>% 
  
  filter(parameter=="result_per_sample_CFUg")

formatC(counts_inoc$value, format = "e", digits = 1)


# CFU wastes --------------------------------------------------------------

counts_wastes <-

counts %>% 
  
  filter(type=="waste") %>% 
  
  group_by(waste,replicate) %>% 
  
  filter(parameter=="result_per_plate_CFUg") %>% 
  
  summarise(n=n(),
            sd=sd(log_value),
            mean=mean(log_value))

formatC(counts_inoc$value, format = "e", digits = 1)

# CFU residues ----------------------------------------------------------

#mean

counts_residue_mean_exp4 <- 

  counts %>% 
    
    filter(type=="residue") %>% 
    
    filter(experiment=="4") %>% 
    
    group_by(experiment,waste,treatment,system) %>% 
    
    filter(parameter=="result_per_sample_CFUg") %>% 
    
    summarise(n=n(),
              sd=sd(log_value,na.rm = TRUE),
              mean=mean(log_value,na.rm = TRUE)) 

counts_residue_mean_exp5 <- 
  
  counts %>% 
  
  filter(type=="residue") %>% 
  
  filter(experiment=="5") %>%
  
  filter(!treatment=="waste_inoculant_digestate") %>% 
  
  group_by(experiment,waste,treatment) %>% 
  
  filter(parameter=="result_per_sample_CFUg") %>% 
  
  summarise(n=n(),
            sd=sd(log_value,na.rm = TRUE),
            mean=mean(log_value,na.rm = TRUE)) 


  # Experiment 4 ------------------------------------------------------------
  
  
  counts_residue_4 <-
  
  counts %>% 
    
    filter(experiment=="4") %>% 
    
    filter(type=="residue") %>% 
    
    group_by(experiment,waste,treatment,waste,system) %>% 
    
    filter(parameter=="result_per_sample_CFUg")
  
  #plot
  
  counts_residue_mean %>% 
    
    filter(experiment=="4") %>% 
    
    ggplot(aes(waste,mean)) +
    
    geom_point(aes(color=system)) +
    
    geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd),width = 0.3) +
    
    geom_jitter(data=counts_residue_4,
                aes(treatment,value,color=system),
                size=3,
                position = position_jitter(0.1),
                alpha=.4) +
    
    labs(y = "CFU/g residue", x="") +
  
    theme_bw(base_size = 11)  +
    theme(axis.text.x = element_text(angle = 0, colour = "black", size = 11))+
    theme(axis.text.y = element_text(angle = 0, colour = "black", size = 11))+
    theme(axis.title.y = element_text(angle = 90, colour = "black", size = 11))+
    theme(axis.title.x = element_text(angle = 0, colour = "black", size = 11))+
    theme(plot.title = element_text(angle = 0, colour = "black", size = 11))+
    theme(text=element_text(family="Arial")) +
    
    scale_x_discrete(limits=c("tomato",
                              "wine",
                              "foodwaste"),
                     labels=c("tomato"="Tomato pomace",
                              "wine" = "White wine pomace",
                              "foodwaste"="Pre-consumer food waste")) +
    
    coord_flip() +
    
    expand_limits(y = c(5, 10))
  
    
    
  #ggsave("Outputs/Plots/exp1_WR.jpeg", height = 8, width = 16, units = 'cm',dpi = 300)

  
  # Experiment 5 ------------------------------------------------------------
  
  
  counts_residue_5 <-
    
    counts %>% 
    
    filter(experiment=="5") %>% 
    
    filter(type=="residue") %>% 
    
    group_by(experiment,waste,treatment,waste,system) %>% 
    
    filter(parameter=="result_per_sample_CFUg")
  
  #plot
  
  counts_residue_mean %>% 
    
    filter(experiment=="5") %>% 
    
    ggplot(aes(treatment,mean)) +
    
    geom_point(aes(color=system)) +
    
    geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd),width = 0.1) +
    
    geom_jitter(data=counts_residue_5,
                aes(treatment,value,color=system),
                size=3,
                position = position_jitter(0.1),
                alpha=.4) +
    
    facet_wrap(~waste) +
    
    labs(y = "CFU/g residue", x="") +
    
    coord_flip() +
    
    theme_bw(base_size = 12)  +
    theme(axis.text.x = element_text(angle = 0, colour = "black", size = 11))+
    theme(axis.text.y = element_text(angle = 0, colour = "black", size = 11))+
    theme(axis.title.y = element_text(angle = 90, colour = "black", size = 11))+
    theme(axis.title.x = element_text(angle = 0, colour = "black", size = 11))+
    theme(plot.title = element_text(angle = 0, colour = "black", size = 11))+
    theme(text=element_text(family="Arial")) +
    
    expand_limits(y = c(5, 10)) +
    
    scale_x_discrete(limits=c("waste_inoculant_digestate",
                              "waste_inoculant",
                              "control"),
                     labels=c("waste_inoculant_digestate"="Food waste inoculant",
                              "waste_inoculant" = "Residue\ninoculant",
                              "control"="Control\n(autoclaved residue inoculant)"))
  
  
  ggsave("Outputs/Plots/exp5_counts.jpeg", height = 8, width = 16, units = 'cm',dpi = 300)
  