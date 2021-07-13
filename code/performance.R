##Analysis performance data

##cleans R brains, empties environment

rm(list=ls())

#load libraries ----

library(openxlsx)
library(tidyverse)
library(stringr)
library(here)
library(ggpubr)
library(gridExtra)
library(moments)
library(rstatix)
library(Hmisc)
library(gtable)
library(lemon)
library(conflicted)
library(ungeviz)

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("intersect", "dplyr")
conflict_prefer("lag", "dplyr")
conflict_prefer("mutate","dplyr")
conflict_prefer("summarise","dplyr")
conflict_prefer("rename","dplyr")

#load data from excel sheets -----

performance_raw <- read.xlsx(xlsxFile =here::here("Data/performance_data.xlsx"),sheet = 1)


#define columns as factors

systems <- c("closed","open")

performance <- performance_raw %>% 
  
  mutate(system=factor(system,levels=systems))


# calculate BCR ----

performance <- performance %>% 
  
  spread(key = parameter,value = value) %>% 
  
  mutate(BCR=(((larval_number*larval_weight_mgDM)/1000)/60)*100) %>% 
  
  select(-larval_number) %>% 
  
  group_by(experiment,waste,treatment,replicate) %>% 
  
  rename(`Larval mass (mg DM)` = larval_weight_mgDM) %>% 
  
  rename(`Substrate reduction (% DM)` = waste_reduction_percDM) %>% 
  
  rename(`Bioconversion rate (% DM)` = BCR) %>% 
  
  rename(`Residue pH (-)`= residue_pH) %>% 
  
  rename(`Residue MC (%)`= residue_MC) %>% 
  
  gather(6:10,key = parameter,value = value,na.rm = FALSE) %>% 

  mutate(parameter=factor(parameter))


# mean all performance metrics/samples ----

 performance_mean <- 

 performance %>% 
   
   group_by(system,waste,parameter,experiment,treatment) %>% 
   
   summarise(n=n(),
             mean=mean(value,na.rm = TRUE),
             sd=sd(value,na.rm=TRUE))


# plot results ----
  
  # mean performance metrics excluding BCR

  performance %>% 
  
    filter(!parameter=="Bioconversion rate (% DM)") %>% 
    
    ungroup(system) %>% 
    
    group_by(parameter,waste,experiment,treatment) %>% 
    
    summarise(n=n(),
              mean=mean(value),
              sd=sd(value))


  # experiment 4 ----


  # set order of facets

  performance$parameter <- factor(performance$parameter, levels = c("Larval mass (mg DM)",
                                                                    "Substrate reduction (% DM)", "Residue MC (%)",
                                                                    "Residue pH (-)", "Bioconversion rate (% DM)"))

  performance_mean$parameter <- factor( performance_mean$parameter, levels = c("Larval mass (mg DM)",
                                                                        "Substrate reduction (% DM)", "Residue MC (%)","Residue pH (-)", "Bioconversion rate (% DM)"))

  # round value on axis
  
  f <- function(x){
    format(round(x, 0), nsmall=0)
  }

  
  # all except for MC ----
  

  # save data for plot  

  performance_4 <- 
  
  performance %>% 
    
    filter(experiment=="4") %>% 
    
    mutate(waste=case_when(waste=="tomato"~"TP",
                           waste=="foodwaste"~"DFW",
                           waste=="wine"~"WWP",
                           TRUE~"")) %>%
    
    filter(!parameter=="Bioconversion rate (% DM)") %>% 
    
    filter(!parameter=="Residue MC (%)")
  
  
  # plot 

  performance_mean %>% 
  
  filter(experiment=="4") %>% 
  
  filter(!parameter=="Bioconversion rate (% DM)") %>% 
  
  filter(!parameter=="Residue MC (%)") %>% 
  
  mutate(waste=case_when(waste=="tomato"~"TP",
                         waste=="foodwaste"~"DFW",
                         waste=="wine"~"WWP",
                         TRUE~"")) %>%
  
  ggplot(aes(waste,mean)) +
  
  scale_shape_manual(values=c(16, 21),name="Rearing\nsystem") +

  geom_hpline(data=performance_4,
              aes(waste, value, shape=system),
              position = position_dodge(0.4),
              width = 0.2, size = 0.7,
              stat = "summary") +
  
  geom_jitter(data=performance_4,
              aes(waste, value, shape=system),
              position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.4),
              size = 3.0
  ) +
    
  stat_summary(data=performance_4,
               aes(waste, value, shape = system),
               fun.data="mean_sdl",  fun.args = list(mult=1), 
               geom = "errorbar",  size = 0.4, width=.2,
               position = position_dodge(0.4)
  ) +
  
  theme_bw(base_size = 15) +
  
  ylab("") +
  
  xlab("") +
  
  facet_wrap(~parameter,scales = "free_y")

  ggsave("output/performance_exp4.jpeg", height = 8, width = 24, units = 'cm',dpi = "print")


  
  # only MC ----
  
  # save data
  
  performance_4_MC <- 
    
    performance %>% 
    
    filter(experiment=="4") %>% 
    
    mutate(experiment="experiment 1") %>% 
    
    mutate(experiment=case_when(experiment=="experiment 1"~"Experiment 1",
                                TRUE~"")) %>% 
    
    mutate(waste=case_when(waste=="tomato"~"TP",
                           waste=="foodwaste"~"DFW",
                           waste=="wine"~"WWP",
                           TRUE~"")) %>%
    
    filter(parameter=="Residue MC (%)")
  
  
  # plot

  performance_mean %>% 
  
  filter(experiment=="4") %>% 
  
  mutate(experiment="experiment 1") %>% 
    
  mutate(experiment=case_when(experiment=="experiment 1"~"Experiment 1",
                                 TRUE~"")) %>% 
  
  filter(parameter=="Residue MC (%)") %>% 
  
  mutate(waste=case_when(waste=="tomato"~"TP",
                         waste=="foodwaste"~"DFW",
                         waste=="wine"~"WWP",
                         TRUE~"")) %>% 
  
  summarise(mean=mean(mean)) %>% 

  
  ggplot(aes(waste,mean)) +
  
    geom_hpline(data=performance_4_MC,
                aes(waste, value, shape=system),
                position = position_dodge(0.3),
                width = 0.2, size = 0.7,
                stat = "summary") +
  
  geom_jitter(data=performance_4_MC, 
              aes(waste, value, shape = system), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.3),
              size = 3.0
  ) +
  
  stat_summary(data=performance_4_MC,
               aes(waste, value, shape = system),
               fun.data="mean_sdl",  fun.args = list(mult=1), 
               geom = "errorbar",  size = 0.4, width=.2,
               position = position_dodge(0.3)
  ) +
  
  theme_bw(base_size = 15) +
  
  ylab("Moisture content (%)") +
  
  xlab("") +
  
  facet_wrap(~experiment) + 
  
  scale_y_continuous(labels=f) +
  
  scale_shape_manual(values=c(16, 1),name="Rearing\nsystem") +
  
  ylim(30,90)

  ggsave("output/performance_exp4_MC.jpeg", height = 8, width = 11, units = 'cm',dpi = "print")


  
  # experiment 5 ----

  #all except for MC

  # save data
  
  performance_5 <- 
    
    performance %>% 
    
    filter(experiment=="5") %>% 
    
    filter(!treatment == "tomato (food waste inoculant)") %>% 
    
    mutate(treatment=case_when(treatment=="wine (sterile inoculant)"~"WWP (o)",
                               treatment=="wine (inoculant)"~"WWP (+)",
                               treatment=="tomato (sterile inoculant)"~"TP (o)",
                               treatment=="tomato (inoculant)"~"TP (+)",
                               TRUE~"")) %>% 
    
    filter(!parameter=="Bioconversion rate (% DM)") %>% 
    
    filter(!parameter=="Residue MC (%)")
  
  # plot
  
  performance_5 <- 
    
    performance_mean %>% 
    
    filter(experiment=="5") %>% 
    
    filter(!parameter=="Bioconversion rate (% DM)") %>% 
    
    filter(!parameter=="Residue MC (%)") %>% 
    
    filter(!treatment == "tomato (food waste inoculant)") %>% 
    
    mutate(treatment=case_when(treatment=="wine (sterile inoculant)"~"WWP (o)",
                               treatment=="wine (inoculant)"~"WWP (+)",
                               treatment=="tomato (sterile inoculant)"~"TP (o)",
                               treatment=="tomato (inoculant)"~"TP (+)",
                               TRUE~"")) %>% 
    
    ggplot(aes(treatment,mean)) +
    
    geom_hpline(data=performance_5,
                aes(treatment, value, shape=system),
                position = position_dodge(0.4),
                width = 0.2, size = 0.7,
                stat = "summary") +
    
    geom_jitter(data=performance_5, 
                aes(treatment, value, shape = system), 
                position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.4),
                size = 3.0
    ) +
    
    stat_summary(data=performance_5,
                 aes(treatment, value, shape = system),
                 fun.data="mean_sdl",  fun.args = list(mult=1), 
                 geom = "errorbar",  size = 0.4, width=.2,
                 position = position_dodge(0.4)
    ) +
    
    theme_bw(base_size = 15) +
    
    ylab("") +
    
    xlab("") +
    
    facet_wrap(~parameter,scales = "free_y",nrow = 2,ncol = 2) +
    
    geom_vline(xintercept = c(2.5), linetype="dashed") +
    
    scale_shape_manual(values=c(16, 21),name="Rearing\nsystem")
  
  
  # for defense
  
  performance_5_slim <- 
    
    performance %>% 
    
    filter(experiment=="5") %>% 
    
    filter(!treatment == "tomato (food waste inoculant)") %>% 
    
    mutate(treatment=case_when(treatment=="wine (sterile inoculant)"~"WWP (o)",
                               treatment=="wine (inoculant)"~"WWP (+)",
                               treatment=="tomato (sterile inoculant)"~"TP (o)",
                               treatment=="tomato (inoculant)"~"TP (+)",
                               TRUE~"")) %>% 
    
    filter(!parameter=="Bioconversion rate (% DM)") %>% 
    
    filter(!parameter=="Residue MC (%)") %>% 
    
    filter(parameter=="Larval mass (mg DM)") %>% 
    
    filter(system=="closed")
  
  
  performance_mean %>% 
    
    filter(experiment=="5") %>% 
    
    filter(!parameter=="Bioconversion rate (% DM)") %>% 
    
    filter(!parameter=="Residue MC (%)") %>% 
    
    filter(!treatment == "tomato (food waste inoculant)") %>% 
    
    filter(parameter=="Larval mass (mg DM)") %>% 
    
    filter(system=="closed") %>% 
    
    mutate(treatment=case_when(treatment=="wine (sterile inoculant)"~"WWP (o)",
                               treatment=="wine (inoculant)"~"WWP (+)",
                               treatment=="tomato (sterile inoculant)"~"TP (o)",
                               treatment=="tomato (inoculant)"~"TP (+)",
                               TRUE~"")) %>% 
    
    ggplot(aes(treatment,mean)) +
    
    geom_point(size=5)+
    
    geom_jitter(data=performance_5_slim, 
                aes(treatment, value), size=2.5) +
    
    stat_summary(data=performance_5_slim,
                 aes(treatment, value),
                 fun.data="mean_sdl",  fun.args = list(mult=1), 
                 geom = "errorbar",  size = 0.4, width=.2,
                 position = position_dodge(0.3)
    ) +
    
    theme_bw(base_size = 15) +
    
    ylab("Larval weight\n[mg DM / larva]") +
    
    xlab("") +
  
    geom_vline(xintercept = c(2.5), linetype="dashed") 
  
  ggsave("output/performance_exp5_defense.jpeg", height = 14, width = 22, units = 'cm',dpi = "print")
  
  
  
    # shift legend to free space in plot
  
  performance_5_legend <- shift_legend2(performance_5)
  
  ggsave("output/performance_exp5.jpeg", height = 16, width = 19, units = 'cm',performance_5_legend)
  
  
  shift_legend2 <- function(p) {
    # check if p is a valid object
    if(!(inherits(p, "gtable"))){
      if(inherits(p, "ggplot")){
        gp <- ggplotGrob(p) # convert to grob
      } else {
        message("This is neither a ggplot object nor a grob generated from ggplotGrob. Returning original plot.")
        return(p)
      }
    } else {
      gp <- p
    }
    
    # check for unfilled facet panels
    facet.panels <- grep("^panel", gp[["layout"]][["name"]])
    empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]), 
                                 USE.NAMES = F)
    empty.facet.panels <- facet.panels[empty.facet.panels]
    
    if(length(empty.facet.panels) == 0){
      message("There are no unfilled facet panels to shift legend into. Returning original plot.")
      return(p)
    }
    
    # establish name of empty panels
    empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
    names <- empty.facet.panels$name
    
    # return repositioned legend
    reposition_legend(p, 'center', panel=names)
  } 
  
  
  # only MC ----
  
  # save data
  
  performance_5_MC <- 
    
    performance %>% 
    
    filter(experiment=="5") %>% 
    
    mutate(experiment="experiment 2") %>% 
    
    mutate(experiment=case_when(experiment=="experiment 2"~"Experiment 2",
                                TRUE~"")) %>% 
    
    filter(!treatment == "tomato (food waste inoculant)") %>% 
    
    mutate(treatment=case_when(treatment=="wine (sterile inoculant)"~"WWP (o)",
                               treatment=="wine (inoculant)"~"WWP (+)",
                               treatment=="tomato (sterile inoculant)"~"TP (o)",
                               treatment=="tomato (inoculant)"~"TP (+)",
                               TRUE~"")) %>% 

    filter(parameter=="Residue MC (%)")
  
  
    # plot

    performance_mean %>% 
    
    filter(experiment=="5") %>% 
      
    mutate(experiment="experiment 2") %>% 
      
      mutate(experiment=case_when(experiment=="experiment 2"~"Experiment 2",
                                  TRUE~"")) %>% 

    filter(parameter=="Residue MC (%)") %>% 
    
    filter(!treatment == "tomato (food waste inoculant)") %>% 
    
    mutate(treatment=case_when(treatment=="wine (sterile inoculant)"~"WWP (o)",
                               treatment=="wine (inoculant)"~"WWP (+)",
                               treatment=="tomato (sterile inoculant)"~"TP (o)",
                               treatment=="tomato (inoculant)"~"TP (+)",
                               TRUE~"")) %>% 
    
    ggplot(aes(treatment,mean)) +
    
      geom_hpline(data=performance_5_MC,
                  aes(treatment, value, shape=system),
                  position = position_dodge(0.3),
                  width = 0.2, size = 0.7,
                  stat = "summary") +
    
    geom_jitter(data=performance_5_MC, 
                aes(treatment, value, shape = system), 
                position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.3),
                size = 3.0
    ) +
    
    stat_summary(data=performance_5_MC,
                 aes(treatment, value, shape = system),
                 fun.data="mean_sdl",  fun.args = list(mult=1), 
                 geom = "errorbar",  size = 0.4, width=.2,
                 position = position_dodge(0.3)
    ) +
    
    theme_bw(base_size = 15) +
    
    ylab("Moisture content (%)") +
    
    xlab("") +
    
    facet_wrap(~experiment,scales = "free_y",nrow = 2,ncol = 2) +
    
    geom_vline(xintercept = c(2.5), linetype="dashed") +
    
    scale_y_continuous(labels=f) +
    
    scale_shape_manual(values=c(16, 1),name="Rearing\nsystem") +
      
      theme(legend.position = "none") +
      
      ylim(30,90)
  
  ggsave("output/performance_exp5_MC.jpeg", height = 8, width = 11, units = 'cm')
  