##Analysis temperature data

#load libraries -----

library(openxlsx)
library(tidyverse)
library(stringr)
library(here)
library(ggpubr)
library(gridExtra)
library(readr)
library(fs)
library(stringr)
library(lubridate)
library(conflicted)


conflict_prefer("summarize", "dplyr")
conflict_prefer("hour", "lubridate")
conflict_prefer("year", "lubridate")
conflict_prefer("month", "lubridate")
conflict_prefer("mutate","dplyr")
conflict_prefer("summarise","dplyr")


#cleans R brains, empties environment -----

rm(list=ls()) 

#import data ----

#temperature labelling syntax to allow smooth importing

  #experiment = Tprod, Wprod, Tinoc, Winoc
  #waste = T, D, W (in Tinoc, D is the the type of inoculant)
  #replicate = 1-4
  #system = open/close
  #type = control/treatment


  experiments <- c("Tprod","Wprod","Tinoc","Winoc")
  wastes <- c("T","D","W")
  replicates <- c("1","2","3","4")
  systems <- c("open","closed")
  types <- c("control","treatment")
  
#set directory where files are

data_dir <- "data/temperature/"

#make list of all file names

csv_temp <- fs::dir_ls(data_dir)

#import all files and save as tibble

temperature_messy <-

  csv_temp %>% 
  
  #add column that shows the name of the files from which the data was imported
  
  map_dfr(read_csv,skip=18,.id="experiment_replicate")


#change the name of the column names

names(temperature_messy) <- c("experiment_replicate","date_time","unit","value")


#seperate the experiment_replicate and data_time columns

  temperature <- 

  temperature_messy %>% 
  
    separate(experiment_replicate, c("experiment","waste","replicate","system","type"),sep="_") %>% 
    
    mutate(experiment = str_replace(experiment, "data/temperature/","")) %>% 
    
    mutate(type = str_replace(type, ".csv","")) %>% 
    
    #convert date_time into POSIXct date-time object
    
    mutate(date_time=parse_date_time(date_time,"%m/%d/%y %I:%M:%S %p")) %>% 
    
    #remove unit - all data in °C
    
    select(-unit) %>% 
    
    mutate(experiment=factor(experiment,levels=experiments)) %>% 
    
    mutate(waste=factor(waste,levels=wastes)) %>% 
    
    mutate(replicate=factor(replicate,levels=replicates)) %>% 
    
    mutate(system=factor(system,levels=systems)) %>% 
    
    mutate(type=factor(type,levels=types)) 
    
  
  
# summarize data by hour ----
  
  # create columns for year, month, day hour
  
  temperature$year <- year(temperature$date_time)
  temperature$month  <- month(temperature$date_time)
  temperature$day  <- day(temperature$date_time)
  temperature$hour <- hour(temperature$date_time)
  
  # group data by year, month, day hour, waste, system, experiment etc.
  
  temperature <- 
    
    temperature %>% 
    
    group_by(experiment,waste,replicate,system,type,year,day,hour,month)
  
  
  # confirm that the grouping works in the right way
  
  tally(temperature)
  
  
  # calculate the mean temperature per hour
  
  temperature <-
    
    summarize(temperature,value=mean(value))
  
  # create time stamp
  
  temperature$date_time <- with(temperature, ymd_h(paste(year, month, day, hour,sep= ' ')))
  

  
  
  
#set beginning of each experiment ----
    
    Tprod_beginning <- c("12/4/19 8:00:00 PM") #noted in lab book
    Wprod_beginning <- c("12/9/19 6:00:00 PM") #note in lab book
    Tinoc_beginning <- c("12/12/19 8:00:00 PM") #estimated basen on trajectory of curve
    Winoc_beginning <- c("12/21/19 10:00:00 PM") #estimated basen on trajectory of curve

    
#change absolute date and timestamp to experimental duration ----
    
    temperature_beginning <- 
    
    temperature %>% 
      
      mutate(date_time_beginning=case_when(experiment=="Tprod"~Tprod_beginning,
                                           TRUE~"")) %>% 
      
      mutate(date_time_beginning=case_when(experiment=="Wprod"~Wprod_beginning,
                                           TRUE~date_time_beginning)) %>% 
    
      mutate(date_time_beginning=case_when(experiment=="Tinoc"~Tinoc_beginning,
                                         TRUE~date_time_beginning)) %>% 
    
      mutate(date_time_beginning=case_when(experiment=="Winoc"~Winoc_beginning,
                                         TRUE~date_time_beginning)) %>% 
    
      mutate(date_time_beginning=parse_date_time(date_time_beginning,"%m/%d/%y %I:%M:%S %p")) %>% 
      
      #only show values that are after the beginning of the experiment
      
      filter(date_time>date_time_beginning) 
     
      
      
      temperature_exp_duration <-   
      
        temperature_beginning %>% 
         
      #calculate experimental duration
      
      mutate(experimental_duration_sec=date_time-date_time_beginning) %>% 
      
      mutate(experimental_duration_day=(as.numeric(experimental_duration_sec)/(60*60*24))) %>% 
      
      select(-experimental_duration_sec) %>% 
      
      #set end of the experiment
      
      filter((experiment=="Tprod" & experimental_duration_day<5.50) |
             (experiment=="Tinoc" & experimental_duration_day<5.55) |
             (experiment=="Wprod" & experimental_duration_day<9.60) |
             (experiment=="Winoc" & experimental_duration_day<8.55))
        


# calculate mean of replicates -----
      
      temperature_rep <- temperature_exp_duration %>% 
        
        select(-date_time_beginning,-experimental_duration_day)
      
      # group data by year, month, day hour, waste, system, experiment etc.
      
      temperature_rep <- 
        
        temperature_rep %>% 
        
        ungroup() %>% 
        
        group_by(experiment,waste,system,type,year,day,hour,month)
      
      
      # confirm that the grouping works in the right way
      
      tally(temperature_rep)
      
      
      # calculate the mean temperature per hour
      
      temperature_rep <-
        
        summarize(temperature_rep,value=mean(value))
      
      # create time stamp
      
      temperature_rep$date_time <- with(temperature_rep, ymd_h(paste(year, month, day, hour,sep= ' ')))
      
      # convert time stamp back to experimental duration
      
      temperature_duration_rep_beginning <- 
      
        temperature_rep %>% 
        
        mutate(date_time_beginning=case_when(experiment=="Tprod"~Tprod_beginning,
                                             TRUE~"")) %>% 
        
        mutate(date_time_beginning=case_when(experiment=="Wprod"~Wprod_beginning,
                                             TRUE~date_time_beginning)) %>% 
        
        mutate(date_time_beginning=case_when(experiment=="Tinoc"~Tinoc_beginning,
                                             TRUE~date_time_beginning)) %>% 
        
        mutate(date_time_beginning=case_when(experiment=="Winoc"~Winoc_beginning,
                                             TRUE~date_time_beginning)) %>% 
        
        mutate(date_time_beginning=parse_date_time(date_time_beginning,"%m/%d/%y %I:%M:%S %p")) %>% 
        
        #only show values that are after the beginning of the experiment
        
        filter(date_time>date_time_beginning) 
      
      
      temperature_duration_rep <-
      
        temperature_duration_rep_beginning %>% 
        
        #calculate experimental duration
        
        mutate(experimental_duration_sec=date_time-date_time_beginning) %>% 
        
        mutate(experimental_duration_day=(as.numeric(experimental_duration_sec)/(60*60*24))) %>% 
        
        select(-experimental_duration_sec) %>% 
        
        #set end of the experiment
        
        filter((experiment=="Tprod" & experimental_duration_day<5.50) |
                 (experiment=="Tinoc" & experimental_duration_day<5.30) |
                 (experiment=="Wprod" & experimental_duration_day<9.60) |
                 (experiment=="Winoc" & experimental_duration_day<8.00))
      
      
   
      
      
#calculate descriptive statistics ----
    
    temperature_descriptive <- 
    
    temperature_exp_duration %>% 
      
      group_by(experiment, waste,system,type) %>% 
      
      summarise(n=n(),
                mean=mean(value),
                median=median(value),
                sd=sd(value),
                max=max(value),
                min=min(value))
    
    write.xlsx(temperature_descriptive, file="Outputs/temperature_descriptive.xlsx")

    
#open vs. closed system ----

    
    # experiment 4 ----
    
    
    f <- function(x){
      format(round(x, 0), nsmall=0)
    }
    
    
    
    temperature_descriptive_exp4 <- 
    
      temperature_descriptive %>% 
    
      mutate(waste=case_when(waste=="T"~"TP",
                             waste=="D"~"DFW",
                             waste=="W"~"WWP",
                             TRUE~"")) %>% 
      
      filter(experiment=="Tprod" | experiment =="Wprod") 
    
      
    
    temperature_exp_duration %>% 
      
      mutate(waste=case_when(waste=="T"~"TP",
                             waste=="D"~"DFW",
                             waste=="W"~"WWP",
                             TRUE~"")) %>% 
      
      filter(experiment=="Tprod" | experiment =="Wprod") %>% 
      
      # filter(waste=="T") %>% 
      
      ggplot(aes(experimental_duration_day,value)) +
      
      # geom_point() +
      
      geom_line(aes(group=interaction(replicate, system),
                    linetype=system)) +
      
      facet_wrap(~waste,scales = "free_x") +
      
      theme_bw(base_size = 15) +
      
      ylab("Residue temperature (°C)") +
      
      xlab("Rearing duration (days)") +
      
       scale_x_continuous(labels=f) + 
    
      scale_linetype_discrete(name = "Rearing\nsystem") +
      
      geom_hline(data=temperature_descriptive_exp4,aes(yintercept=median))
    
    
    ggsave("output/temperature_exp4.jpeg", height = 10, width = 20, units = 'cm',dpi = "print")  
    
   
    temperature_descriptive %>% 
      
      filter(experiment=="Tprod" | experiment =="Wprod") %>% 
      
      ggplot(aes(system,mean)) + 
      
      geom_point() +
      
      facet_wrap(~waste) +
      
      geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd),width = 0.3) +
      
      theme_bw(base_size = 15)
    
    
    # experiment 5 ----
  
      # calculate descriptive statics for exp 5
    
    
    temperature_descriptive_exp5 <- 
      
      temperature_descriptive %>% 
      
      filter(!waste=="D") %>% 
      
      mutate(waste=case_when(waste=="T"~"TP",
                             waste=="W"~"WWP",
                             TRUE~"")) %>% 
      
      mutate(type=case_when(type=="control"~"sterile inoculant",
                            type=="treatment"~"inoculant",
                            TRUE~"")) %>% 
      
      filter(experiment=="Tinoc" | experiment =="Winoc")
    
    
    temperature_duration_rep %>% 
      
      filter(!waste=="D") %>% 
      
      mutate(waste=case_when(waste=="T"~"TP",
                             waste=="W"~"WWP",
                             TRUE~"")) %>% 
      
      mutate(type=case_when(type=="control"~"sterile inoculant",
                             type=="treatment"~"inoculant",
                             TRUE~"")) %>% 
      
      filter(experiment=="Tinoc" | experiment =="Winoc") %>% 
      
      # filter(waste=="T") %>% 
      
      ggplot(aes(experimental_duration_day,value)) +
      
      # geom_point() +
      
      geom_line(aes(group=interaction(system, type),
                    color=type,
                    linetype=system)) +
      
      facet_wrap(~waste,scales = "free_x") +
      
      theme_bw(base_size = 15) +
      
      ylab("Residue temperature (°C)") +
      
      xlab("Rearing duration (days)") +
      
      scale_x_continuous(labels=f) + 
      
      scale_linetype_discrete(name = "Rearing\nsystem") +
      
      scale_color_discrete(name = "Inoculant type")
      
      # geom_hline(data=temperature_descriptive_exp5,aes(yintercept=median))
    
    ggsave("output/temperature_exp5.jpeg", height = 10, width = 16, units = 'cm',dpi = 300)  
    
    
    temperature_descriptive %>% 
      
      filter(experiment=="Tinoc" | experiment =="Winoc") %>% 
      
      filter(type=="control") %>% 
      
      ggplot(aes(system,mean)) + 
      
      geom_point(aes(color=type)) +
      
      facet_wrap(~waste) +
      
      geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd),width = 0.3) +
      
      theme_bw(base_size = 15) +
      
      labs(y = "mean temperature °C", x="")
 
    