#script for reading in TRIP data & reorganize and manipulate
#version to read in multiple files for multiple constructs/experiments
#loop through all files and carry out normalizing/scaling procedures

library(tidyr)
library(dplyr)
library(tidyverse)
library(plotrix)
library(ggplot2)
library(viridis)
library(forcats)
library(scales)
library(gplots)


#input file names and experiments
#read in csv file with
#1. columns listing csv files (header: "File"),
#2. column listing conditions (header: "Experiment")
file_names<-read_csv("TRIP-layout.csv")

#run through loop carrying out all the analysis for each file

#create empty tibbles for data_storage
#
TRIP_long_HPGscaled<-tibble()

#test loop
#for (i in file_names$File) {
#  print(i)
#  experi=file_names %>% filter(File == i) %>% select(Experiment) %>% as.character()
#  print(experi)
#}

#loop through individual files and analyze
for (i in file_names$File) {
    
  #extract experiment names
  experi=file_names %>% filter(File == i) %>% select(Experiment) %>% as.character()
  
  #Read data
  data<-read_csv(i)
  
  # Select relevant columns
  pd_prot_selected <- data %>%
    select(
      "Accession",
      "Description",
      "Gene Symbol",
      starts_with("Abundance:"), # select all columns that start with "Abundance:"
    )
  
  # Convert into long format
  pd_prot_long <- pd_prot_selected %>%
    pivot_longer(
      cols = starts_with("Abundance:"),
      names_to = "sample",
      values_to = "intensity"
    )
  
  #read in list of Tg TRIP interactors (based on comparion to (-)Biotin comparison)
  interactors<-read_csv("../Tg-TRIP-interactors.csv")
  #convert to all uppercase
  interactors$interactor<-toupper(interactors$interactor)
  
  #filter based on prior Tg TRIP interactors
  pd_prot_long<-pd_prot_long %>% filter(`Gene Symbol` %in% c(interactors$interactor))
  
  #Separate Sample annotations into separate columns
  pd_prot_long_sep <- pd_prot_long %>% separate_wider_delim(sample, delim = ",", names = c("TMT","type","rep","condition"))
  
  #Subset only 0-3 Hr course data
  pd_prot_long_TC <- pd_prot_long_sep %>%
    filter(condition %in% c(" 0_0Hr"," 0_5Hr", " 1_0Hr", " 1_5Hr"," 2_0Hr", " 3_0Hr"))
  
  #Tg bait normalization - individual per rep
  # Subset Tg row to obtain mean reference values
  rep_means <- pd_prot_long_TC %>% filter(Accession == "P01266") %>% select(condition, rep, intensity) %>%
    group_by(rep) %>%
    summarize(rep_mean = mean(intensity))
  
  # Subset Tg row to obtain reference values
  conditions_Tg <- pd_prot_long_TC %>% filter(Accession == "P01266") %>% select(condition, rep, intensity) %>%
    group_by(rep) 
  
  # Add Tg reference intensities
  pd_prot_long_norm <- pd_prot_long_TC %>%
    # Group by conditions and rep
    group_by(condition, rep) %>%
    # Add the mean Tg for each rep
    mutate(Tg_mean_int_rep = rep_means$rep_mean[match(rep, rep_means$rep)])
    # Add the Tg intensity for each rep and conditions
  
  #Normalize each rep separately
  pd_prot_long_norm<-left_join(pd_prot_long_norm,conditions_Tg,by=c('rep'='rep', 'condition'='condition')) %>%
    # Apply the normalization factor to the value column
    mutate(normalized_value = intensity.x * Tg_mean_int_rep / intensity.y)
  
  #pivot to wide format and export to csv
  pd_prot_norm_wide<-pd_prot_long_norm %>%
    #order based on rep
    arrange(rep) %>%
    #pivot to wide format
    pivot_wider(
      id_cols=c("Accession","Description","Gene Symbol"),
      names_from = c("TMT","condition","rep"),
      values_from="normalized_value")
  filename_TgNorm=paste0(experi,"-TRIP-PD2-Sequest-TgNorm.csv")
  write_csv(pd_prot_norm_wide, filename_TgNorm)
  
  #Generate merged table of TRIP Time-course data normalized to Tg, and unnnormalized (-)Hpg, (-)Biotin 
  #extract tibble with control values, expand tibble to same variables
  pd_prot_control<-filter(pd_prot_long_sep, condition %in% c(" (-)Hpg"," (-)BiotinPD")) %>%
    mutate(Tg_mean_int_rep=c(NA)) %>%
    mutate(intensity.y=c(NA)) %>%
    mutate(normalized_value=intensity) %>%
    rename(intensity.x=intensity)
    
  pd_prot_long_TgNorm_full<-full_join(pd_prot_long_norm, pd_prot_control)
  
  #pivot to wide format and export to csv
  pd_prot_TgNorm_full_wide<-pd_prot_long_TgNorm_full %>%
    #order based on rep
    arrange(rep) %>%
    #pivot to wide format
    pivot_wider(
      id_cols=c("Accession","Description","Gene Symbol"),
      names_from = c("TMT","condition","rep"),
      values_from="normalized_value")
  filename_TgNorm_Full=paste0(experi,"-TRIP-PD2-Sequest-TgNorm-Full.csv")
  write_csv(pd_prot_TgNorm_full_wide, filename_TgNorm_Full)
  
  #Carry out normalization of TRIP Time Course data to (-)HPG background:
  # Subset (-)Hpg control conditions
  conditions_noHPG <- pd_prot_long_sep %>% filter(condition ==" (-)Hpg") %>%
    select(Accession, rep, intensity)
  #join add (-) HPG data to TRIP data tibble as extra column
  pd_prot_long_Bgnorm<-left_join(
    pd_prot_long_norm,conditions_noHPG,
    by=c('rep'='rep',"Accession"="Accession")) %>%
    rename(noHPG_intensity=intensity) %>%
    #Add column with noHPG normalized intensities
    mutate(HPG_norm_value = normalized_value / noHPG_intensity) %>%
    #replace HPG normalized ratios < 1 (below background) with 1
    mutate(HPG_norm_value = ifelse(HPG_norm_value < 1, 1, HPG_norm_value))
  
  #Scale the (-)HPG normalized data from 0-1 within each replicate:
  #Find the max enrichment over (-HPG) for each replicate
  rep_max_tib <- pd_prot_long_Bgnorm %>% select(Accession, condition, rep, HPG_norm_value) %>%
    group_by(Accession, rep) %>%
    summarize(rep_max = max(HPG_norm_value))
  # Add max enrichment intensities
  pd_prot_long_Bgnorm2 <- pd_prot_long_Bgnorm %>%
    # Add the max intensity for each rep
    left_join(rep_max_tib, by=c('rep'='rep',"Accession"="Accession")) %>%
    #Drop values to NA where max intensity <=1 (no enrichment above background at any time point)
    mutate(rep_max = ifelse(rep_max <= 1, NA, rep_max)) %>%
    #Add column with scaled intensities:
    mutate(scaled_intensity = log2(HPG_norm_value) / log2(rep_max))
  
  #pivot to wide format and export to csv
  pd_prot_wide_HPGscaled<-pd_prot_long_Bgnorm2 %>%
    #order based on rep
    arrange(condition) %>%
    #pivot to wide format
    pivot_wider(
      id_cols=c("Accession","Description","Gene Symbol"),
      names_from = c("TMT","condition","rep"),
      values_from=scaled_intensity)
  filename_HPGscaled_Full=paste0(experi,"-TRIP-PD2-Sequest-HPGscaled-Full.csv")
  write_csv(pd_prot_wide_HPGscaled, filename_HPGscaled_Full)
  
  
  #join data tibbles
  #if the prior tibble is empty, replace the tibble and add extra column with experiment descriptor
  if (nrow(TRIP_long_HPGscaled) == 0) {
    TRIP_long_HPGscaled<-pd_prot_long_Bgnorm2 %>% mutate(experiment=experi)
  } else {
    # if the existing tibble contains data, bind the new tibble and add experiment descriptor
    joined_tibble <- bind_rows((TRIP_long_HPGscaled),
                               mutate(pd_prot_long_Bgnorm2, experiment = experi))
    TRIP_long_HPGscaled<-joined_tibble
  }

}

#Calculate summary statistics for TRIP time course
TRIP_HPGscaled_sum <- TRIP_long_HPGscaled %>%
  group_by(Accession,Description,`Gene Symbol`,condition,experiment) %>%
  summarize(mn_scaled_intensity = mean(scaled_intensity, na.rm=TRUE),
            mn_scaled_sem = std.error(scaled_intensity, na.rm=TRUE)
  )

#pivot HPG_scaled data to wide format and export to csv
TRIP_HPGscaled_wide<-TRIP_HPGscaled_sum %>%
  #order based on experiment
  arrange(experiment) %>%
  #pivot to wide format
  pivot_wider(
    id_cols=c("Accession","Gene Symbol"),
    names_from = c("experiment","condition"),
    values_from=c(mn_scaled_intensity, mn_scaled_sem))
#export to csv file
write_csv(TRIP_HPGscaled_wide, "TRIP-combined-HPGscaled.csv")

#=======================================================================================
#Subset sample output tibble for heatmap
gene_list=c("HSPA5","HSP90B1","HYOU1","DNAJB11","DNAJC10","DNAJC3","SDF2L1")
gene_list2=c("ERP29","ERP44","P4HB","PDIA3","PDIA4","PDIA6","TXNDC12","TXNDC5")
gene_list3=c("CCDC47","EMC1","EMC4")
gene_list4=c("CALR","CANX","STT3A","UGGT1")
gene_list5=c("ATL3","CCPG1","CTSB","CTSD","PGRMC1","RTN3")
gene_list6=c("COLGALT1","NID1","P4HA1","PLOD1","PLOD3","PPIA","PPIB","SERPINH1")
gene_list7=c("EDEM3","HSPA8","OS9","SEL1L","VCP")
gene_list8=c("HSPA5","HSP90B1","HYOU1","DNAJB11","DNAJC10")
pd_sample<-TRIP_HPGscaled_sum %>%
  filter(`Gene Symbol`%in% gene_list6) %>%
  arrange(factor(`Gene Symbol`, levels=gene_list)) %>%
  #unite Gene Symbol and experiment column
  unite(col='Gene_exp', c(`Gene Symbol`, experiment), sep='_') 

#Plot heatmap 
ggplot(pd_sample, aes(x = condition, y = factor(`Gene_exp`), fill = mn_scaled_intensity)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 0.9, hjust = 1)) +
  scale_fill_viridis(limits = c(0,1), option="magma", na.value = "grey95") +
  xlab("Sample") +
  ylab("Gene") +
  labs(fill = "Relative Enrichment vs No Hpg") +
  ggtitle("WT") + 
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank()
  )

#=======================================================================================
#plotting heatmaps for individual experiments
pd_sample<-TRIP_HPGscaled_sum %>%
  filter(`Gene Symbol`%in% gene_list8) %>%
  filter(experiment=='C1264R_3') %>%
  arrange(factor(`Gene Symbol`, levels=gene_list8))

#Plot heatmap 
ggplot(pd_sample, aes(x = condition, y = factor(`Gene Symbol`), fill = mn_scaled_intensity)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 0.9, hjust = 1)) +
  scale_fill_viridis(limits = c(0,1), option="magma", na.value = "grey95") +
  xlab("Sample") +
  ylab("Gene") +
  labs(fill = "Relative Enrichment vs No Hpg") +
  ggtitle("WT") + 
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank()
  )

#=======================================================================================
#Subset individual protein sample output in tibble for time-course x-y plots
gene_highlight="ERP44"
TRIP_individual<-TRIP_HPGscaled_sum %>%
  filter(`Gene Symbol` == gene_highlight)
#conversion table to float time points
TC_conversion=tibble(condition = c(" 0_0Hr"," 0_5Hr", " 1_0Hr", " 1_5Hr"," 2_0Hr", " 3_0Hr"),
                     time=c(0.0,0.5,1.0,1.5,2.0,3.0))
#merge time points as floats into tibble
TRIP_individual2<-merge(x=TRIP_individual, y=TC_conversion, by=("condition"))

#plot individual time course
ggplot(TRIP_individual2, aes(x = time, y = mn_scaled_intensity, color = experiment)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = mn_scaled_intensity - mn_scaled_sem, ymax = mn_scaled_intensity + mn_scaled_sem, fill = experiment), alpha = 0.3, linetype = 0) +
  scale_color_manual(values = c("red", "darkgreen", "grey40")) +
  scale_fill_manual(values = c("red", "darkgreen", "grey40")) +
  labs(x = "Time", y = "Scaled Intensity", color = "Experiment") +
  ylim(0, 1.0)+
  theme_classic() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16)) +
  ggtitle(gene_highlight) 

#=======================================================================================
#Subset list of protein output (e.g. PDIs) in tibble for time-course x-y plots
gene_list2=c("ERP29","ERP44","P4HB","PDIA3","PDIA4","PDIA6","TXNDC12","TXNDC5")
treatment<-c("WT")
TRIP_prot_list<-TRIP_HPGscaled_sum %>%
  filter(experiment %in% treatment) %>%
  filter(`Gene Symbol`%in% gene_list2) %>%
  arrange(factor(`Gene Symbol`, levels=gene_list))
#conversion table to float time points
TC_conversion=tibble(condition = c(" 0_0Hr"," 0_5Hr", " 1_0Hr", " 1_5Hr"," 2_0Hr", " 3_0Hr"),
                     time=c(0.0,0.5,1.0,1.5,2.0,3.0))
#merge time points as floats into tibble
TRIP_prot_list2<-merge(x=TRIP_prot_list, y=TC_conversion, by=("condition"))

#plot time course for list of proteins
ggplot(TRIP_prot_list2, aes(x = time, y = mn_scaled_intensity, color = `Gene Symbol`)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = mn_scaled_intensity - mn_scaled_sem, ymax = mn_scaled_intensity + mn_scaled_sem, fill = `Gene Symbol`), alpha = 0.1, linetype = 0) +
  labs(x = "Hours", y = "Relative Enrichment", color = "Protein") +
  ylim(0, 1.0)+
  theme_classic() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16)) +
  ggtitle(treatment) 

#=======================================================================================
#Subset aggregate pathways (e.g. Chaperones, PDIs) in tibble for time-course x-y plots
#select gene_list from pathway annotations
pathway_focus<-pathway_level[10]
gene_list<-pathway_annot_long %>% filter(Pathway==pathway_focus) %>% pull(Protein)
#subset data by pathway and summarize  
TRIP_pathway<-TRIP_HPGscaled_sum %>%
  filter(`Gene Symbol`%in% gene_list) %>%
  group_by(experiment,condition) %>%
  #summarize to median, 1st, and 3rd quartile
  summarize(med_pathway=median(mn_scaled_intensity,na.rm=TRUE),
            qrt1_pathway=quantile(mn_scaled_intensity,0.25,na.rm=TRUE),
            qrt3_pathway=quantile(mn_scaled_intensity,0.75,na.rm=TRUE))
#conversion table to float time points
TC_conversion=tibble(condition = c(" 0_0Hr"," 0_5Hr", " 1_0Hr", " 1_5Hr"," 2_0Hr", " 3_0Hr"),
                     time=c(0.0,0.5,1.0,1.5,2.0,3.0))
#merge time points as floats into tibble
TRIP_pathway2<-merge(x=TRIP_pathway, y=TC_conversion, by=("condition"))

#plot time course by pathways
ggplot(TRIP_pathway2, aes(x = time, y = med_pathway, color = experiment)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin =  qrt1_pathway, ymax = qrt3_pathway, fill = experiment), alpha = 0.3, linetype = 0) +
  scale_color_manual(values = c("red", "darkgreen", "grey40")) +
  scale_fill_manual(values = c("red", "darkgreen", "grey40")) +
  labs(x = "Time", y = "Scaled Intensity", color = "Experiment") +
  ylim(0, 1.0)+
  theme_classic() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16)) +
  ggtitle(pathway_focus)

#=======================================================================================
#plotting full heatmaps of all interactor time courses (for supplemental figures)

#read in pathway annotations for interactors
pathway_annot<-read_csv("PathwaySortedInteractors.csv")
#convert to matching table in long format
pathway_annot_long<-pathway_annot %>%
  mutate(rowid = row_number()) %>%
  pivot_longer(cols = -rowid,
               names_to = "Pathway",
               values_to = "Protein") %>%
  filter(Protein != 'NA')

#add pathway descriptors to master data frame
TRIP_HPGscaled_sum <-left_join(
  TRIP_HPGscaled_sum,pathway_annot_long,by=c('Gene Symbol'='Protein'))

#subset datasets based on interactor annotations and experiment
pd_sample<-TRIP_HPGscaled_sum %>%
  filter(`Gene Symbol`%in% pathway_annot_long$Protein)

%>%
  filter(experiment=="WT")

#Order of pathways
pathway_level<-colnames(pathway_annot)
#ordered protein list
ordered_proteins<-pathway_annot_long %>%
  filter(Protein %in% pd_sample$`Gene Symbol`) %>%
  arrange(factor(Pathway,levels=c(pathway_level)),rowid) %>% 
  pull(Protein)

pd_sample$`Gene Symbol`<-factor(pd_sample$`Gene Symbol`,levels=rev(ordered_proteins))
#pd_sample$Pathway<-factor(pd_sample$Pathway,levels=c(seq(length(pathway_level))))

#Plot heatmap 
ggplot(pd_sample, aes(x = condition, y = `Gene Symbol`,fill = mn_scaled_intensity)) +
  geom_tile() +
  #geom_tile(aes(y = `Gene Symbol`, fill = ), width = 0.2, height = 0.8) +
  scale_x_discrete(expand = c(0, 0.1)) +
  scale_y_discrete(expand = c(0.2, 0)) +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 0.9, hjust = 1),
        axis.text.y = element_text(size = 2)) +
  scale_fill_viridis(limits = c(0,1), option="magma", na.value = "grey95") +
  xlab("Sample") +
  ylab("Gene") +
  labs(fill = "Relative Enrichment vs No Hpg") +
  ggtitle("A2234D") + 
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank()
  )

#plot heatmap of st. errors
ggplot(pd_sample, aes(x = condition, y = `Gene Symbol`,fill = mn_scaled_sem)) +
  geom_tile() +
  #geom_tile(aes(y = `Gene Symbol`, fill = ), width = 0.2, height = 0.8) +
  scale_x_discrete(expand = c(0, 0.1)) +
  scale_y_discrete(expand = c(0.2, 0)) +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 0.9, hjust = 1),
        axis.text.y = element_text(size = 2)) +
  scale_fill_viridis(limits = c(0,0.75), na.value = "grey95") +
  xlab("Sample") +
  ylab("Gene") +
  labs(fill = "Relative Enrichment vs No Hpg") +
  ggtitle("WT") + 
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank()
  )

#Plot of std. error distribution
ggplot(pd_sample, aes(x = experiment, y = mn_scaled_sem, fill = condition)) +
  geom_boxplot() +
  labs(title = "Distribution of SEM",
       x = "Variable",
       y = "SEM") +
  #scale_y_log10() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_classic()
  

#=======================================================================================
#conversion of filtered data to wide format for hierarchical clustering
pd_sample.wide<-pd_sample %>% 
  filter(mn_scaled_intensity!="NaN") %>%
  pivot_wider(
    id_cols=c(Accession,`Gene Symbol`),
    names_from = "condition",
    values_from= mn_scaled_intensity) %>%
  unite(col="rowID", c(`Gene Symbol`, Accession), sep='_') %>%
  column_to_rownames("rowID")

#Hierarchical clustering of time course data and ordered heatmap
pd_sample.matrix<-as.matrix(pd_sample.wide)
pd_sample.dist<-dist(pd_sample.matrix, method="euclidean")
pd_sample.hist<-hclust(pd_sample.dist, "ward.D")

#Define the sub-clusters:
mycl <- cutree(pd_sample.hist, h=max(pd_sample.hist$height/9))
#Check how many sub-clusters are generated:
max(mycl)

#Define a color palette to highlight the sub-clusters on the heatmap:
clusterCols <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]

#plot in heatmap
heatmap.2(pd_sample.matrix,
          main="Hierarchical Cluster",
          Rowv=as.dendrogram(pd_sample.hist),
          Colv=NA, dendrogram="row",
          scale="none",
          density.info="none",
          trace="none",
          col= magma(100),
          RowSideColors= myClusterSideBar)

#bind the sub-cluster assignment with the original data matrix:
foo <- cbind(pd_sample.matrix, clusterID=mycl) %>%
  rownames_to_column()

#Write the sorted data to a new csv file:
write.csv(foo[pd_sample.hist$order,],file='clustered_A2234D_2D_wardD.csv')
#This produces a csv file, which can be opened in Excel




