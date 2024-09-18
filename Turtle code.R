#See tutorial https://www.datacamp.com/tutorial/r-data-import-tutorial
#### set working directory ####
setwd("/Users/korri/Desktop/R files") 
#### 
####

## Install packages
BiocManager::install("microbiome")
devtools::install_github("vmikk/metagMisc")

 if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
 }
 BiocManager::install("phyloseq")
install.packages("devtools")

### Load packages
library(phyloseq)
library(microbiome)
library(data.table)
library(vegan)
library(metagMisc)
library(knitr)
library(ggplot2)
library(dplyr)


#### IMPORT DATA ####

###Import meta data as 'data.frame' **83 obs. of  8 variables
read.csv("turt_meta.csv")
Meta<-read.csv("turt_meta.csv")
str(Meta)  #83 obs. of  6 variables ##note: now says 7 var... including X column
write.csv(Meta, file="turt_meta.csv")


###Import raw data as Formal class 'phyloseq'
turt_data<-readRDS("Turtle_500_filter_no_chmiera.rds")
str(turt_data)    #Formal class 'phyloseq' [package "phyloseq"] with 5 slots
#@Data [1:85, 1:930]

###View data ***85 samples and 930 taxa
summarize_phyloseq(turt_data)
sample_variables(turt_data)
nsamples(turt_data)    #85
ntaxa(turt_data)  #930

### Exploring the object ###
# turt_data
# str(turt_data)  #Formal class 'phyloseq' [package "phyloseq"] with 5 slots
# str(turt_data@tax_table)  #Formal class 'taxonomyTable' [package "phyloseq"] with 1 slot
# turt_data@tax_table
# turt_data@tax_table@.Data["Chloroplast", "Order"]
# str(turt_data@tax_table@.Data)
# test <- as.data.frame(turt_data@tax_table@.Data)
# dim(test)  # 930   6
# Orders_in_data.vec <- test$Order
# length(Orders_in_data.vec)  #930
### /END/ Exploring the object ###



#### Extracting ASVs ####

#See tutorial https://micca.readthedocs.io/en/latest/phyloseq.html
#and https://rdrr.io/bioc/phyloseq/man/assign-otu_table.html
# see https://www.rdocumentation.org/packages/phyloseq/versions/1.16.2/topics/otu_table for use of otu_table
#otu_table() is a phyloseq function which extract the OTU table from the phyloseq object
ASVs<-otu_table(t(turt_data), taxa_are_rows = TRUE) #[930 taxa and 85 samples]
str(ASVs) #Formal class 'otu_table' [package "phyloseq"] with 2 slots   #[1:930, 1:85]
head(ASVs)  #OTU Table: [6 taxa and 85 samples]
ASV_sums <- rowSums(ASVs)
# ?rowSums
hist(ASV_sums, breaks= 200)
length(which(ASV_sums <= 500)) #how many ASVs are there with N reads or fewer? (N=50 in example)   #1 ASV



#### Filter out the Chloroplast #####

#Note: Have to use a filter on the Class of chloroplast rather than the Order, 
#    because NAs exist in the Order column, and these are also dropped when using the subset command by default
Turt_data_no_chloro <-   subset_taxa(turt_data, Class != "Cyanobacteriia") 
Turt_data_no_chloro  #phyloseq-class experiment-level object. OTU Table: [ 929 taxa and 85 samples ]

# Checking number of chpl ASVs 
#Turt_data_only_chloro <- subset_taxa(turt_data, Order =="Chloroplast")

##Filter out mitochondria   ***ISSUE with code,subset_taxa function removes ASVs with Family ID or NA in addition 
#   Mitochondria as specified ... need to remove mito by accessing Class, Phylum, or Kindom that Mito belongs to
# Note: Don't see any mito in taxonomy file ... skip step?
# Turt_data_no_chloro_no_mit<-subset_taxa(Turt_data_no_chloro,Family !="Mitochondria")
# Turt_data_no_chloro_no_mit



#### RAREFY DATA ####

#See rarefaction info on pg.19 at https://www.evolution.unibas.ch/walser/bacteria_community_analysis/2015-02-10_MBM_tutorial_combined.pdf
#Rarefying was first recommended for microbiome counts in order to moderate differences in the presence of rare OTUs.
#Goals fo rarefaction: 1)  Standardize unequal sequencing effort. 2)Enable similarity comparisons along a
#range of samples or a gradient. 3) Enable comparison of different runs or replicates.
#Rarefaction curves represent the diversity as a function of sequencing depth.
#in ecology, rarefaction is repeated sampling procedure to assess species richness
library(vegan)
#see https://rdrr.io/cran/vegan/man/rarefy.html for use of rarecurve. Rarefaction of species richness.
rare<- rarecurve(t(otu_table(Turt_data_no_chloro)), step=50, cex=0.5) #not working...
# ?rarecurve   #Rarefied species richness for community ecologists.

# plot_bar(rare, fill="Rank2")

#see tutorial https://rdrr.io/bioc/phyloseq/man/rarefy_even_depth.html for use of rarefy_even_depth. 
#Resample and OTU table so all samples have same library sizes.
# rarefy_even_depth(physeq, sample.size = min(sample_sums(physeq)),rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
# Arguments phyloseq =(Required). 
#A phyloseq-class object that you want to trim/filter.
#Sample.size =(Optional) a single integer value equal to the number of reads being stimulated AKA depth, and equal to each value returned by sample_sums on the output. 

# Exploring the object
str(Turt_data_no_chloro)  #Formal class 'phyloseq' [package "phyloseq"] with 5 slots
str(Turt_data_no_chloro@otu_table)
str(Turt_data_no_chloro@otu_table@.Data)
Turt_data_no_chloro@otu_table@.Data[1:5,1:5]
rowSums(Turt_data_no_chloro@otu_table@.Data)

Rare_3000 <- rarefy_even_depth(Turt_data_no_chloro, sample.size=3000) # if reads per sample are fewer than sample.size, the sample will be dropped
?set.seed
rowSums(Rare_3000@otu_table@.Data)
summarize_phyloseq(Rare_3000) #Min. number of reads = 3000  #Max. number of reads = 3000
#Total number of reads = 249000
# ?rarefy_even_depth
#2 samples removed because they contained fewer reads than `sample.size`. Up to first five removed samples are: 
# Sample-30
# Sample-33
# Note might need to change rarefication to get core? see alphaT <- microbiome::alpha on line 187
Rare_3000    #[ 929 taxa and 83 samples ]

rare<- rarecurve((otu_table(Rare_3000)), step=50, cex=0.5) #error 
rare
str(rare)    #List of 83

saveRDS(Rare_3000, file ="Turtle_data_Rare3k.rds")



#### Make a taxonomy file ####

Taxa<-tax_table(Rare_3000)
Taxa
write.csv(Taxa, file="Taxa.csv")
str(Taxa)   #Formal class 'taxonomyTable' [package "phyloseq"] with 1 slot

#See https://www.programmingr.com/tutorial/rowsums-in-r/ for sure of RowSums. 
#   rowSums is a Function to sum up the row set in large data set
#   x=data.frame
str(turt_data)  ##read in turt_data from phyloseq to data.frame?
ASVsums<-rowSums(turt_data)
dim(turt_data@otu_table)
turt_data@otu_table[1:5, 1:5]
turt_data.mat <- turt_data@otu_table
turt_data.mat <- turt_data.mat@.Data
str(turt_data.mat)
turt_data.mat <- as.matrix(turt_data.mat)
str(turt_data.mat)
turt_data.mat[1:5,1:5]
ASVsums <- rowSums(turt_data.mat)
ASVsums <- colSums(turt_data.mat)
#error in rowsums: ‘x’ must be an array of at least two dimensions. 
# Which occurs when you feed a vector (single dimensional series of values) into a function 
# which expects to look at an array.
Rare_3000@otu_table
turt_data@sam_data
Rare_3000@sam_data@.Data



#### META DATA ####
  ## Meta does not have updated sample types
read.csv(file="Meta.csv")
#see tutorial https://rdrr.io/bioc/phyloseq/man/sample_data-methods.html
# sample_data is a method for both constructing and accessing a table of sample-level variables (sample_data-class), 
# which in the phyloseq-package is represented as a special extension of the data.frame-class. 
# Note: use str() function to check data format
# Meta1<-sample_data(turt_data)   #Note: samples that were rem. by rarefaction, Sample-30 was soil, Sample-33 was egg.
Meta<-sample_data(Rare_3000) 
Meta_data<-meta(Rare_3000)
head(Meta_data)
head(Meta) #6 samples by 5 sample variables
str(Meta)   #'data.frame':	83 obs. of  5 variables
#   formal class 'sample_data' [package "phyloseq"] with 4 slots
write.csv(Meta, file="Turt_meta.csv")

###raw data analysis file:
meta<-import_qiime_sample_data("Meta_turtle.txt")

head(meta)
Rare_turtle <- phyloseq(otu_table(Rare_3000, taxa_are_rows=FALSE), tax_table(Rare_3000))
ps1 <-merge_phyloseq(ps,meta)

#The summarize_phyloseq function will give information on whether data is compositional (i.e. containing only relative rather than absolute data/quantitative description of parts of some whole) or not, reads (min. max, median, average), sparsity, presence of singletons and sample variables.
summarize_phyloseq((turt_data)) # see basic info
summarize_phyloseq(Rare_3000)


##### Summarizing amplicon sequence variants (ASVs) #####

ASVs<-otu_table(t(Rare_3000), taxa_are_rows = TRUE) # build an OTU table
ASVs
dim(ASVs) #dim() function retrieves or sets the dimension of an object.  ##929  83
ASVs[1:5,1:5] # view ASVs data

#the rowSums() function returns the sums of each row in the data set. 
ASV_sums <- rowSums(ASVs) # Total reads per ASV
ASV_sums

hist(ASV_sums, breaks= 200) # Full histogram
hist(ASV_sums[ASV_sums < 50000], breaks= 20) # Remove outlier ASV_sums
min(ASV_sums)   #8
which(ASV_sums > 200000) #named integer(0)
length(which(ASV_sums <= 1000)) # how many ASVs are there with N reads or fewer? (N=1000 in example) #883



#### ALPHA DIVERSITY ####

library(microbiome)
library(knitr)
prevalence(Rare_3000)
alphaT <- microbiome::alpha(Rare_3000, index =  "all")   # Warning: Core is empty with the given abundance and prevalence 
# The function core_abundance(), which is run as one of the indices run within the alpha() function, when run at the settings internal to alpha
#   results in no output "Core is empty with the given abundance and prevalence"
#   it is possible to run on its own, if need to change the parameters, e.g., 
#core_abundance(x = Rare_3000)
str(Rare_3000) #Formal class 'phyloseq' [package "phyloseq"] with 5 slots
#Detection prevalence is the number of predicted positive events (both true positive and false positive) divided by the total number of predictions.
#core_abundance (x,...) function Calculates the community core abundance index.
#core_abundance(Rare_3000, detection = 0.2, prevalence = /100)
dim(alphaT)  #83 22
colnames(alphaT)
alphaT[1:5,1:5]
str(alphaT)
colnames(alphaT)

write.csv(alphaT, file="turtle_div.csv") 
turt_div <- read.csv("turtle_div.csv")
head(turt_div)
rownames(turt_div) <- turt_div$X
str(turt_div)   #'data.frame':	83 obs. of  23 variables:

# alpha_d <- read.csv(file= "turtle_div.csv")
alpha_d <- turt_div
head(alpha_d)
str(alpha_d)


#### MERGE ALPHA DIVERSITY AND META DATA ####

# Import the metadata
Meta_n <- read.csv(file="Meta.csv")
head(Meta_n)
nrow(Meta_n)  #83
nrow(alpha_d)  #83
union(Meta_n$X, alpha_d$X) # checking to see if the sample is in both dataframes
length(union(Meta_n$X, alpha_d$X)) # checking to see if the sample is in both dataframes #83

colnames(Meta_n)
colnames(alpha_d)

# Merge meta n and alpha d
Combined <- merge(x = Meta_n, y = alpha_d, by = "X")
write.csv(Combined, file = "combinedTURT.csv")

Combined<-read.csv("combinedTURT.csv")


Combined<-merge(Meta_n,alpha_d, by = "X")
write.csv(Combined, file = "combinedTURT.csv")
Combined <- read.csv("combinedTURT.csv")
str(Combined)  #data.frame: 83 obs. of  29 variable

####Subset phyloseq object to the OTUs that are present in all samples. 
# turt_data <- read.csv("turt_data.csv")
# turt_data
str(turt_data)
# head(turt_data)
str(Combined)
View(Combined)
library(devtools)

library(metagMisc)

#install.packages("remotes")
#remotes::install_github("vmikk/metagMisc")

#Extract OTUs present in all samples
library(metagMisc)


#turt_data
#str(turt_data)  #Formal class 'phyloseq' [package "phyloseq"] with 5 slots
#turt_data@otu_table[1:5,1:5]



#### PLOT COMPOSITION ####

#eg:  ps2 <- phyloseq_extract_shared_otus(turt_data, samp_names = c("B", "C"))
##  phyloseq_extract_shared_otus(turt_data, samp_names = sample_names(turt_data))
?phyloseq_extract_shared_otus
#use phyloseq_extract_shared_otus function to subset phyloseq object to the OTUs that are present in all samples. 

###Cloaca vs egg vs soil   ###should rarefied data rather than turt_data be used here?
str(turt_data)  #Formal class 'phyloseq'
str(Combined)   #'data.frame'
str(Rare_3000)  #Formal class 'phyloseq'
ps.c.e.s <- phyloseq_extract_shared_otus(Rare_3000, samp_names = c("Sample-104", "Sample-113", "Sample-28", "Sample-32", "Sample-46", "Sample-53", "Sample-62", "Sample-64", "Sample-83", "Sample-116", "Sample-121", "Sample-31", "Sample-51", "Sample-59", "Sample-65", "Sample-70", "Sample-94", "Sample-103", "Sample-123", "Sample-125", "Sample-54")) 
ps.c.e.s
otu_table(ps.c.e.s)    #[9 taxa and 21 samples]
#use aggregate_rare function to combine rare taxa ##at phylum level
ps.c.e.s_phylum<-aggregate_rare(ps.c.e.s, level = "Phylum", detection = 1/100, prevalence = 50/100) 
ps.c.e.s_phylum   #[ 5 taxa and 21 samples ]
Comp.c.e.s_phylum<- plot_composition(ps.c.e.s_phylum, group_by = "Type1") 
Comp.c.e.s_phylum
##at genus level
ps.c.e.s_genus<-aggregate_rare(ps.c.e.s, level = "Genus", detection = 1/100, prevalence = 50/100) 
ps.c.e.s_genus
Comp.c.e.s_genus<- plot_composition(ps.c.e.s_genus, group_by = "Type1") 
Comp.c.e.s_genus


### Egg vs colaca  ##at phylum level
ps.e.c_phylum <- phyloseq_extract_shared_otus(Rare_3000, samp_names = c("Sample-116", "Sample-121", "Sample-31", "Sample-51", "Sample-59", "Sample-65", "Sample-70", "Sample-94", "Sample-104", "Sample-113", "Sample-28", "Sample-32", "Sample-46", "Sample-53", "Sample-62", "Sample-64", "Sample-83"))
ps.e.c_phylum<-aggregate_rare(ps.e.c_phylum, level = "Phylum", detection = 1/100, prevalence = 50/100)
Comp.e.c_phylum <- plot_composition(ps.e.c_phylum, group_by = "Type1") + theme_classic() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ theme(legend.title = element_text(size=16), legend.text = element_text(size=14), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16))
Comp.e.c_phylum
##at genus level
ps.e.c_genus <- phyloseq_extract_shared_otus(Rare_3000, samp_names = c("Sample-116", "Sample-121", "Sample-31", "Sample-51", "Sample-59", "Sample-65", "Sample-70", "Sample-94", "Sample-104", "Sample-113", "Sample-28", "Sample-32", "Sample-46", "Sample-53", "Sample-62", "Sample-64", "Sample-83"))
ps.e.c_genus<-aggregate_rare(ps.e.c_genus, level = "Genus", detection = 1/100, prevalence = 50/100)
Comp.e.c_genus <- plot_composition(ps.e.c_genus, group_by = "Type1") + theme_classic() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ theme(legend.title = element_text(size=16), legend.text = element_text(size=14), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16))
Comp.e.c_genus

#Combined$Type1 = factor(combined$Type1, levels = c("egg, cloaca"))


### Egg vs Soil
##at phylum level
ps.e.s_phylum <- phyloseq_extract_shared_otus(Rare_3000, samp_names = c("Sample-116", "Sample-121", "Sample-31", "Sample-51", "Sample-59", "Sample-65", "Sample-70", "Sample-94", "Sample-103", "Sample-123", "Sample-125", "Sample-54"))
ps.e.s_phylum 
ps.e.s_phylum<-aggregate_rare(ps.e.s_phylum, level = "Phylum", detection = 1/100, prevalence = 50/100)
ps.e.s_phylum
Comp.e.s_phylum<- plot_composition(ps.e.s_phylum, group_by = "Type1") + theme_classic() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ theme(legend.title = element_text(size=16), legend.text = element_text(size=14), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16))
Comp.e.s_phylum

##at genus level
ps.e.s_genus <- phyloseq_extract_shared_otus(Rare_3000, samp_names = c("Sample-116", "Sample-121", "Sample-31", "Sample-51", "Sample-59", "Sample-65", "Sample-70", "Sample-94", "Sample-103", "Sample-123", "Sample-125", "Sample-54"))
ps.e.s_genus 
ps.e.s_genus<-aggregate_rare(ps.e.s_genus, level = "Phylum", detection = 1/100, prevalence = 50/100)
ps.e.s_genus
Comp.e.s_genus<- plot_composition(ps.e.s_genus, group_by = "Type1") + theme_classic() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ theme(legend.title = element_text(size=16), legend.text = element_text(size=14), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16))
Comp.e.s_genus


####Subsetting only one sample type

# New_object_name_cloaca<- subset_samples(Rare_3000, Type1 == "cloaca")
# New_object_name_egg <- subset_samples(Rare_3000, Type1 == "egg")
# New_object_name_soil <- subset_samples(Rare_3000, Type1 == "soil")
# New_object_name_final_turt <- subset_samples(Rare_3000, Type1 == "final turt" )


#### Hatchling vs Soil
ps.h.s<-phyloseq_extract_shared_otus(Rare_3000, samp_names = c("Sample-61", "Sample-57", "Sample-35", "Sample-3", "Sample-26", "Sample-124", "Sample-115", "Sample-105", "Sample-102", "Sample-103", "Sample-123", "Sample-125", "Sample-54"))
ps.h.s
ps.h.s_phylum<-aggregate_rare(ps.h.s, level = "Phylum", detection = 1/100, prevalence = 50/100)
Comp.h.s<- plot_composition(ps.h.s_phylum, group_by = "Type1") + theme_classic() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ theme(legend.title = element_text(size=16), legend.text = element_text(size=14), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16))
Comp.h.s

#### Hatchling vs cloaca
ps.h.c<-phyloseq_extract_shared_otus(Rare_3000, samp_names = c("Sample-61", "Sample-57", "Sample-35", "Sample-3", "Sample-26", "Sample-124", "Sample-115", "Sample-105", "Sample-102","Sample-104", "Sample-113", "Sample-28", "Sample-32", "Sample-46", "Sample-53", "Sample-62", "Sample-64", "Sample-83"))
ps.h.c
ps.h.c_phylum<-aggregate_rare(ps.h.c, level = "Phylum", detection = 1/100, prevalence = 50/100)
Comp.h.c<- plot_composition(ps.h.c, group_by = "Type1") + theme_classic() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ theme(legend.title = element_text(size=16), legend.text = element_text(size=14), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16))
Comp.h.c

#### Hatchling vs final turt
ps.h.t<-phyloseq_extract_shared_otus(Rare_3000, samp_names = c("Sample-61", "Sample-57", "Sample-35", "Sample-3", "Sample-26", "Sample-124", "Sample-115", "Sample-105", "Sample-102","Sample-27", "Sample-52", "Sample-58", "Sample-91", "Sample-95"))
ps.h.t
ps.h.t_phylum<-aggregate_rare(ps.h.t, level = "Phylum", detection = 1/100, prevalence = 50/100)
Comp.h.t<- plot_composition(ps.h.t_phylum, group_by = "Type1") + theme_classic() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ theme(legend.title = element_text(size=16), legend.text = element_text(size=14), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16))
Comp.h.t


#Error in x[seq_len(n)] : object of type 'S4' is not subsettable 
#.... Rare_3000 is a phyloseq objusct, might need a data.frame



#### SUBSAMPLE DATA BY SAMPLE TYPE ####

###Sub-sample only cloaca, egg, and soil sample types
Sub_samp<-subset_samples(Rare_3000, Type1 %in% c("cloaca", "egg", "soil"))
sub_samp_ag <- aggregate_rare(Sub_samp, level = "Phylum", detection = 1/100, prevalence = 50/100)

###Turtle ONLY
Turt_sub_samp <-subset_samples(Rare_3000, Type1 %in% c("hatchling", "One_year_wild_cloaca", "One_year_wild_body", "final_turt"))
sub_samp_Turt_Only <- aggregate_rare(Turt_sub_samp, level = "Phylum", detection = 1/100, prevalence = 50/100)

###Host ONLY
Host_sub_samp <-subset_samples(Rare_3000, Type1 %in% c("hatchling", "One_year_wild_cloaca", "One_year_wild_body","final_turt", "egg", "mid egg", "final egg", "in egg"))
sub_samp_Host_Only <- aggregate_rare(Host_sub_samp, level = "Phylum", detection = 1/100, prevalence = 50/100)


#### all categories
All_samp<-subset_samples(Rare_3000, Type1 %in% c("cloaca", "egg", "in egg", "mid egg", "final egg", "hatchling", "final turt", "One_year_wild_cloaca", "One_year_wild_body", "soil", "mid sub", "final sub"))
All_samp_ag <- aggregate_rare(All_samp, level = "Phylum", detection = 1/100, prevalence = 50/100)
All_samp_ag   #[ 9 taxa and 83 samples ]


#### Naming ASVs #### 
ASV_phyla <- microbiome::transform(Rare_3000, "compositional")
ASV_phyla <- aggregate_rare(ASV_phyla, level = "Phylum", detection = 1/100, prevalence = 30/100)
#The psmelt function changes the data from a phyloseq format to regular format - this makes it easier to make plots 
ASV_phyla<- psmelt(ASV_phyla)
str(ASV_phyla)  #'data.frame':	581 obs. of  10 variables:
ASV_phyla

 ##note: Rare_3000_data does not have updated meta

#### Make plot with ggplot with better ASV names. you can also use any of your meta data to collapse samples into groups. 
# head(Combined)
Rare_3000_data <- psmelt(Rare_3000)
Combined$Type1 = factor(Combined$Type1, levels = c("cloaca", "egg", "in egg", "mid egg", "final egg", "hatchling", "final turt", "One_year_wild_cloaca", "One_year_wild_body", "soil", "mid sub", "final sub"))
Best_plot<- ggplot(Rare_3000_data, aes(x=Type1, y = Abundance, fill= Phylum)) +geom_bar(stat="summary") + theme_classic() + xlab("Sample Category") + theme(legend.title = element_text(size=16), legend.text = element_text(size=14)) + theme(axis.title.x = element_text(size=16), axis.text.x = element_text(size=9, face="bold"), axis.title.y = element_text(size=16))
Best_plot



plot_composition(Rare_3000, group_by = "Type1") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())

Comp_plot <- plot_composition(Rare_3000, group_by = "Type1", color = "SampleID") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank() + xlab("Sample")) 
Comp_plot

###order categories
Combined$Type1 = factor(Combined$Type1, levels = c("cloaca", "egg", "in egg", "mid egg", "final egg", "hatchling", "final turt", "One_year_wild_cloaca", "One_year_wild_body", "soil", "mid sub", "final sub"))



#### Try new package 
#devtools::install_github("david-barnett/microViz")
library("microViz")
library(patchwork)

comp_barplot(Rare_3000,
             tax_level = "Phylum", n_taxa = 8,
             bar_outline_colour = NA,
             sample_order = "bray",
             bar_width = 0.7,
             palette = distinct_palette(8, pal = "kelly"),
) + coord_flip() %>%
  patchwork::wrap_plots(nrow = 3, guides = "collect")
#NAs detected in phyloseq tax_table: Consider using tax_fix() to make taxa uniquely identifiable
#Error: Can only wrap ggplot and/or grob objects or a list of them
Rare_3000 %>% tax_fix_interactive()




##### Beta diversity. BRAY AND JACCARD DISTANCE #####

### Ordination Plots

MDS1 <- ordinate(sub_samp_ag, method = "MDS", distance = "bray", weighted = TRUE)

MDS_Bray <- plot_ordination(sub_samp_ag, MDS1, color = "Type1", shape = "Sample_cat", title = "Bray-Curtis Distance")  + stat_ellipse(type = "norm", linetype = 1) + theme_classic()+ guides(size = "none") + geom_point(size=4.5)+ theme(legend.title = element_text(size=18), legend.text = element_text(size=16)) + theme(axis.title.x = element_text(size=18), axis.title.y = element_text(size=18))

MDS1 <- ordinate(Rare_3000, method = "MDS", distance = "bray", weighted = TRUE)
MDS_Bray<- plot_ordination(Rare_3000, MDS1, color = "Type1", shape = "Sample_cat", title = "Bray-Curtis Distance") + stat_ellipse(type = "norm", linetype = 1) + theme_classic() + theme_classic() + geom_point(size=4) + theme(legend.title = element_text(size=16), legend.text = element_text(size=14) + theme(axis.title.x = element_text(size=16), axis.text.x = element_text(size=9), axis.title.y = element_text(size=16)))
MDS_Bray


MDS2<- ordinate(Rare_3000, method = "MDS", distance = "jaccard")
MDS_Jaccard<- plot_ordination(Rare_3000, MDS2, color = "Type1", title = "Jaccard Distance")+ stat_ellipse(type = "norm", linetype = 1) + theme_classic() + theme_classic() + geom_point(size=4) + theme(legend.title = element_text(size=16), legend.text = element_text(size=14) + theme(axis.title.x = element_text(size=16), axis.text.x = element_text(size=9), axis.title.y = element_text(size=16)))
MDS_Jaccard



#### Making figures into panel. This is easy with grid.arrange. 

library(gridExtra)
library(phyloseq)
library(vegan)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
#install.packages("devtools")
#install.packages("cluster")
library(cluster)
library(ggplot2)

grid.arrange(MDS_Bray, MDS_Jaccard)



###Let's do  a test to ask if these samples have the same centroid (mean in n-dimensional space). This is similar to ANOVA, but for distance data. Not only does it look at centroid, but also dispersion (variation from centriod). 
#This is an adonis test- also called a PERMANOVA. Make sure you save the results.

### Beta-diversity Tests

Bray_dist<- phyloseq::distance(Rare_3000, method = "bray", weighted = TRUE)
Sample_star <- data.frame(sample_data(Rare_3000))
Test_bray<-adonis(Bray_dist ~ Type1, data = Sample_star)
Test_bray

pairwise.adonis(Bray_dist,Sample_star$Type1, Sample_star$Type1)


Ja_dist<- phyloseq::distance(Rare_3000, method = "jaccard")
Sample_star <- data.frame(sample_data(Rare_3000))
Test_Ja<-adonis(Ja_dist ~ Type1, data = Sample_star)
Test_Ja

pairwise.adonis(Ja_dist,Sample_star$Type1)



####Bray and Jaccard for only egg, cloaca, soil subsetted data 
Sub_samp_e.c.s<-subset_samples(Rare_3000, Type1 %in% c("cloaca", "egg", "in egg", "mid egg", "final egg", "hatchling", "final turt", "One_year_wild", "soil", "mid sub", "final sub"))
sub_samp_ <- aggregate_rare(Sub_samp, level = "Phylum", detection = 1/100, prevalence = 50/100)

MDS1 <- ordinate(sub_samp_ag, method = "MDS", distance = "bray", weighted = TRUE)
MDS_Bray_sub <- plot_ordination(sub_samp_ag, MDS1, color = "Type1", title = "Bray-Curtis Distance") + theme_classic() + stat_ellipse(type = "norm", linetype = 1) + theme_classic() + theme_classic() + geom_point(size=4) + theme(legend.title = element_text(size=16), legend.text = element_text(size=14) + theme(axis.title.x = element_text(size=16), axis.text.x = element_text(size=9), axis.title.y = element_text(size=16)))
MDS_Bray_sub


MDS2 <- ordinate(sub_samp_ag, method = "MDS", distance = "jaccard", weighted = TRUE)
MDS_Jac_sub <- plot_ordination(sub_samp_ag, MDS2, color = "Type1", title = "Jaccard Distance") + stat_ellipse(type = "norm", linetype = 1) + theme_classic()  + theme_classic() + geom_point(size=4) + theme(legend.title = element_text(size=16), legend.text = element_text(size=14) + theme(axis.title.x = element_text(size=16), axis.text.x = element_text(size=9), axis.title.y = element_text(size=16)))
MDS_Jac_sub


grid.arrange(MDS_Bray_sub, MDS_Jac_sub) 


### Test
Bray_dist_sub<- phyloseq::distance(sub_samp_ag, method = "bray", weighted = TRUE)
Sample_star <- data.frame(sample_data(sub_samp_ag))
Test_bray_sub<-adonis(Bray_dist_sub ~ Type1, data = Sample_star)
Test_bray_sub

pairwise.adonis(Bray_dist_sub,Sample_star$Type1, Sample_star$Type1)


Ja_dist_sub<- phyloseq::distance(sub_samp_ag, method = "jaccard")
Sample_star <- data.frame(sample_data(sub_samp_ag))
Test_Ja_sub<-adonis(Ja_dist_sub ~ Type1, data = Sample_star)
Test_Ja_sub

pairwise.adonis(Ja_dist_sub,Sample_star$Type1)


##### Ordination Plots ####

##Env vs host
MDS1 <- ordinate(sub_samp_ag, method = "MDS", distance = "bray", weighted = TRUE)

MDS_Bray_env_host <- plot_ordination(Rare_3000, MDS1, color = "Sample_cat", title = "Bray-Curtis Distance") + stat_ellipse(type = "norm", linetype = 1) + theme_classic()+ guides(size = "none") + geom_point(size=4) + theme(legend.title = element_text(size=18), legend.text = element_text(size=16)) + theme(axis.title.x = element_text(size=18), axis.title.y = element_text(size=18))

MDS1 <- ordinate(Rare_3000, method = "MDS", distance = "bray", weighted = TRUE)
MDS_Bray_env_host <- plot_ordination(turt_data, MDS1, color = "Sample_cat", title = "Bray-Curtis Distance") + stat_ellipse(type = "norm", linetype = 1) + theme_classic() + theme_classic() + geom_point(size=4) + theme(legend.title = element_text(size=16), legend.text = element_text(size=14) + theme(axis.title.x = element_text(size=16), axis.text.x = element_text(size=9), axis.title.y = element_text(size=16)))
MDS_Bray_env_host

MDS2<- ordinate(Rare_3000, method = "MDS", distance = "jaccard")

MDS_Jaccard_env_host<- plot_ordination(Rare_3000, MDS2, color = "Sample_cat", title = "Jaccard Distance") + stat_ellipse(type = "norm", linetype = 1) + theme_classic()  + theme_classic() + geom_point(size=4) + theme(legend.title = element_text(size=16), legend.text = element_text(size=14) + theme(axis.title.x = element_text(size=16), axis.text.x = element_text(size=9), axis.title.y = element_text(size=16)))
MDS_Jaccard_env_host

grid.arrange(MDS_Bray_env_host, MDS_Jaccard_env_host)



###Tests
Bray_dist_env_host<- phyloseq::distance(Rare_3000, method = "bray", weighted = TRUE)
Sample_star <- data.frame(sample_data(Rare_3000))
Test_bray_env_host<-adonis(Bray_dist ~ Sample_cat, data = Sample_star)
Test_bray_env_host

pairwise.adonis(Bray_dist_env_host,Sample_star$Sample_cat, Sample_star$Type1)


Ja_dist_env_host<- phyloseq::distance(Rare_3000, method = "jaccard")
Sample_star <- data.frame(sample_data(Rare_3000))
Test_Ja_env_host<-adonis(Ja_dist ~ Type1, data = Sample_star)
Test_Ja_env_host

pairwise.adonis(Ja_dist_env_host,Sample_star$Sample_cat)


####Bray and Jaccard for only egg, cloaca, soil subsetted data
#   consider removing outlier cloaca sample --figure out which sample is the outlier
MDS1 <- ordinate(sub_samp_ag, method = "MDS", distance = "bray", weighted = TRUE)
MDS_Bray_sub <- plot_ordination(sub_samp_ag, MDS1, color = "Type1", title = "Bray-Curtis Distance") + stat_ellipse(type = "norm", linetype = 1)+ ylab("NMDS2") + xlab("NMDS1") + theme_classic() + ylab("NMDS2") + xlab("NMDS1") + theme_classic() + geom_point(size=3) + theme(legend.title = element_text(size=16), legend.text = element_text(size=14) + theme(axis.title.x = element_text(size=16), axis.text.x = element_text(size=9), axis.title.y = element_text(size=16)))
MDS_Bray_sub


MDS2 <- ordinate(sub_samp_ag, method = "MDS", distance = "jaccard", weighted = TRUE)
MDS_Jac_sub <- plot_ordination(sub_samp_ag, MDS2, color = "Type1", title = "Jaccard Distance")+ stat_ellipse(type = "norm", linetype = 1) + ylab("NMDS2") + xlab("NMDS1") + theme_classic() + ylab("NMDS2") + xlab("NMDS1") + theme_classic() + geom_point(size=3) + theme(legend.title = element_text(size=16), legend.text = element_text(size=14) + theme(axis.title.x = element_text(size=16), axis.text.x = element_text(size=9), axis.title.y = element_text(size=16)))
MDS_Jac_sub

grid.arrange(MDS_Bray_sub, MDS_Jac_sub) 

### Test
Bray_dist_sub<- phyloseq::distance(sub_samp_ag, method = "bray", weighted = TRUE)
Sample_star <- data.frame(sample_data(sub_samp_ag))
Test_bray_sub<-adonis(Bray_dist_sub ~ Type1, data = Sample_star)
Test_bray_sub

pairwise.adonis(Bray_dist_sub,Sample_star$Type1, Sample_star$Type1)


Ja_dist_sub<- phyloseq::distance(sub_samp_ag, method = "jaccard")
Sample_star <- data.frame(sample_data(sub_samp_ag))
Test_Ja_sub<-adonis(Ja_dist_sub ~ Type1, data = Sample_star)
Test_Ja_sub

pairwise.adonis(Ja_dist_sub,Sample_star$Type1)



###host through development:

sub_samp_host<-subset_samples(Rare_3000, Type1 %in% c("egg", "mid egg", "final egg", "hatchling", "final turt", "One_year_wild_cloaca", "One_year_wild_cloaca"))

MDS1 <- ordinate(sub_samp_host, method = "MDS", distance = "bray", weighted = TRUE)
MDS_Bray_host<- plot_ordination(sub_samp_host, MDS1, color = "Type1", title = "Bray-Curtis Distance") + stat_ellipse(type = "norm", linetype = 1)+ geom_point(size=3) + theme_classic()+ guides(size = "none") + theme(legend.title = element_text(size=18), legend.text = element_text(size=16)) + theme(axis.title.x = element_text(size=18), axis.title.y = element_text(size=18))
MDS_Bray_host

######## subsetting for turtles

sub_samp_turtles<-subset_samples(Rare_3000, Type1 %in% c("hatchling", "final turt", "One_year_wild_cloaca", "One_year_wild_body"))

MDS1 <- ordinate(sub_samp_turtles, method = "MDS", distance = "bray", weighted = TRUE)
MDS_Bray_turtles<- plot_ordination(sub_samp_turtles, MDS1, color = "Type1", title = "Bray-Curtis Distance") + theme_classic() + stat_ellipse(type = "norm", linetype = 1) + geom_point(size=3) + guides(size = "none") + theme(legend.title = element_text(size=18), legend.text = element_text(size=16)) + theme(axis.title.x = element_text(size=18), axis.title.y = element_text(size=18))
MDS_Bray_turtles


######## subsetting for eggs and environment through time
sub_samp_egg.env<-subset_samples(Rare_3000, Type1 %in% c("egg", "soil", "mid egg", "mid sub", "final egg", "final sub"))

MDS2 <- ordinate(sub_samp_egg.env, method = "MDS", distance = "bray", weighted = TRUE)
MDS_Bray_egg.env<- plot_ordination(sub_samp_egg.env, MDS2, color = "Type1", title = "Bray-Curtis Distance") + theme_classic() + stat_ellipse(type = "norm", linetype = 1) + geom_point(size=3) + guides(size = "none") + theme(legend.title = element_text(size=18), legend.text = element_text(size=16)) + theme(axis.title.x = element_text(size=18), axis.title.y = element_text(size=18))
MDS_Bray_egg.env



##### Graph design #### 

+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())

+ theme_classic() + geom_point(size=3) + theme(legend.title = element_text(size=16), legend.text = element_text(size=14) + theme(axis.title.x = element_text(size=16), axis.text.x = element_text(size=9), axis.title.y = element_text(size=16)))

#add elipses
+ stat_ellipse(type = "norm", linetype = 1)


###order categories
Combined$Type1 = factor(Combined$Type1, levels = c("cloaca", "egg", "in egg", "mid egg", "final egg", "hatchling", "final turt", "One_year_wild", "soil", "mid sub", "final sub"))






##### VENN DIAGRAM ####

if (!require(devtools)) install.packages("devtools")
devtools::install_github("gaospecial/ggVennDiagram")
install.packages("ggvenn")
install.packages("ggVennDiagram")

library("ggVennDiagram")
library("ggvenn")


install.packages("VennDiagram")
#load Venn diagram package
library("VennDiagram")


# move to new plotting page
#grid.newpage()

# create pairwise Venn diagram
#draw.pairwise.venn(area1=20, area2=45,cross.area=10, category=c("soil", "egg"),fill=c("Red","Yellow"))
# 3 way
#draw.triple.venn(area1=40, area2=15, area3=10, n12=5, n23=12, n13=4, n123=2, category=c("soil","egg","cloaca"),col="Red",fill=c("Green","Yellow","Blue"))

# use data frame as input
#M <-tibble(value=c(1,3,2,7,5),'cloaca'=c(TRUE,FALSE,TRUE,FALSE,FALSE),'soil'=c(TRUE,TRUE,FALSE,FALSE,TRUE), 'egg'=c(TRUE,FALSE,FALSE,TRUE,TRUE))

# create Venn diagram and display all sets
#ggvenn(M)

# use list as input
#D <-list('egg'=c("Sample-116", "Sample-121", "Sample-31", "Sample-33", "Sample-51", "Sample-59", "Sample-65", "Sample-70", "Sample-94"),'cloaca'=c("Sample-104", "Sample-113", "Sample-28", "Sample-32", "Sample-46", "Sample-53", "Sample-62", "Sample-64", "Sample-83"),'soil'=c("Sample-103", "Sample-123", "Sample-125", "Sample-54"))

# creating venn diagram for four sets
# and displaying only two sets
#ggvenn(D,c("egg","cloaca", "soil"),show_percentage=FALSE, fill_color=c("blue","yellow", "green"))

#get_vennlist(Combined, ...)

## S4 method for signature 'phyloseq'
#get_vennlist(obj, factorNames, ...)

## S4 method for signature 'data.frame'
#get_vennlist(o)

####New code
#get_vennlist(combined, ...)

## S4 method for signature 'phyloseq'
#get_vennlist(obj, factorNames, ...)

## S4 method for signature 'data.frame'
#get_vennlist(combined, sampleinfo = NULL, factorNames = NULL, ...)


####
#data(turt_data)
#vennlist <- get_vennlist(turt_data, factorNames="group")
#vennlist
#library(VennDiagram)
#venn.diagram(vennlist, height=5, width=5, filename = "./test_venn.pdf", alpha = 0.85, fontfamily = "serif", fontface = "bold",cex = 1.2, cat.cex = 1.2, cat.default.pos = "outer", cat.dist = c(0.22,0.22,0.12,0.12), margin = 0.1, lwd = 3, lty ='dotted', imagetype = "pdf"




##### use this code for venn diagram
library(devtools)
#nstall.packages("Git")
#install_github("Russel88/MicEco")
library("MicEco")
#ps_venn(Rare_3000, group="Type1")

#usage:
ps_venn(
  ps,
  group,
  fraction = 0,
  weight = FALSE,
  relative = TRUE,
  plot = TRUE,
  ...
)

ps_venn(Sub_samp, "Type1", font = 14, labels = list(cex = 1.3), fraction = 0.1)



##### ALPHA DIVERSITY ####

#We need the rarefied data, use "Rare_3000"
#Open combined file
Combined<-read.csv("combinedTURT.csv")

## PLOTS
library(phyloseq)
Sub_samp<-subset_samples(Rare_3000, Type1 %in% c("cloaca", "egg", "soil"))
sub_samp_ag <- aggregate_rare(Sub_samp, level = "Phylum", detection = 1/100, prevalence = 50/100)

subset_samples(physeq = Rare_3000, )


### There are many alpha diversity measurements. Lets' focus on a  few: Observed *richness) Shannon, and evenness_simpson and make some plots. 

#Observed: How many "species" are observed. 

#Shannon: How difficult it is to predict the identity of a randomly chosen individual.Takes account of richness and evenness. 
#Simpson: The probability that two randomly chosen individuals are the same species. Only looks at evenness. 
library(vegan)
library(ggplot2)
library(gridExtra)
diversity(Combined, shannon)    #Error in diversity(Combined, shannon) : input data must be numeric

#Colour code for host vs env, reorder sequence by time in development/sample point
Combined$Type1 = factor(Combined$Type1, levels = c("soil","cloaca","egg","mid sub","mid egg","final sub","final egg","in egg", "hatchling","final turt","One_year_wild")) 
Richness_shannon<-ggplot(Combined) +geom_boxplot(aes(x=Type1, y =diversity_shannon, color=Sample_cat)) + geom_jitter(aes(x=Type1, y =diversity_shannon)) + theme_classic() + ylab("Richness") + xlab("Sample Type") + theme(axis.title.x = element_text(size=16), axis.text.x = element_text(size=11), axis.title.y = element_text(size=16))
Richness_shannon  + scale_x_discrete(limits = c("soil","cloaca","egg","mid sub","mid egg","final sub","final egg","in egg", "hatchling","final turt","One_year_wild")) 


Richness_simpson<-ggplot(Combined) +geom_boxplot(aes(x=Type1, y =evenness_simpson)) + geom_jitter(aes(x=Type1, y =evenness_simpson)) + theme_classic() + ylab("Richness") + xlab("Sample Type") + theme(axis.title.x = element_text(size=16), axis.text.x = element_text(size=11), axis.title.y = element_text(size=16))
Richness_simpson + scale_x_discrete(limits = c("soil","cloaca","egg","mid sub","mid egg","final sub","final egg","in egg", "hatchling","final turt","One_year_wild")) 

grid.arrange(Richness_shannon, Richness_simpson)

## TESTS
install.packages("FSA")
library(FSA)

dunnTest(observed~Type1, data=Combined)    #Warning message:  Type1 was coerced to a factor.

#Stats for richness 
kruskal.test(observed~Type1, data= Combined)

dunnTest(observed~Type1, data= Combined)    #Warning message:  Type1 was coerced to a factor.






######### PLAN FOR FIGIRES IN PAPER ###########

#   Fig.1 Alpha and beta diversities ####Note: use pd whole tree for alpha and unifrac for beta div
grid.arrange(Richness_shannon, wt.unifrac)
#   Fig.2 Core bubble plot
#   Fig.3 NTI (Nearest taxon index) ##see Van Veelen et al., 2017 for eg.
#   Fig.4 Source and Venn diagram


#### FIGURE 1. Alpha Diversity ####  

# Add shannon div to figure1
#Shannon diversity: How difficult it is to predict the identity of a randomly chosen individual.Takes account of richness and evenness.
Combined$Type1 = factor(Combined$Type1, levels = c("cloaca","soil","egg","mid sub","mid egg","final egg","in egg","final sub", "hatchling","final turt","One_year_wild")) 
Richness_shannon <- ggplot(Combined) +geom_boxplot(aes(x=Type1, y =diversity_shannon, color=Sample_cat)) + scale_x_discrete(limits = c("cloaca","soil","egg","mid egg","mid sub","final egg","in egg", "final sub","hatchling","final turt","One_year_wild")) + geom_jitter(aes(x=Type1, y =diversity_shannon)) + theme_classic() + ylab("Richness") + xlab("Sample Type") + theme(axis.title.x = element_text(size=16), axis.text.x = element_text(size=11), axis.title.y = element_text(size=16))


## see paper for example https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0371-6#Fig6
#add pairwise contarcats --eg:a, b Letters represent pairwise contrasts (P < 0.01) of sample type means of woodlarks (lower case grey) and skylarks (capital red)

##Add ANOVA statistics --Pairwise Tukey-Kramer contrasts among sample types


##---PD whole tree  (i.e. Phylogenetic Diversity) 
#  Alpha diversity is the diversity within a sample or group of samples
# Phylogenetic diversity (PD) – takes into consideration the phylogeny of microbes 
# to estimate diversity across a tree.
##See alpha diversity info at https://www.evolution.unibas.ch/walser/bacteria_community_analysis/2015-02-10_MBM_tutorial_combined.pdf

##see tutorial https://kembellab.ca/r-workshop/biodivR/SK_Biodiversity_R.html
#install.packages("picante", dependencies=TRUE)
library(picante)
library(ape)
library(vegan)
library(permute)
library(lattice)
library(nlme)
Rare_3000
##Load community data ---Ecological community data consist of observations of the (relative) 
# abundance of species in different samples. The format for community data is a data.frame 
# with samples in the rows and species in the columns. 
Rare_3000 #phyloseq-class experiment-level object
str(ASVs) #'data.frame'
class(Rare_3000_data) # "data.frame"
# get the dimension of the community object (rows x columns)
dim(Rare_3000_data) #  77107    14
#check to make sure our rows and columns have reasonable-looking names.
rownames(Rare_3000_data)
# take a peek at the data (just the first five rows/columns)
Rare_3000_data[1:5, 1:5]
head(colnames(Rare_3000_data)) # "OTU"       "Sample"    "Abundance" "SampleID"  "Type"      "Type1" 

head(colnames(ASVs))
# check total abundance in each sample
apply(Rare_3000_data, 1, sum) # Error in FUN(newX[, i], ...) : invalid 'type' (character) of argument
Taxa
ASV_phyla
ASVs
# colnames(Combined)
# names(Combined)
# obs <- Combined$observed
# sample <- Combined$SampleID
# type <- Combined$Type1
# phy_tree(Combined)  #Error in access(physeq, "phy_tree", errorIfNULL) : phy_tree slot is empty.
# str(Combined)  #data.frame':	83 obs. of  29 variables

str(ASV_phyla)  #Formal class 'phyloseq' [package "phyloseq"] with 5 slots
colnames(Rare_3000)   #NULL
names(Rare_3000)      #NULL
phy_tree <- phy_tree(Rare_3000)   #Phylogenetic tree with 929 tips and 928 internal nodes.
#Tip labels: ASV1, ASV2, ASV3, ASV4, ASV5, ASV6, ...
#Rooted; includes branch lengths.
str(phy_tree)  #List of 4

Rare_3000_data<-psmelt(Rare_3000)
str(Rare_3000_data)  #'data.frame':	77107 obs. of  14 variables:
colnames(Rare_3000_data)
row.names(Rare_3000_data$SampleID)
ASVs<-Rare_3000_data$OTU
abun<-Rare_3000_data$Abundance
type<-Rare_3000_data$Type1
sample <- Rare_3000_data$SampleID

PD_Tree <- pd(Rare_3000, phy_tree, include.root=TRUE )  #Error in samp > 0 : comparison (6) is possible only for atomic and list types
str(Rare_3000@phy_tree)

## Calculate Faith's Phylogenetic Diversity 
#---Calculate the sum of the total phylogenetic branch length for one or multiple samples.
#See tutorial https://search.r-project.org/CRAN/refmans/picante/html/pd.html
library(picante)
#Usage:   pd(samp, tree, include.root=TRUE)
#Arguments: samp = commynuity data matrix; tree = A phylo tree object; 
#include.root	= Should the root node be included in all PD calculations (default = TRUE)
#Returns a dataframe of the PD and species richness (SR) values for all samples
pd(Rare_3000, phy_tree, include.root=TRUE) #Error in samp > 0 : comparison (6) is possible only for atomic and list types

## Try another method for pd whole tree #see https://rdrr.io/github/twbattaglia/btools/man/estimate_pd.html#heading-1
# estimate_pd function to estimate phylogenetic diversity from phyloseq object with an OTU table and phylogenetic tree slot.
str(Rare_3000)  # Formal class 'phyloseq' [package "phyloseq"] with 5 slots
# Estimate the Faiths phylogenetic diverstiy from an OTU table and phylogenetic tree
estimate_pd(Rare_3000) #Error in estimate_pd(Rare_3000) : could not find function "estimate_pd" 
#Not sure which package this function is in... maybe phyloseq. Also on website: twbattaglia/btools: A suite of R function for all types of microbial diversity analyses
library(phyloseq)
install.packages("twbattaglia/btools")
install.packages("twbattaglia")

# Check package microbiome for pd whole tree

#phyloMDA：An R package for phylogeny-aware microbiome data analysis. See https://github.com/liudoubletian/phyloMDA
install.packages("MGLM")
library(MGLM)
library(plyr)
install.packages("casper")
library(caper)
install.packages("genlasso")
library(genlasso)
library(magrittr)
library(foreach)
library(ape)
install.packages("miLineage")
library(miLineage)
library(ggplot2)
library(dplyr)
library(readxl)
library(methods)
library(BiocManager)
library(phyloseq)
# install.packages("ggtree")
#if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("ggtree")
library(ggtree)
# install.packages("adaANCOM")
# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install("adaANCOM")
library(adaANCOM) # Not avail.

install.packages("devtools")  
devtools::install_github("liudoubletian/phyloMDA") 
library(phyloMDA)
library(phyloseq); packageVersion("phyloseq")
library(phyloMDA); packageVersion("phyloMDA")

plot_tree(Rare_3000, "treeonly", nodeplotblank, label.tips="taxa_names")
tree <- phy_tree(Rare_3000)


##Use ggtree function in ggplot2 package to visualize tree. See https://guangchuangyu.github.io/ggtree-book/chapter-ggtree.html
phy_tree <- Rare_3000@phy_tree
phy_tree # Phylogenetic tree with 929 tips and 928 internal nodes.
ggtree(phy_tree)
ggtree(phy_tree, layout="circular")
ggtree(phy_tree, branch.length='none')
ggtree(phy_tree, layout="fan", open.angle=120)

ggplot(phy_tree) + geom_tree() + theme_tree()
plot_tree(Rare_3000)


### Try calculatePD: Calculate Faith's Phylogenetic Diversity 
#see https://rdrr.io/github/RIVM-IIV-Microbiome/biomeUtils/man/calculatePD.html
#see also https://rdrr.io/github/RIVM-IIV-Microbiome/biomeUtils/src/R/calculatePD.R
# A wrapper around picante for phyloseq objects to calculate and add Phylogenetic Diversity and Species Richness to sample data.
#Usage: calculatePD(x, justDF = FALSE, include_root = TRUE)
# devtools::install_github("RIVM-IIV-Microbiome/biomeUtils")
library(biomeUtils)
calculatePD(Rare_3000, justDF=FALSE, include_root=TRUE)

# Faiths phylo div in one sample: see http://scikit-bio.org/docs/0.4.2/generated/generated/skbio.diversity.alpha.faith_pd.html



## See package abdiv at https://www.quantargo.com/help/r/latest/packages/abdiv/0.2.0
# install.packages("abdiv")
library(abdiv)
#Usage: faith_pd(x, tree, x_labels = NULL)
# x	= A numeric vector of species counts or proportions, or a logical vector of species presence/absence.
# tree = A phylogenetic tree object.
# tree_labels = A character vector of species labels for x.
phy_tree # Phylogenetic tree with 929 tips and 928 internal nodes.

# Faith's phylogenetic diversity for whole tree is equal to the sum of the
# branch lengths.
sum(phy_tree$edge.length) # 42.25105
faith_pd(c(), phy_tree) #  Length of x does not match number of tips in tree.

# Can use named vector or additional argument to match species to tree.
phy_tree$tip.label
faith_pd(c(0, 0, 0, 10, 12), phy_tree)
faith_pd(c(d=10, e=12), phy_tree)
faith_pd(c(10, 12), phy_tree, c("d", "e"))
faith_pd(ASVs, phy_tree, x_labels = NULL)
str(phy_tree)
str(Rare_3000) # Formal class 'phyloseq' [package "phyloseq"] with 5 slots
treetips <- Rare_3000@phy_tree@tip.label # Error: trying to get slot "tip.label" from an object (class "phylo") that is not an S4 object

phy_tree$tip.label
faith_pd(929, phy_tree)





## See tutroial: https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/alpha-diversities.html

# Check how much data you have after rarifying
print(Rare_3000)  #929 taxa and 83 samples

# quick check for sampling depth. his plot is only for a quick and dirty check for reads per samples.
barplot(sample_sums(Rare_3000), las =2)

p.rar <- plot_taxa_prevalence(Rare_3000, "Phylum")
# Error ... look for taxa prevalence plot from previous section of the tutorial.

library(ggplot2)
library(microbiome) # data analysis and visualisation
library(phyloseq) # also the basis of data object. Data analysis and visualisation
# install.packages("microbiomeutilities")
# if (!requireNamespace("BiocManager", quietly = TRUE))
library(microbiomeutilities) # some utility tools ##error
library(RColorBrewer) # nice color options
# install.packages("ggpubr")
# if (!requireNamespace("BiocManager", quietly = TRUE))
library(ggpubr) # publication quality figures, based on ggplot2
library(DT) # interactive tables in html and markdown
library(data.table) # alternative to data.frame
library(dplyr) # data handling  

## 5.2.2 Phylogenetic diversity
# Phylogenetic diversity is calculated using the picante package.
library(picante)
rar.asvtab <- as.data.frame(Rare_3000@otu_table)
rar.asvtab
rar.tree <- Rare_3000@phy_tree
rar.tree
# hmp.meta from previous code chunks
# We first need to check if the tree is rooted or not 
Rare_3000@phy_tree #Rooted; includes branch lengths.
# it is a rooted tree....
df.pd <- pd(t(rar.asvtab), rar.tree, include.root=T) 
# t(ou_table) transposes the table for use in picante and the tre file comes from the first code chunk we used to read tree file (see making a phyloseq object section).
## Error in UseMethod("is.rooted") : 
#no applicable method for 'is.rooted' applied to an object of class "NULL"
#  In drop.tip(tree, treeabsent) : drop all tips of the tree: returning NULL

##Try treeio package to re-root tree, S3 method for treedata: https://docs.ropensci.org/treeio/reference/root-method.html
library(treeio)
library(tidytree)
library(dplyr)
# S3 method for treedata
#Usage: root(phy, outgroup, node = NULL, edgelabel = TRUE, ...)
root(rar.tree, rar.asvtab, node = NULL, edgelabel = TRUE)
# Phylogenetic tree with 929 tips and 928 internal nodes.
#Tip labels:
#  ASV1, ASV2, ASV3, ASV4, ASV5, ASV6, ...
#Rooted; includes branch lengths.
#Re-try pd. Still same error msg.
pd <- pd(t(rar.asvtab), rar.tree, include.root=T) 
#why is object class NULL? check str
str(rar.tree) #List of 4 --might be the issue?
str(rar.asvtab) #'data.frame':	83 obs. of  929 variables:
#Try converting tree to different format (see https://yulab-smu.top/treedata-book/chapter2.html)
tbl.tree <- as_tibble(rar.tree)
tbl.tree #gives table with parent/node/branch.length/label
#The tbl_tree object can be converted back to a phylo object using the as.phylo() method.
phylo_tree <- as.phylo(tbl.tree)
phylo_tree #Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'#Also defined by ‘tidytree’
str(phylo_tree) #list of 4.... still same format
#2.1.2 The treedata object
# The tidytree package defines treedata class to store a phylogenetic tree with associated data. After mapping external data to the tree structure, the tbl_tree object can be converted to a treedata object.
pd.tree <- as.treedata(rar.tree)
pd.tree # 'treedata' S4 object'.
#re-try with treedata format
pd <- pd(t(rar.asvtab), pd.tree, include.root=T) #Error in tree$edge.length : $ operator not defined for this S4 class
str(Rare_3000)
pd <- pd(t(Rare_3000$otu_table), Rare_3000$phy_tree, include.root=T) #Still get same error msg: $ operator not defined for this S4 class





#### FIGURE 1B. Beta Diversity ####
##---Unifrac (Unique Fraction Metric)

##UniFrac is a β-diversity measure that uses phylogenetic information to compare environmental samples. 
#UniFrac, coupled with standard multivariate statistical techniques including principal coordinates analysis (PCoA), 
#identifies factors explaining differences among microbial communities. (https://www.nature.com/articles/ismej2010133)

##see pg 35 https://www.evolution.unibas.ch/walser/bacteria_community_analysis/2015-02-10_MBM_tutorial_combined.pdf
#UniFrac distances are based on the fracGon of branch length shared between two communities within a phylogenetic 
#tree constructed from the 16S rRNA gene sequences from all communities being compared.
#  With unweighted UniFrac, only the presence or absence of lineages are considered (community membership).
#  With weighted UniFrac, branch lengths are weighted based on the relaGve abundances of lineages within 
#communities (community structure).

##see tutorial https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/beta-diversity-metrics.html
#Beta-diversity: Measures for differences between samples from different groups 
#to identify if there are differences in the overall community composition and structure.

library(microbiome) #data analysis and visualization
library(phyloseq) ##data analysis and visualization
library(dplyr) #data handling
library(RColorBrewer) #nice color options
# library(ggpubr) #publication quality figures for ggplot2




### Use a colour fade gradient for sample types
library(RColorBrewer)
par(mar=c(3,4,2,2))
display.brewer.all() #Look at colour palette options
# Usage: + scale_fill_brewer(palette="Dark2") #add to very end of plot code outside of brackets for other specifications
## Issue... all the fade colour scheme options in RColorBrewer pkg only allow up to n=9, need colours for n=11
##Try gradient clours for scatter plots...
#usage: +scale_color_gradient(low="blue", high="red") ----Error: Discrete value supplied to continuous scale
##look into viridis package for more options...

## Make groups an ordered factor: 1 = cloaca, 2=egg, 3=mid egg..... unsure how to do this in S4 class system
### Phyloseq ordering code from Marissa
#pseq@sam_data$Treatment <- factor(pseq@sam_data$Treatment,
#                                  levels = c("Control", "Probiotics","Probiotics + HT", "High temperature", "NA"))
#plot_ordination(pseq, ord, color = "Treatment", shape = "Age") + geom_point(size = 4)

Rare_3000@sam_data$Type1 <- factor(Rare_3000@sam_data$Type1, levels = c("cloaca","soil", "egg", "mid egg","mid sub", "final egg","in egg",  "final sub", "hatchling", "final turt", "One_year_wild_body", "One_year_wild_cloaca"))
###IT WORKED!!
plot_ordination(Rare_3000, ord, color = "Type1", shape = "Sample_cat") + geom_point(size = 4)



summarize_phyloseq(Rare_3000)
#Unweighted UniFrac only considers the presence/absence of taxa between sample pairs.
ordu.unwt.uni <- ordinate(Rare_3000, "PCoA", "unifrac", weighted=F)
unwt.unifrac <- plot_ordination(Rare_3000, ordu.unwt.uni, color="Type1") + scale_fill_brewer(palette="YlGnBu") + ggtitle("Unweighted UniFrac") + stat_ellipse(type = "norm", linetype = 1) + theme_classic() + theme_classic() + geom_point(size=4) + theme(legend.title = element_text(size=16), legend.text = element_text(size=14) + theme(axis.title.x = element_text(size=16), axis.text.x = element_text(size=9), axis.title.y = element_text(size=16))) + scale_color_brewer(palette="YlGnBu")
print(unwt.unifrac)
unwt.unifrac + scale_color_brewer(palette="YlGnBu")  ##Warning: n too large, allowed maximum for palette YlGnBu is 9
unwt.unifrac + scale_color_brewer(palette="Paired") 

#Weighted Unifrac does consider the abundances of different taxa.
ps1.rel <- microbiome::transform(Rare_3000, "compositional")
ps1.rel  #phyloseq-class experiment-level object.  [ 929 taxa and 83 samples ]
ordu.wt.uni <- ordinate(ps1.rel , "PCoA", "unifrac", weighted=T)
ordu.wt.uni #"No correction was applied to the negative eigenvalues"
# check for Eigen values 
# barplot(ordu.unwt.uni$values$Eigenvalues[1:10])
wt.unifrac <- plot_ordination(ps1.rel,ordu.wt.uni, color="Type1", shape = "Sample_cat") + ggtitle("Weighted UniFrac") + stat_ellipse(type = "norm", linetype = 1) + theme_classic() + theme_classic() + geom_point(size=4) + theme(legend.title = element_text(size=16), legend.text = element_text(size=14) + theme(axis.title.x = element_text(size=16), axis.text.x = element_text(size=9), axis.title.y = element_text(size=16))) + scale_color_brewer(palette="Spectral") 
print(wt.unifrac)



####Calculate UniFrac Distance for all sample-pairs in a phyloseq-class object.
##See tutorial https://www.rdocumentation.org/packages/phyloseq/versions/1.16.2/topics/UniFrac
#UniFrac() function accesses the abundance (otu_table-class) and a phylogenetic tree (phylo-class) data within an experiment-level (phyloseq-class) object. 
#phyloseq-class objuct must contain at minimum a phylogenetic tree (phylo-class) and contingency table (otu_table-class).
#Weighted-UniFrac takes into account the relative abundance of species/taxa shared between samples, whereas unweighted-UniFrac only considers presence/absence.
# normalized (argument)--should the output be normalized such that values range from 0 to 1 independent of branch length values? Default is TRUE. Note that (unweighted) UniFrac is always normalized by total branch-length, and so this value is ignored when weighted == FALSE.
# parallel ---should execute calculation in parallel, using multiple CPU cores simultaneously? This can dramatically hasten the computation time for this function. However, it also requires that the user has registered a parallel ``backend'' prior to calling this function. Default is FALSE. If FALSE, UniFrac will register a serial backend so that foreach::%dopar% does not throw a warning.
# fast ---Do you want to use the ``Fast UniFrac'' algorithm? Implemented natively in the phyloseq-package. TRUE is only supported option
Weighted_UniFrac <- UniFrac(Rare_3000, weighted=FALSE, normalized=TRUE, parallel=FALSE, fast=TRUE)
UniFrac_table <- phyloseq(ASVs, phy_tree)  #phyloseq-class experiment-level object (that has been pruned and comprises the minimum arguments necessary for UniFrac().)
#otu_table()   OTU Table:         [ 929 taxa and 83 samples ]
#phy_tree()    Phylogenetic Tree: [ 929 tips and 928 internal nodes ]



###UniFrac with subset of data

library(phyloseq)
#subset_samples(physeq, ...)
###Sub-sample only cloaca, egg, and soil sample types *Initial sources of egg microbiome
Sub_samp<-subset_samples(Rare_3000, Type1 %in% c("cloaca", "soil", "egg"))
Sub_samp  #Includes phy_tree() and is subset [ 929 taxa and 21 samples ]
sub_samp_ag <- aggregate_rare(Sub_samp, level = "Phylum", detection = .5/100, prevalence = 20/100) ##ISSUE is aggregate_rare() function 
# removes phy_tree from phyloseq object ...try to use sub_samp for UniFrac

Sub_samp  #phyloseq-class experiment-level object [ 929 taxa and 21 samples ]
ps2.rel <- microbiome::transform(Sub_samp, "compositional")
ps2.rel  #phyloseq-class experiment-level object.  [ 929 taxa and 21 samples ]
ordu.wt.uni_sub <- ordinate(ps2.rel , "PCoA", "unifrac", weighted=T)
ordu.wt.uni_sub #"No correction was applied to the negative eigenvalues"
# check for Eigen values 
# barplot(ordu.unwt.uni_sub$values$Eigenvalues[1:10])
wt.unifrac_sub <- plot_ordination(ps2.rel,ordu.wt.uni_sub, color="Type1") + ggtitle("Egg and Microbial Sources") + stat_ellipse(type = "norm", linetype = 1) + theme_classic() + geom_point(size=4) + theme(legend.title = element_text(size=16), legend.text = element_text(size=14) + theme(axis.title.x = element_text(size=16), axis.text.x = element_text(size=9), axis.title.y = element_text(size=16)))
wt.unifrac_sub 
#   consider removing outlier cloaca sample --figure out which sample is the outlier

###ADD HATCHLING --when eggs laid look more like soil, then more like turtles once hatched. 
# ADD in egg. interior of egg seems unaffected by exterior of egg/env. Add final turt --pattern continues, turtles have specific microbiome from env, but egg is affected by env microbes.
# ADD post release turtles. post rel begin to shift toward cloaca... wild turtles have different microbes than captive.. or old vs young turtles have different microbes
#ADD final sub and final egg. environment seems to determine host microbiome ---egg looks like soil, final egg and hatchling look like final sub, and post release shifts from final sub towards cloaca/wild maternal turtles.
early_sub_samp <-subset_samples(Rare_3000, Type1 %in% c("hatchling", "egg", "cloaca", "soil","mid egg", "mid sub", "final egg", "final sub", "One_year_wild"))
early_sub_samp #[ 929 taxa and 16 samples ]
ps1.5.rel <- microbiome::transform(early_sub_samp, "compositional")
ps1.5.rel  #phyloseq-class experiment-level object.  [ 929 taxa and 21 samples ]
ordu.wt.uni_early <- ordinate(ps1.5.rel , "PCoA", "unifrac", weighted=T)
ordu.wt.uni_early #"No correction was applied to the negative eigenvalues"
wt.unifrac_early <- plot_ordination(ps1.5.rel, ordu.wt.uni_early, color="Type1") + ggtitle("Early Development") + stat_ellipse(type = "norm", linetype = 1) + theme_classic() + theme_classic() + geom_point(size=4) + theme(legend.title = element_text(size=16), legend.text = element_text(size=14) + theme(axis.title.x = element_text(size=16), axis.text.x = element_text(size=9), axis.title.y = element_text(size=16)))
wt.unifrac_early


###Turtle ONLY *Wild vs captive
Turt_sub_samp <-subset_samples(Rare_3000, Type1 %in% c("final egg", "hatchling", "One_year_wild", "final turt", "cloaca"))
Turt_sub_samp #[ 929 taxa and 16 samples ]
ps3.rel <- microbiome::transform(Turt_sub_samp, "compositional")
ps3.rel  #phyloseq-class experiment-level object.  [ 929 taxa and 21 samples ]
ordu.wt.uni_turt <- ordinate(ps3.rel , "PCoA", "unifrac", weighted=T)
ordu.wt.uni_turt #"No correction was applied to the negative eigenvalues"
wt.unifrac_turt <- plot_ordination(ps3.rel, ordu.wt.uni_turt, color="Type1") + ggtitle("Wild vs Captive Turtle") + stat_ellipse(type = "norm", linetype = 1) + theme_classic() + geom_point(size=4) + theme(legend.title = element_text(size=16), legend.text = element_text(size=14) + theme(axis.title.x = element_text(size=16), axis.text.x = element_text(size=9), axis.title.y = element_text(size=16)))
wt.unifrac_turt ##new added final

###Host ONLY *WPT microbiome succession 
Host_sub_samp <-subset_samples(Rare_3000, Type1 %in% c("hatchling", "One_year_wild", "final turt", "egg", "mid egg", "final egg", "in egg"))
Host_sub_samp # [ 929 taxa and 54 samples ]
ps4.rel <- microbiome::transform(Host_sub_samp, "compositional")
ps4.rel  #phyloseq-class experiment-level object.  [ 929 taxa and 21 samples ]
ordu.wt.uni_host <- ordinate(ps4.rel , "PCoA", "unifrac", weighted=T)
ordu.wt.uni_host #"No correction was applied to the negative eigenvalues"
wt.unifrac_host <- plot_ordination(ps4.rel, ordu.wt.uni_host, color="Type1") + ggtitle("Host Through Development") + stat_ellipse(type = "norm", linetype = 1) + theme_classic() + theme_classic() + geom_point(size=4) + theme(legend.title = element_text(size=16), legend.text = element_text(size=14) + theme(axis.title.x = element_text(size=16), axis.text.x = element_text(size=9), axis.title.y = element_text(size=16)))
wt.unifrac_host

###Egg ONLY
Egg_sub_samp <-subset_samples(Rare_3000, Type1 %in% c("egg", "mid egg", "final egg", "in egg"))
Egg_sub_samp # [ 929 taxa and 33 samples ]
ps5.rel <- microbiome::transform(Egg_sub_samp, "compositional")
ps5.rel  #phyloseq-class experiment-level object.  [ 929 taxa and 21 samples ]
ordu.wt.uni_egg <- ordinate(ps5.rel , "PCoA", "unifrac", weighted=T)
ordu.wt.uni_egg #"No correction was applied to the negative eigenvalues"
wt.unifrac_egg <- plot_ordination(ps5.rel, ordu.wt.uni_egg, color="Type1") + ggtitle("Egg beta diverity") + stat_ellipse(type = "norm", linetype = 1) + theme_classic() + geom_point(size=4) + theme(legend.title = element_text(size=16), legend.text = element_text(size=14) + theme(axis.title.x = element_text(size=16), axis.text.x = element_text(size=9), axis.title.y = element_text(size=16)))
wt.unifrac_egg

##Egg and environment though incubation ## IMPORTANT; succession
sub_samp_egg.env<-subset_samples(Rare_3000, Type1 %in% c("egg", "soil", "mid egg", "mid sub", "final egg", "final sub","in egg"))
sub_samp_egg.env
ps6.rel <- microbiome::transform(sub_samp_egg.env, "compositional")
ps6.rel  #phyloseq-class experiment-level object.  [ 929 taxa and 21 samples ]
ordu.wt.uni_egg.env <- ordinate(ps6.rel , "PCoA", "unifrac", weighted=T)
ordu.wt.uni_egg.env #"No correction was applied to the negative eigenvalues"
wt.unifrac_egg.env <- plot_ordination(ps6.rel, ordu.wt.uni_egg.env, color="Type1", shape = "Sample_cat") + ggtitle("Egg and Environment") + stat_ellipse(type = "norm", linetype = 1) + theme_classic() + theme_classic() + geom_point(size=4) + theme(legend.title = element_text(size=16), legend.text = element_text(size=14) + theme(axis.title.x = element_text(size=16), axis.text.x = element_text(size=9), axis.title.y = element_text(size=16)))
wt.unifrac_egg.env
# consider adding shape = time. add time as a new column in data ##need to find out from andrea timeline, days post collection


####FUll FIG 1 A and B####
## See turtorial for using grid.arrange function: https://cran.r-project.org/web/packages/egg/vignettes/Ecosystem.html
#note: still missing alpha pd, add in when created
library(vegan)
library(ggplot2)
library(gridExtra)
 install.packages("gridExtra")
 
grid.arrange(Richness_shannon) ##ADD alpa_pd
wt.unifrac
grid.arrange(wt.unifrac_sub,wt.unifrac_egg.env, wt.unifrac_host, wt.unifrac_turt, nrow = 2)

###Inval tutorial ****
install.packages("indicspecies")
library(indicspecies)

#load data
Rare_3000
#load objects
OTU = Rare_3000@otu_table
Tax = Rare_3000@tax_table
Metadata = Rare_3000@sam_data
Tree = Rare_3000@phy_tree

#Extract abundance matrix ----
#from the phyloseq object using phyloseq
OTU1 = as(OTU, "matrix")
Metadata = as(Metadata, "matrix")
write.csv(OTU1, file="turt_meta.cvs",row.names=TRUE)
write.csv(Metadata, file="turt_meta.cvs",row.names=TRUE ) 

write.table(OTU1,file="data_table.csv",sep=",",dec = " ")
write.table(Metadata,file="data_table_Meta.csv",sep=",",dec = " ")
####Format to example data and reload below for actual test 
#reload edited table
pc_FUN = read.csv("data_table.csv", header= TRUE)
View(pc_FUN)


### use colour brewer set1 for palette

####  FIGURE 2 Core bubble plot ####
## Maybe put bubble plot next to taxa table
print(Taxa)
phy_tree(Rare_3000)
## Maybe only include 10 most abundant genera in plot

## microbiome package, command to make phyloseq object of core bacteria ---make object with core, then bubble plot
#see https://r-graph-gallery.com/320-the-basis-of-bubble-plot.html
#A bubble plot is a scatterplot with an added third dimension: the value of an additional numeric variable is represented 
# through the size of the dots. With ggplot2, bubble chart are built using the geom_point() function. 
# At least three variable must be provided to aes(): x, y and size. The legend will automatically be built by ggplot2.
## Ask Danielle for bubble plot. Choose which bacteria ---inval tutorial?
library(ggplot2)
library(dplyr)
# Most basic bubble plot
#ggplot(data, aes(x=gdpPercap, y=lifeExp, size = pop)) + geom_point(alpha=0.7)
# Most basic bubble plot
ggplot(Rare_3000_data, aes(x = Type1 , y = Phylum, size = Abundance)) + geom_point(alpha=0.7)
View(Rare_3000_data)
#Phylum level
Rare_3000_data %>%
  arrange(desc(Abundance))  %>%
  ggplot(aes(x= Type1, y= Phylum, size = Abundance)) +
  geom_point(alpha=0.5) +  xlab("Sample Type") +
  scale_size(range = c(.1, 24), name="ASV Abundance")


View(gen_abun)

install.packages("hrbrthemes")
library(hrbrthemes)

#Order plus Phylum 
Rare_3000_data %>% 
  ggplot(aes(x= Type1, y= Order, size = Abundance, color = Phylum)) + 
  geom_point(aes(size = Abundance, fill = Phylum), alpha = 0.75, shape = 21) + theme_classic() + xlab("Sample Type") +
  scale_size(limits = c(1, 1145),range = c(.1, 24), breaks = c(1, 200, 400, 600, 800, 1000), name="ASV Abundance") + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(size = 12, angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(size = 9), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 12), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right") + scale_x_discrete(limits = c("egg","mid egg","final egg","in egg", "hatchling","final turt","One_year_wild", "cloaca","soil","mid sub","final sub") +  scale_colour_ipsum()) 
#Warning message: Removed 61700 rows containing missing values (`geom_point()`).... this is from removing 0 values using scale_size(limits = c() function



### top taxa only using MarissaLag Github bubble plot repository https://github.com/MarissaLag/Phyloseq-and-microbiome-analysis/blob/main/inval_tutorial/INVAL_TUTORIAL_FOR_MICROBIOME_2023.R 
pseq_gen <- microbiome::aggregate_rare(Rare_3000, level = "Genus", detection = 50/100, prevalence = 70/100)
gen_abun <- psmelt(pseq_gen)


top5G.names = sort(tapply(taxa_sums(Rare_3000), tax_table(pseq)[, "Genus"], sum), TRUE)[1:5]
top5G = subset_taxa(pseq, Genus %in% names(top5G.names))
top5G <- psmelt(top5G)

top3P.names = sort(tapply(taxa_sums(Rare_3000), tax_table(pseq)[, "Phylum"], sum), TRUE)[1:3]
top3P = subset_taxa(pseq, Phylum %in% names(top3P.names))
top3P <- psmelt(top3P)

top10O.names = sort(tapply(taxa_sums(Rare_3000), tax_table(pseq)[, "Order"], sum), TRUE)[1:10]
top10O = subset_taxa(pseq, Order %in% names(top10O.names))
top10O <- psmelt(top10O)

top5O.names = sort(tapply(taxa_sums(Rare_3000), tax_table(pseq)[, "Order"], sum), TRUE)[1:5]
top5O = subset_taxa(pseq, Order %in% names(top5O.names))
top5O <- psmelt(top5O)

##Top orders bubble plot
top5O %>% 
  ggplot(aes(x= Type1, y= Order, size = Abundance, color = Phylum)) + 
  geom_point(aes(size = Abundance, fill = Phylum), alpha = 0.75, shape = 21) + theme_classic() + xlab("Sample Type") +
  scale_size(limits = c(1, 1145),range = c(.1, 24), breaks = c(1, 200, 400, 600, 800, 1000), name="ASV Abundance") + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(size = 12, angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(size = 9), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 12), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right") + scale_x_discrete(limits = c("egg","mid egg","final egg","in egg", "hatchling","final turt","One_year_wild", "cloaca","soil","mid sub","final sub") +  scale_colour_ipsum()) 

#Create objects ----

OTU = pseq@otu_table
Tax = pseq@tax_table
Metadata = pseq@sam_data
Tree = pseq@phy_tree

pseq_fam <- microbiome::aggregate_rare(pseq, level = "Family", detection = 50/100, prevalence = 70/100)

#convert to compositional data

pseq.fam.rel <- microbiome::transform(pseq_fam, "compositional")

pseq.core <- core(pseq.fam.rel, detection = .1/100, prevalence = 50/100)

pseq.core <- microbiome::transform(pseq.core, "compositional")

gen_abun %>%
  ggplot(aes(x= Type1, y= Genus, size = Abundance, color = Genus)) +
  geom_point(aes(size = Abundance, fill = Genus), alpha = 0.75, shape = 21) + theme_classic() + xlab("Sample Type") +
  scale_size(limits = c(1, 1145),range = c(.1, 24), breaks = c(1, 200, 400, 600, 800, 1000), name="ASV Abundance") + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(size = 12, angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(size = 9), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 12), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right") + scale_x_discrete(limits = c("egg","mid egg","final egg","in egg", "hatchling","final turt","One_year_wild", "cloaca","soil","mid sub","final sub")) 


#Trying to get genus and phyla... error... not merging well
pseq_phyla <- microbiome::aggregate_rare(Rare_3000, level = "Phylum", detection = 50/100, prevalence = 70/100)
phyla_abun <- psmelt(pseq_phyla)

phyla_abun
gen_abun

Combined_abun <- merge(gen_abun, phyla_abun, by = "OTU")
write.csv(Combined_abun, file = "combinedABUN.csv")

Combined_abun<-read.csv("combinedABUN.csv")
Combined_abun<-merge(gen_abun,phyla_abun, by = "Sample")

Combined_abun %>%
  ggplot(aes(x= Type1, y= Genus, size = Abundance, color = Phylum)) +
  geom_point(aes(size = Abundance, fill = Genus), alpha = 0.75, shape = 21) + theme_classic() + xlab("Sample Type") +
  scale_size(limits = c(1, 1145),range = c(.1, 24), breaks = c(1, 200, 400, 600, 800, 1000), name="ASV Abundance") + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(size = 12, angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(size = 9), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 12), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right") + scale_x_discrete(limits = c("egg","mid egg","final egg","in egg", "hatchling","final turt","One_year_wild", "cloaca","soil","mid sub","final sub")) 

#for formatting code see: https://jkzorz.github.io/2019/06/05/Bubble-plots.html
#A few classic improvement:
#use of the viridis package for nice color palette
#use of theme_ipsum() of the hrbrthemes package
#custom axis titles with xlab and ylab
#add stroke to circle:change shape to 21 and specify color (stroke) and fill
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)
#hrbrthemes::import_roboto_condensed() 

### show top 10*** Goal: https://www.researchgate.net/figure/Bubble-plot-showing-the-relative-abundance-of-the-top-25-most-abundant-taxa-found-across_fig2_361270406
#see https://github.com/joey711/phyloseq/issues/1575
###use relative % genus instead of ASV
devtools::install_github("gmteunisse/fantaxtic")
require("fantaxtic")
require("phyloseq")
data(Rare_3000_data)
top <- top_taxa(Rare_3000_data, 
                tax_level = "Order", 
                n_taxa = 10,
                grouping = "Type1")
plot_nested_bar(top$ps_obj, top_level = "Phylum", nested_level = "Order")

##see https://jkzorz.github.io/2019/06/05/Bubble-plots.html
#might need to import taxa data to wide format first then change to long format? upload csv file with species as columns and samples as rows

#https://rdrr.io/github/microbiome/microbiome/man/top_taxa.html
top25 <- top_taxa(Rare_3000, n=25) ##lists top ASVs... not sure how to get top Orders
str(top25) #ps_obj ...need data.frame structure for bubble plot data
#Use psmelt function to change the data format from a phyloseq to a data.frame 
top25_taxa<- psmelt(top25) #Error in access(object, "otu_table", errorIfNULL) : otu_table slot is empty.


##Try tax_glom



#### FIGURE 3 NTI ####
#**SWAP FOR INDICATOR ANALYSIS?
##  (Nearest taxon index) 
## see tutoriual https://pedrohbraga.github.io/CommunityPhylogenetics-Workshop/CommunityPhylogenetics-Workshop.html
## Phylogenetic trees (i.e., evolutionary tree or cladogram) are branching diagrams illustrating the evolutionary relationships among taxa.
#The nearest taxon index (NTI) is a standardized measure of the mean phylogenetic distance to the nearest taxon in each sample/community (MNTD). 
# Check for needed packages, and install the missing ones
required.libraries <- c("ape", "picante", 
                        "pez", "phytools",
                        "vegan", "adephylo", 
                        "phylobase", "geiger", 
                        "mvMORPH", "OUwie", 
                        "hisse", "BAMMtools",
                        "phylosignal", "Biostrings",
                        "devtools","ggplot2", 
                        "kableExtra", "betapart", "gridExtra",
                        "reshape2")

needed.libraries <- required.libraries[!(required.libraries %in% installed.packages()[,"Package"])]
if(length(needed.libraries)) install.packages(needed.libraries)

# Load all required libraries at once
lapply(required.libraries, require, character.only = TRUE)

### Install ggtree from BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree")

##use picante pkg
library(picante)

phy_tree

### Make a tree file?  
PhylogeneticTree<-phy_tree(Rare_3000)
PhylogeneticTree
write.csv(PhylogeneticTree, file="PhylogeneticTree.csv")
##Error in as.data.frame.default(x[[i]], optional = TRUE, stringsAsFactors = stringsAsFactors) : cannot coerce class ‘"phylo"’ to a data.frame

get_NRI_NTI(obj, ...)

## See tutorial for calculating related phylogenetic alpha metric https://rdrr.io/github/YuLab-SMU/MicrobiotaProcess/man/get_NRI_NTI-methods.html

## S4 method for signature 'matrix'
get_NRI_NTI(
  Rare_3000,
  mindepth,
  sampleda,
  tree,
  metric = c("PAE", "NRI", "NTI", "PD", "HAED", "EAED", "IAC", "all"),
  abundance.weighted = FALSE,
  force = FALSE,
  seed = 123,
  ...
)  ##Error in get_NRI_NTI(Rare_3000, mindepth, sampleda, tree, metric = c("PAE",  : could not find function "get_NRI_NTI"
#cant find anything on internet abouty this function...

##see https://chiliubio.github.io/microeco_tutorial/model-based-class.html

#see https://rdrr.io/cran/iCAMP/man/NTI.p.html
install.packages("iCAMP")
library(iCAMP)
NTI.p(Taxa, dis, nworker = 4, memo.size.GB = 50,
      weighted = c(TRUE, FALSE), rand = 1000,
      check.name = TRUE, output.MNTD = c(FALSE, TRUE),
      sig.index=c("SES","NTI","Confidence","RC","all"),
      silent=FALSE)

### INDICATOR ANALYSIS
#see tutorial: https://cran.r-project.org/web/packages/indicspecies/vignettes/IndicatorSpeciesAnalysis.html
install.packages("indicspecies")
library(indicspecies)





#### FIGURE 4A. Composition plots of shared taxa from microbial sources ####
## Remove sample ID on x axis... maybe replace with sample type
## Look into changing colour schemes for bar plots
## make taxonomic composition plots of complete taxa rather than only shared 
comp_plot_all <- plot_composition(Rare_3000, group_by = "Type1", average_by = "Type1", otu.sort = "abundance", sample.sort = "Rare_3000_ordered")  + theme_classic() + scale_color_brewer(palette="Spectral") + theme(axis.text.x = element_text(angle = 45, vjust=0.5))
comp_plot_all 

###organize by categories
Combined$Type1 = factor(Combined$Type1, levels = c("cloaca", "egg", "in egg", "mid egg", "final egg", "hatchling", "final turt", "One_year_wild", "soil", "mid sub", "final sub"))
Comp_plot <- plot_composition(Rare_3000, group_by = "Type1") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank() + xlab("Sample")) 
Rare_3000_ordered = c("cloaca","soil", "egg", "mid egg","mid sub", "final egg","in egg",  "final sub", "hatchling", "final turt", "One_year_wild")

summarize_phyloseq(Rare_3000)

# see https://rdrr.io/github/microbiome/microbiome/man/plot_composition.html for code usage and argument options
#try to use sample.sort argument to order?
# make x axis titles vertical
## Plot all sample categories
All_phylum <- aggregate_rare(Rare_3000, level = "Phylum", detection = 1/100, prevalence = 50/100) 
All_phylum   #[ 5 taxa and 21 samples ]
Comp_plot_phylum <- plot_composition(All_phylum, group_by = "Type1", average_by = "Type1", otu.sort = "abundance", sample.sort = "Rare_3000_ordered") + theme_classic() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ theme(legend.title = element_text(size=16), legend.text = element_text(size=14), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16))
Comp_plot_phylum


#See tutorial https://david-barnett.github.io/microViz/articles/web-only/compositions.html
#devtools::install_github("david-barnett/microViz")
library("microViz")
library(patchwork)
comp_barplot(Rare_3000,)

#use phyloseq_extract_shared_otus function to subset phyloseq object to the OTUs that are present in all samples. 

###shared ASVs between eggs and microbe sources
ps.c.e.s <- phyloseq_extract_shared_otus(Rare_3000, samp_names = c("Sample-104", "Sample-113", "Sample-28", "Sample-32", "Sample-46", "Sample-53", "Sample-62", "Sample-64", "Sample-83", "Sample-116", "Sample-121", "Sample-31", "Sample-51", "Sample-59", "Sample-65", "Sample-70", "Sample-94", "Sample-103", "Sample-123", "Sample-125", "Sample-54")) 
ps.c.e.s
otu_table(ps.c.e.s)    #[9 taxa and 21 samples]
#use aggregate_rare function to combine rare taxa 
##at phylum level
ps.c.e.s_phylum<-aggregate_rare(ps.c.e.s, level = "Phylum", detection = 1/100, prevalence = 50/100) 
ps.c.e.s_phylum   #[ 5 taxa and 21 samples ]
Comp.c.e.s_phylum<- plot_composition(ps.c.e.s_phylum, group_by = "Type1") + theme_classic() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ theme(legend.title = element_text(size=16), legend.text = element_text(size=14), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16))
Comp.c.e.s_phylum
##at genus level
ps.c.e.s_genus<-aggregate_rare(ps.c.e.s, level = "Genus", detection = 1/100, prevalence = 50/100) 
ps.c.e.s_genus
Comp.c.e.s_genus<- plot_composition(ps.c.e.s_genus, group_by = "Type1")
Comp.c.e.s_genus

### Egg vs colaca  ##at phylum level
ps.e.c_phylum <- phyloseq_extract_shared_otus(Rare_3000, samp_names = c("Sample-116", "Sample-121", "Sample-31", "Sample-51", "Sample-59", "Sample-65", "Sample-70", "Sample-94", "Sample-104", "Sample-113", "Sample-28", "Sample-32", "Sample-46", "Sample-53", "Sample-62", "Sample-64", "Sample-83"))
ps.e.c_phylum<-aggregate_rare(ps.e.c_phylum, level = "Phylum", detection = 1/100, prevalence = 50/100)
Comp.e.c_phylum <- plot_composition(ps.e.c_phylum, group_by = "Type1") + theme_classic() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ theme(legend.title = element_text(size=16), legend.text = element_text(size=14), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16))
Comp.e.c_phylum
##at genus level
ps.e.c_genus <- phyloseq_extract_shared_otus(Rare_3000, samp_names = c("Sample-116", "Sample-121", "Sample-31", "Sample-51", "Sample-59", "Sample-65", "Sample-70", "Sample-94", "Sample-104", "Sample-113", "Sample-28", "Sample-32", "Sample-46", "Sample-53", "Sample-62", "Sample-64", "Sample-83"))
ps.e.c_genus<-aggregate_rare(ps.e.c_genus, level = "Genus", detection = 1/100, prevalence = 50/100)
Comp.e.c_genus <- plot_composition(ps.e.c_genus, group_by = "Type1") + theme_classic() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ theme(legend.title = element_text(size=16), legend.text = element_text(size=14), axis.title.x = element_text(size=16), axis.title.y = element_text(size=16))
Comp.e.c_genus


##Find way to include phyla info on genera plots



#### FIGURE 4B. Venn Diagram of shared ASVs ####
##Maybe make pannel of venn diagrams like Vanveelan et al 2017 https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0371-6#Fig6
##### use this code for venn diagram
library(devtools)
#install.packages("Git")
install_github("Russel88/MicEco")
library("MicEco")
#ps_venn(Rare_3000, group="Type1")

#usage:
# ps_venn(
#  ps,
#  group,
#  fraction = 0,
#  weight = FALSE,
#  relative = TRUE,
#  plot = TRUE,
#  ...
# )
##Shared ASVs for cloaca, egg, and soil
Venn_c.e.s <- ps_venn(Sub_samp, "Type1", font = 14, labels = list(cex = 1.3), fraction = 0.1)
Venn_c.e.s

##Sharred ASVs turtle only
Venn_turt <- ps_venn(Turt_sub_samp, "Type1", font = 14, labels = list(cex = 1.3), fraction = 0.1)
Venn_turt 

##Shared ASVs egg only
Venn_egg <- ps_venn(Egg_sub_samp, "Type1", font = 14, labels = list(cex = 1.3), fraction = 0.1)
Venn_egg

##Shared ASVs babies through development
sub_samp_babies<-subset_samples(Rare_3000, Type1 %in% c("hatchling", "final turt", "One_year_wild"))
sub_samp_babies  #[ 929 taxa and 29 samples ]
Venn_babies <- ps_venn(sub_samp_babies, "Type1", font = 14, labels = list(cex = 1.3), fraction = 0.1)
Venn_babies

##environment 
sub_samp_env<-subset_samples(Rare_3000, Type1 %in% c("soil", "mid sub", "final sub"))
sub_samp_env #[ 929 taxa and 29 samples ]
Venn_env <- ps_venn(sub_samp_env, "Type1", font = 14, labels = list(cex = 1.3), fraction = 0.1)
Venn_env



#### HEATMAP ####
#see Wang et al 2020; package pheatmap; heatmap with hierarchical clustering and dendrogram on axes
# https://cran.r-project.org/web/packages/pheatmap/pheatmap.pdf
# https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/
install.packages(pheatmap)
library(pheatmap)
# install DESeq if necessary
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# install and load package
BiocManager::install("DESeq")
library("DESeq") #Error in library("DESeq") : there is no package called ‘DESeq’

turt_data.mat <- turt_data@otu_table
turt_data.mat <- turt_data.mat@.Data
str(turt_data.mat)

pheatmap(turt_data, color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                              "RdYlBu")))(100), kmeans_k = NA, breaks = NA, border_color = "grey60",
         cellwidth = NA, cellheight = NA, scale = "none", cluster_rows = TRUE,
         cluster_cols = TRUE, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", clustering_method = "complete",
         clustering_callback = identity2, cutree_rows = NA, cutree_cols = NA,
         treeheight_row = ifelse((class(cluster_rows) == "hclust") || cluster_rows,
                                 50, 0), treeheight_col = ifelse((class(cluster_cols) == "hclust") ||
                                                                   cluster_cols, 50, 0), legend = TRUE, legend_breaks = NA,
         legend_labels = NA, annotation_row = NA, annotation_col = NA,
         annotation = NA, annotation_colors = NA, annotation_legend = TRUE,
         annotation_names_row = TRUE, annotation_names_col = TRUE,
         drop_levels = TRUE, show_rownames = T, show_colnames = T, main = NA,
         fontsize = 10, fontsize_row = fontsize, fontsize_col = fontsize,
         angle_col = c("270", "0", "45", "90", "315"), display_numbers = F,
         number_format = "%.2f", number_color = "grey30", fontsize_number = 0.8
         * fontsize, gaps_row = NULL, gaps_col = NULL, labels_row = NULL,
         labels_col = NULL, filename = NA, width = NA, height = NA,
         silent = FALSE, na_col = "#DDDDDD")
#Error in as.vector(data) : no method for coercing this S4 class to a vector

#https://joey711.github.io/phyloseq/plot_heatmap-examples.html for heatmap without dendrogram hiearchical clustering.
library("phyloseq")
library("ggplot2")

genera_turt <- subset_taxa(Rare_3000, level ="Genus")
genera_turt <- prune_taxa(names(sort(taxa_sums(genera_turt),TRUE)[1:300]), genera_turt)
plot_heatmap(genera_turt, sample.label="Type1")

heatmap(otu_table(Rare_3000))

###https://www.r-bloggers.com/2013/02/from-otu-table-to-heatmap/
#otu_table() is a phyloseq function which extract the OTU table from the phyloseq object
ASVs<-otu_table(t(turt_data), taxa_are_rows = TRUE) #[930 taxa and 85 samples]
str(ASVs)

taxa.names <- ASVs$taxonomy
dat2 <- as.matrix(ASV)
                  
