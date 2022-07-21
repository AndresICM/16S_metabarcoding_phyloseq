# 16S_metabarcoding_phyloseq
#Analysis of Mothur Processed Reads in RStudio 
This tutorial was adapted from http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html

# Introduction
This tutorial is a continuation of the 16S-rRNA-Metabarcoding-analysis using Mothur. 
Sequences were generated using data from a hydrocarbon bioremediation project. Two treatments were selected for this tutorial, bioaugmentation with *Acinetobacter, Pseudomonas* and *Rhodococcus* strains, and a control. 
Both treatments were inoculated with a high concentration of diesel before the begining of the experiment and were periodically turned over for aeration. Temperature, pH, total petroleum hydrocarbons (TPH) and other physicochemical parameters were monitored. 
The fundamental question of the experiment was to observe the bacterial communities' changes across the experiment, and evaluate how they change while the TPH concentration decreases.

#Download Files
If you already did the "Raw Reads Analysis Tutorial Using Mothur" (https://github.com/AndresICM/16S-rRNA-Metabarcoding-analysis), then you already have two files needed for this tutorial in your "Contis" folder. If so, copy "Contigs/tutorial.an.shared" and "Contigs/tutorial.taxonomy" to your working folder, otherwise download these files from this repository. 
Additionally, there is a file called Map_file.csv. This file is a table with information from each sample. The first column tells us the group that Mothur assigned to each sample, then there is the original name assigned to them; third, the bioremediation strategy of each sample (bioaugmentation or control); fourth, a number assigned to each sample, just to make our lives easier when making graphics. Finally, there is the real time (in weeks) of when the sample was obtained, the Total Petroleum Hydrocarbons (TPH) in ppm, pH and Humidity (%) of each sample.

```
#Sample,Original_name,Treatment,Number,Time,TPH,pH,Humidity
#Group_0,BA3T01,Bioaugmentation,1,0,17830.29863,6.97,8.8
#Group_1,BA3T03,Bioaugmentation,2,2,15965.66405,8.42,7.18
#Group_2,BA3T07,Bioaugmentation,3,6,6851.652868,9.17,6.22
#Group_3,BA3T11,Bioaugmentation,4,10,7646.375894,9.28,10.02
#Group_4,LF3T01,Control,1,0,18478.85447,6.9,5
#Group_5,LF3T03,Control,2,2,19130.27379,6.84,6.64
#Group_6,LF3T07,Control,3,6,13256.93642,6.6,3.64
#Group_7,LF3T11,Control,4,10,9246.770358,6.9,5.13
```
Now we have everything we need to start this tutorial

# Load Libraries
We will need these R packages to process our 16S rRNA metabarcoding data
```
library(ggplot2)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)
library(ggpubr)
```

# Import Files to RStudio

I will set as a working environment a folder called Tutorial

```
setwd("~/Tutorial")
```
Then, we can assign some variables for the files that we will use

```
# Assign variables for imported data
sharedfile = "Contigs/tutorial.an.shared"
taxfile = "Contigs/tutorial.taxonomy"
mapfile = "Map_file.csv"
```
We can now import the map file and the Mothur data as a phyloseq object

```
# Import mothur data
mothur_data <- import_mothur(mothur_shared_file = sharedfile,
                             mothur_constaxonomy_file = taxfile)
                             
# Import sample metadata
map <- read.csv(mapfile)
map <- sample_data(map)
```
Now to merge the map file with Mothur data, we need that the rownames from the map file matches the sample names in the shared and taxonomy files

```
# Assign rownames to be Sample 
rownames(map) <- map$Sample

# Merge mothurdata object with sample metadata
moth_merge <- merge_phyloseq(mothur_data, map)
moth_merge
```
Now we have a phyloseq object called moth_merge. We can take a look at the column names of that file.

```
colnames(tax_table(moth_merge))

#[1] "Rank1" "Rank2" "Rank3" "Rank4" "Rank5" "Rank6"
```
So, we better assign some taxonomik names to that file and inspect how our column names looks like then

```
colnames(tax_table(moth_merge)) <- c("Kingdom", "Phylum", "Class", 
                                     "Order", "Family", "Genus")
                                     
colnames(tax_table(moth_merge))

#[1] "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"
```
Now we will assign a different name to our phyloseq file, and remove mitochondria and chloroplast OTUs from our files, even though we already have donde that using Mothur

```
Valp <- moth_merge %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Family  != "mitochondria" &
      Class   != "Chloroplast"
  )
Valp
```
#


```
```

```
```
