# 16S_metabarcoding_phyloseq
#Analysis of Mothur Processed Reads in Rstudio 

# Introduction
This tutorial is a continuation of the 16S-rRNA-Metabarcoding-analysis using Mothur. 
Sequences were generated using data from a hydrocarbon bioremediation project. Two treatments were selected for this tutorial, bioaugmentation with *Acinetobacter, Pseudomonas* and *Rhodococcus* strains, and a control. 
Both treatments were inoculated with a high concentration of diesel before the begining of the experiment and were periodically turned over for aeration. Temperature, pH, total petroleum hydrocarbons (TPH) and other physicochemical parameters were monitored. 
The fundamental question of the experiment was to observe the bacterial communities' changes across the experiment, and evaluate how they change while the TPH concentration decreases.

#Download Files
If you already did the "Raw Reads Analysis Tutorial Using Mothur" (https://github.com/AndresICM/16S-rRNA-Metabarcoding-analysis), then you already have two files needed for this tutorial in your "Contis" folder. If so, copy "Contigs/tutorial.an.shared" and "Contigs/tutorial.taxonomy" to your working folder, otherwise download this files from this repository. 
Additionally, there is file called Map_file.csv. This file is a table with information from each sample. The first column tells us the name that Mothur assigned to each sample, then there is the original name assigned to each sample. Third, the treatment of each sample, if bioaugmentation or control, fourth, a number assigned to each sample for easier graphs. Finally, there is the real time in weeks of when the sample was taken, the Total Petroleum Hydrocarbons (TPH) in ppm, pH and Humidity (%) of each sample.

```
Sample,Original_name,Treatment,Number,Time,TPH,pH,Humidity
Group_0,BA3T01,Bioaugmentation,1,0,17830.29863,6.97,8.8
Group_1,BA3T03,Bioaugmentation,2,2,15965.66405,8.42,7.18
Group_2,BA3T07,Bioaugmentation,3,6,6851.652868,9.17,6.22
Group_3,BA3T11,Bioaugmentation,4,10,7646.375894,9.28,10.02
Group_4,LF3T01,Control,1,0,18478.85447,6.9,5
Group_5,LF3T03,Control,2,2,19130.27379,6.84,6.64
Group_6,LF3T07,Control,3,6,13256.93642,6.6,3.64
Group_7,LF3T11,Control,4,10,9246.770358,6.9,5.13
```

```
```

```
```

```
```

```
```

```
```

```
```

```
```

```
```
