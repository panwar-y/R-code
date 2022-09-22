##Script to process python script: aid_pv_batching.py's output to make table having PV count in each 1000bp bin region, get no. of unique samples and cohort frequency.
#Input : Dynamic Bins obtained from python script
#Output: Data-set with Cohort frequency and No. of samples and phased variants occurring the specific bin.


#Setting up work directory####
setwd("~/Dropbox (Partners HealthCare)/MCL_share/MCL_clean")

#Loading required libraries
library(dplyr)
library(data.table)
library(plyr)
library(readxl)
library(GenomicRanges)

## Taking Input ######
input_file <- read.csv("~/Dropbox (Partners HealthCare)/MCL_share/MCL_clean/NONORM_dynamic_bins.csv")
View(input_file)
input_file <- input_file[,c("chr2", "Region.Start","Region.End","Tumor_Sample_Barcode","cytoband2","Start_Position","End_Position","ClosestGene")]
input_file <- input_file[!duplicated(input_file),] ## Removing Duplicates.

## counting no. of potential PVs per region #########
Input_Variant_count <- plyr::count(input_file,c("chr2", "Region.Start","Region.End","ClosestGene"))
head(Input_Variant_count)
names(Input_Variant_count)[5] <- "Variant_Count"
View(Input_Variant_count)

## Counting no. of unique samples the PVs come from#######
input_file2 <- input_file[,c("chr2", "Region.Start","Region.End","Tumor_Sample_Barcode")]
input_file2 <- input_file2[!duplicated(input_file2),] ## Removing Duplicates.
Input_Sample_count <- plyr::count(input_file2,c("chr2", "Region.Start","Region.End"))
names(Input_Sample_count)[4] <- "Sample_Count"


## Merging Sample count and PV count table just created.##########
#input Sample count and Input 
input_count <- merge(Input_Variant_count,Input_Sample_count,by=c("chr2", "Region.Start","Region.End"),all.x=T)

totalsamplecount <- length(unique(input_file$Tumor_Sample_Barcode))
totalsamplecount
input_count$cohort_freq <- input_count$Sample_Count/totalsamplecount ## Adding & calculating Cohort frequency column

## Taking aid data as input and then overlapping it with the table above to see if the  region is in AID region or not.##########

aid <- as.data.frame(read_excel("~/Dropbox (Partners HealthCare)/MCL_share/References/NIHMS1714594-supplement-Supplementary_Tables_1-9.xlsx", sheet = "Table S2", skip = 3))
aid_g <- GRanges(aid$Chromosome, IRanges(aid$`Region Start`, aid$`Region End`))


input_count_g <- GRanges(input_count$chr2,IRanges(input_count$Region.Start,input_count$Region.End))

aid_input_overlap <- as.data.frame(findOverlaps(input_count_g, aid_g))

input_count$Overlap_aid_region <- "No"
input_count$Overlap_aid_region[aid_input_overlap$queryHits] <- "Yes"
input_count$Bin_size <- (input_count$Region.End - input_count$Region.Start) #adding the bin size.
View(input_count)

## Writing the input_count file #######
write.csv(input_count, "~/Dropbox (Partners HealthCare)/MCL_share/MCL_clean/MCL_PV_dynamic275_bins.csv", row.names = F)

View(input_count[input_count$Overlap_aid_region == "Yes",]) ## Viewing Regions that overlap with AID regions only.


