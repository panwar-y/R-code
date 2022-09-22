#Script to 1) Make Genomic ranges object and overlapping cytobands with phased variants
#          2) Making a table with AID regions (yes/no)
#          3) Getting no. of PV(s) per AID region and Samples per AID region

#Initial Input: Phased Variants exluding normal samples and having filer:Mutect2 Filter = Pass & clustered events (except germline.) and cytobands_hg19.
#Output 1) Data-set of Phased Variants overlapped with Cytoband and (Yes/No) AID Column
#       2) Table for no. of PV(s) per AID region and Samples per AID region

#Setting Working Directory#################################################################################################################################################
setwd("~/Dropbox (Partners HealthCare)/MCL_share/MCL_clean")

## Loading Libraries ######################################################################################################################################################
library(magrittr)
library(dplyr)
library(data.table)
library(plyr)
library(readxl)
library(GenomicRanges)
library(rio)
library(tidyverse)

## Making Genomic Ranges Object ###############################################################################################################################################
#1st Object: Phased Variants exluding normal samples and having filer:Mutect2 Filter = Pass & clustered events (except germline.) 

Clus_df <- read.delim2("~/Dropbox (Partners HealthCare)/MCL_share/MCL_clean/MCL_phased_variants_clustered_filtered_NONORM.txt", sep=",")
#View(Clus_df)

granges_clustered <- makeGRangesFromDataFrame(Clus_df,
                                              keep.extra.columns=FALSE,
                                              ignore.strand=FALSE,
                                              seqinfo=NULL,
                                              seqnames.field=c("seqnames", "seqname",
                                                               "Chromosome", "chrom",
                                                               "chr", "chromosome_name",
                                                               "seqid"),
                                              start.field="Start_Position",
                                              end.field=c("End_position", "stop"),
                                              strand.field="Strand",
                                              starts.in.df.are.0based=FALSE)

#2nd Object: Cytobands
cyto_df <- read.delim("~/Dropbox (Partners HealthCare)/MCL_share/MCL_clean/cytoBand_hg19.txt")
#View(cyto_df)

granges_cytoBand_hg19 <- makeGRangesFromDataFrame(cyto_df,
                                                  keep.extra.columns=FALSE,
                                                  ignore.strand=FALSE,
                                                  seqinfo=NULL,
                                                  seqnames.field="chr",
                                                  start.field="start",
                                                  end.field="stop",
                                                  strand.field="Strand",
                                                  starts.in.df.are.0based=FALSE)

## Merging cytobands with Clus_df making overlap_gr_count2 ##############################################################################################################
overlap_gr <- as.data.frame(findOverlaps(granges_cytoBand_hg19, granges_clustered))

overlap_gr_count <- plyr::count(overlap_gr$queryHits) #adding x and freq col
View(overlap_gr_count)

overlap_gr_count2 <- merge(overlap_gr_count, cyto_df, by.x="x", by.y=0, all.x=TRUE)
View(overlap_gr_count2)

#Merging the overlap_gr_count2, overlap_gr and Clus_df.
overlap_clustered <- merge(Clus_df, overlap_gr, by.x = 0, by.y="subjectHits", all.x=TRUE)
View(overlap_clustered)

overlap_clustered2 <- merge(overlap_clustered, overlap_gr_count2, by.x= "queryHits", by.y ="x", all.x=TRUE)
View(overlap_clustered2)

overlap_clustered_new <- merge(overlap_gr_count2, overlap_clustered, by.x= "x", by.y ="queryHits", all.x=TRUE)
View(overlap_clustered_new)

## Sub-setting the data.
Final_table <- overlap_clustered2[,c("Tumor_Sample_Barcode", "freq", "cytoband2", "Chromosome", "Start_Position", "End_Position", "start", "stop","t_depth")]


## Making a table with AID regions (yes/no) with Final Table#################################################################################################################
Final_table$chr2 <- paste("chr", Final_table$Chromosome,sep = "")
aid <- as.data.frame(read_excel("../References/NIHMS1714594-supplement-Supplementary_Tables_1-9.xlsx", sheet = "Table S2", skip = 3))
aid_g <- GRanges(aid$Chromosome, IRanges(aid$`Region Start`, aid$`Region End`))

Final_table_pv_g <- GRanges(Final_table$chr2, IRanges(Final_table$Start_Position, Final_table$End_Position))

Final_table_pv_aid_overlap <- as.data.frame(findOverlaps(Final_table_pv_g,aid_g)) #Final_table_pv_aid_overlap = FTPAO
head(Final_table_pv_aid_overlap)


maf2_aid_overlap2 <- merge(Final_table, Final_table_pv_aid_overlap, by.x=0, by.y="queryHits", all.x=TRUE)
maf2_aid_overlap2 <- merge(maf2_aid_overlap2, aid, by.x = "subjectHits", by.y = 0, all.x = T)

maf2_aid_overlap2$covers_aid_region <- "No"
maf2_aid_overlap2$covers_aid_region[maf2_aid_overlap2$subjectHits %in% Final_table_pv_aid_overlap$subjectHits] <- "Yes"
View(maf2_aid_overlap2) #Final table to create dynamic bins.

##  Writing the Final table that will be used to create dynamic bins ####################################################################################################################################################
write.csv(maf2_aid_overlap2,"~/Desktop/maf2_aid_overlap2_NONORM.csv", row.names = FALSE)


## Getting no. of PV(s) per AID region and Samples per AID region ###############################################################################################################
maf2_aid_overlap2 <- maf2_aid_overlap2[!duplicated(maf2_aid_overlap2),] #removing duplicates.
View(maf2_aid_overlap2)
names(maf2_aid_overlap2)[14] <- "Region_Start"
names(maf2_aid_overlap2)[15] <- "Region_End"

foo <- plyr::count(maf2_aid_overlap2$Region_Start)
#head(foo)
names(foo)[1] <- "Region_Start"
#head(foo)
names(foo)[2] <- "No_of_PV"

#creating function to group and summarise data
Count_Samples <- function(df) {
      df %>%
      group_by(Region_Start, Region_End, chr2, Tumor_Sample_Barcode) %>%
      dplyr::summarise(Count = n()) }

new_count_s_PV <- Count_Samples(maf2_aid_overlap2)
new_Count_table <- new_count_s_PV %>% drop_na(Region_Start)
#View(new_Count_table)

Total_Sample_count_per_AID<- plyr::count(new_Count_table, c("Region_Start", "Region_End", "chr2"))
#View(Total_Sample_count_per_AID)

names(Total_Sample_count_per_AID)[4] <- "No_of_Samples"
#View(Total_Sample_count_per_AID)

#Merging both counts
PV_Sample_per_AID <- merge(Total_Sample_count_per_AID, foo, by= "Region_Start", all.x=T)
total_samplecount <- length(unique(maf2_aid_overlap2$Tumor_Sample_Barcode))

PV_Sample_per_AID$AID_cohort_freq <- PV_Sample_per_AID$No_of_Samples/total_samplecount
#View(PV_Sample_per_AID)

#Writing final table for no. of PV(s) per AID region and Samples per AID region
write.csv(PV_Sample_per_AID,"~/Desktop/PV_Sample_per_AID_NONORM.csv", row.names = FALSE)

