#Script to Identify Phased Variants(PV) and PVs in AID region
#Input : MAF files obtained after Variant calling and Annotation on Terra Cloud using .wdl scripts
#Output : Data-set of Phased Variants.



#Setting working directory 
setwd("~/Dropbox (Partners HealthCare)/Bioinformatics and Data Science Group/MCL_share/MCL_clean/Samples_maf
      ")
#Loading Library
library(data.table)

#Creating Function for excluding data based on filter
'%!like%' <- function(x,y)!('%like%'(x,y))

#Finding MAF files in Directory
files <- list.files(pattern=".maf")

#Merging All MAF files and filtering 
mergedmaf <- data.frame()
for(i in 1:length(files)){
  m1 <- read.delim2(files[i], skip = 1)
  m1 <- m1[m1$Variant_Type == "SNP",]
  m1 <- m1[m1$FILTER == "PASS" | 
             m1$FILTER %like% "clustered_events",]
  m1 <- m1[m1$FILTER %!like% "germline",]
  
  mergedmaf <- rbind(mergedmaf, m1)
}


#Getting unique samples.
samples <- as.character(unique(mergedmaf$Tumor_Sample_Barcode))

#Identifying Putative Phased Variants.
phasedvars <- data.frame()
for(i in 1:length(samples)){
  s=samples[i]
  smaf <- mergedmaf[mergedmaf$Tumor_Sample_Barcode == s,]
  for(m in 1:nrow(smaf)){
    st=smaf$Start_Position[m]-170
    ed=smaf$Start_Position[m]+170
    ch=smaf$Chromosome[m]
    
    smaf2 <- smaf[smaf$Chromosome == ch &
                    smaf$Start_Position >= st &
                    smaf$Start_Position <= ed,]
    
    if(nrow(smaf2)>1){
      phasedvars <- rbind(phasedvars, smaf2)
    } else next
    
    
  }
}

phasedvars <- phasedvars[!duplicated(phasedvars),] #Removing Duplicates


write.table(phasedvars, "MCL_phased_variants.txt", row.names = F, sep = "\t", quote = F)

#phasedvars_sub <- phasedvars[phasedvars$Tumor_Sample_Barcode %in% samples[1:8],]
#write.table(phasedvars_sub, "Phased_Variants_subset.txt", row.names = F, sep = "\t", quote = F)
