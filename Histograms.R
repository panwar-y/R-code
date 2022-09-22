# Making Histograms.
# In this script, Histogram is created by removing outlier values and setting any frequency above 100 at 100 to correctly visualize histogram.
# To remove outliers, 2 functions have been created named detect_oulier and remove_outlier.

#Loading libraries #########################

library(plyr)
library(ggplot2)

## Making density plot/ Histogram : No. of PV per Sample###########################################################################################################################

boo <- plyr::count(maf2_aid_overlap2$Tumor_Sample_Barcode)
View(boo)

detect_outlier <- function(x) {
  
  # calculate first quantile
  Quantile1 <- quantile(x, probs=.25)
  
  # calculate third quantile
  Quantile3 <- quantile(x, probs=.75)
  
  # calculate inter quantile range
  IQR = Quantile3-Quantile1
  
  # return true or false
  x > Quantile3 + (IQR*1.5) | x < Quantile1 - (IQR*1.5)
}

# create remove outlier function
remove_outlier <- function(dataframe,
                           columns=names(dataframe)) {
  
  # for loop to traverse in columns vector
  for (col in columns) {
    
    # remove observation if it satisfies outlier function
    dataframe <- dataframe[!detect_outlier(dataframe[[col]]), ]
  }
  
  # return dataframe
  print("Remove outliers")
  print(dataframe)
}

new_boo <-remove_outlier(boo,"freq")
new_boo$adj_freq <- paste(new_boo$freq)
head(new_boo)
new_boo$freq2 <- new_boo$freq

new_boo$freq2[new_boo$freq2 >= 100] <- 100
head(new_boo)
View(new_boo)
new_boo$adj_freq <- NULL

#write.csv(boo, file = "~/Desktop/PV_per_sample.csv", row.names = FALSE)

ggplot(data = new_boo, aes(x = freq2), binwidth = 20) +
  geom_histogram()

#Histogram for just  AID region PV Samples

AID_subset <- maf2_aid_overlap2[maf2_aid_overlap2$covers_aid_region == "Yes",]
View(AID_subset)
AID_sample_count <- plyr::count(AID_subset$Tumor_Sample_Barcode)
View(AID_sample_count)

new_AID_sample_count <- remove_outlier(AID_sample_count,"freq")
ggplot(data = new_AID_sample_count, aes(x = freq), binwidth = 20) +
  geom_histogram()