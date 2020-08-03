## ----setup, include=FALSE----------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = FALSE)
rm(list = ls())
library(tidyverse)
library(dplyr)
library(plyr)
library(readr)
library(tibble)
library(kableExtra)

mydir <- "RHICCA/FOLLOWUP_6MONTHS"




## ----fileprep, include=FALSE-------------------------------------------------------------------------------------------

# read in .txt file of uncertain IDs as list
uncertain <- scan("C:\\Users\\Olivia Swarthout\\Desktop\\uncertain_6months.txt", what = character(), sep = ",")

# read in all .csv file paths from directory
plates = list.files(path = mydir, pattern = "*.csv", full.names=TRUE, recursive = TRUE, 
                    include.dirs = TRUE ) 


#create an empty data frame to hold all merged data
merged <- data.frame("SampleID"=character(), "PlateID" = character(), "CD163" = character(),
                     "ICAM1"=character(), "ILb"=character(), "IL6"=character())

# get rid of "questionable results"
plates <- Filter(function(x) !any(grepl("Questionable", x)), plates)



## ----plates, echo=FALSE------------------------------------------------------------------------------------------------
print(plates)



## ----merge-------------------------------------------------------------------------------------------------------------

 for (x in plates[1:22]){
  current <- read.csv(x)
  
  #converts entire table to type character  
  current[] <- lapply(current, as.character)
  
  #some tables had slight discrepancies in column names so this will account for that by renaming all columns
  if(ncol(current) == 6){
    colnames(current) <- c("SampleID", "PlateID", "CD163", "ICAM1", "ILb", "IL6")
  }
  if(ncol(current) == 7){
    colnames(current) <- c("SampleID", "PlateID", "CD163", "ICAM1", "ILb", "IL6", "metadata")
  }
  

  #a regex to extract the subdirectory name from the full path
  platenum <- str_extract(x, regex("PLATE \\-*\\d+\\.*\\d*", ignore_case = TRUE))
 
  #take just plate number from subdirectory name, store in plate ID column
  platenum <- gsub(platenum, pattern = "PLATE ", replacement = "")
  
  #get the second element of the path (to get timepoint)
  time <- (strsplit(x, "\\/"))[[1]][2]
  
  current$PlateID <- paste(time, platenum, sep = "_")
  
  #add new column containing timepoint
  current$timepoint <- time
  
  merged <- bind_rows(merged, current)
 }


## ----QC----------------------------------------------------------------------------------------------------------------
s <- which(merged$SampleID == "QC 1"| merged$SampleID == "QC 2" | merged$SampleID == "QC1"|merged$SampleID == "QC2")
merged <- merged[-s,]


## ----uncertain---------------------------------------------------------------------------------------------------------
#flag IDs that were in red text (i.e. validity is uncertain)
#add a metadata column if one does not already exist
if(ncol(merged)<8){
  merged$metadata <- NA
}
p <- which(!is.na(match(merged$SampleID, uncertain)) | merged$SampleID == "")
merged$metadata[p] <- "uncertain/missing ID"



## ----onlimit-----------------------------------------------------------------------------------------------------------
t <- grep("(?<!<\\ )2.44140625", merged$ILb, perl = TRUE, value = FALSE)
if(length(t)==0){
  print("No ILb values lie on the lower detection limit")
}else{
  print("The following samples have ILb value at detection limit:")
  merged[t,]
}
u <- grep("(?<!<\\ )2.44140625", merged$IL6, perl = TRUE, value = FALSE)
if(length(u)==0){
  print("No IL6 values lie on the lower detection limit")
}else{
  print("The following samples have IL6 value at detection limit:")
  merged[u,]
}
v <- grep("(?<!<\\ )0.2441406", merged$CD163, perl = TRUE, value = FALSE)
if(length(v)==0){
  print("No CD163 values lie on the lower detection limit")
}else{
  print("The following samples have CD163 value at lower detection limit:")
  merged[v,]
}
x <- grep("(?<!>\\ )1000", merged$CD163, perl = TRUE, value = FALSE)
if(length(x)==0){
  print("No CD163 values lie on the upper detection limit")
}else{
  print("The following samples have CD163 value at upper detection limit:")
  merged[x,]
}
y <- grep("(?<!>\\ )250000", merged$ICAM1, perl = TRUE, value = FALSE)
if(length(y)==0){
  print("No ICAM1 values lie on the upper detection limit")
}else{
  print("The following samples have ICAM1 value at upper detection limit:")
  merged[y,]
}



## ----imputelimits------------------------------------------------------------------------------------------------------

a <-  grepl("<",merged$CD163)
merged$CD163Limit[a] <- "below"
cat(sum(a, na.rm = TRUE), "samples out of", nrow(merged), "have CD163 levels below detection limit")
merged$CD163[a] <- 0.2441406

b <-  grepl(">",merged$CD163)
merged$CD163Limit[b] <- "above"
cat(sum(b, na.rm = TRUE), "samples out of", nrow(merged), "have CD163 levels above detection limit")
merged$CD163[b] <- 1000

c <- grepl(">",merged$ICAM1)
merged$ICAM1Limit[c] <- "above"
cat(sum(c, na.rm = TRUE), "samples out of", nrow(merged), "have ICAM1 levels above detection limit")
merged$ICAM1[c] <- 250000

d <- grepl("<",merged$ILb)
merged$ILblimit[d] <- "below"
cat(sum(d, na.rm = TRUE), "samples out of", nrow(merged), "have ILb levels below detection limit")
merged$ILb[d] <- 2.44140625

e <-  grepl("<",merged$IL6)
merged$IL6Limit[e] <- "below" 
cat(sum(e, na.rm = TRUE), "samples out of", nrow(merged), "have IL6 levels below detection limit")
merged$IL6[e] <- 2.44140625

merged[merged == "N/A"] <- NA



## ----asnumeric---------------------------------------------------------------------------------------------------------
#converts all biomarker values to type numeric
merged$CD163 <- as.numeric(merged$CD163, na.rm = TRUE)
merged$ICAM1 <- as.numeric(merged$ICAM1, na.rm = TRUE)
merged$IL6 <- as.numeric(merged$IL6, na.rm=TRUE)
merged$ILb <- as.numeric(merged$ILb, na.rm=TRUE)



## ----findduplicates----------------------------------------------------------------------------------------------------
#find number of missing sample IDs
missing <- sum(merged$SampleID == "")
cat(missing, "samples have missing ID")

#create table with non-missing IDs that are repeated
merged$SampleID[merged$SampleID == ""] <- NA

#create frequency table for SampleID columns, subset to only frequencies greater than 1
n_occur <- data.frame(table(merged$SampleID))
duplicates <- n_occur[n_occur$Freq > 1,]
duplicates$Var1 <- factor(duplicates$Var1)
  kable(duplicates, caption = "Sample IDs With More Than One Occurrence") %>%
  kable_styling()


## ----removeduplicates--------------------------------------------------------------------------------------------------

#find plate IDs of repeat plates, reconstruct timepoint identifier
repeats <- grep("REPEATS", plates, value = TRUE)
repeatplates  <- str_extract(repeats, regex("PLATE \\-*\\d+\\.*\\d*", ignore_case = TRUE))
repeatplates <- paste(time, repeatplates, sep = "_")
repeatplates <- gsub(repeatplates, pattern = "PLATE ", replacement = "")

#iterate through duplicates
for (item in duplicates$Var1){
  #for each duplicated ID, find all matching rows
  current <- which(merged$SampleID == item)
  #set indicator boolean, iterate through all matches of the id
  isrep <- FALSE
    for (ind in current){
      if (merged$PlateID[ind] %in% repeatplates){
        #indicate true if one of the matches is from a repeats plate
        isrep <- TRUE
      }
    }
    #if one is from the repeats plate, throw out the other ones that are not repeats
    if(isrep == TRUE){
    for (ind in current) {
      if((merged$PlateID[ind] %in% repeatplates) == FALSE)
        cat("sample with ID", merged$SampleID[ind], "found on plate", merged$PlateID[ind], "removed")
        
        merged <- merged[-ind,]
    }
}
}



## ----------------------------------------------------------------------------------------------------------------------
kable(head(merged)) %>%
  kable_styling()


## ----------------------------------------------------------------------------------------------------------------------
filepath <- paste(time, "_MERGED.csv", sep="")
write.csv(merged, filepath, row.names = FALSE)

cat(filepath, "with", nrow(merged), "lines of data created in working directory")


