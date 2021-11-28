rm(list=ls())
#-----------------------------
#Libraries
#-----------------------------
library(FSA)

#-----------------------------
#Set your choices
#------------------------------
#Set directory
setwd("C:/Users/arenasd/Desktop/MPS_112821")
#What two groups you want to analyze
#In example: "HHV8", "Not_HHV8", "PTLD", "Not_PTLD", "iMCD", "Not_iMCD", ...
G1 <- "SerumiMCD"
G2 <- "SerumHD"

#----------------------------------------------------
#Source the separate file that contains the functions
#----------------------------------------------------
source("functions.R") #R likes "/"


#---------------------------------------------------------------------------
#Choose Percentile cutoff and significance of the binomial proportion test
#----------------------------------------------------------------------------
upper_c <- 0.05
significance = 0.05

#--------------------------------------------------------------------------------------------------
#Read PathoChip data
#--------------------------------------------------------------------------------------------------
dm1 = read.csv("1_Original_data_organized.csv", header=TRUE, stringsAsFactors = FALSE)

#----------------------------------
#Set data to analyze
#-----------------------------------
dm1$avg.subj <- dm1[[G1]]
dm1$avg.ctrl <- dm1[[G2]]

#-------------------------------------------------------------------------------------
#Calculate the percentile cutoffs using data from ALL probes
#The cutoffs will define the top percentile probes
#-------------------------------------------------------------------------------------
c_subj <- calc_cutoff_from_data(dm1$avg.subj, upper_c)
dm1$cutoff.subj <- rep(c_subj, length(dm1$avg.subj))
c_ctrl <- calc_cutoff_from_data(dm1$avg.ctrl, upper_c)
dm1$cutoff.ctrl <- rep(c_ctrl, length(dm1$avg.ctrl))
print(summary(dm1$cutoff.ctrl))

#==========================================
#Create factors of all possible accessions
#==========================================
accessions <- as.factor(dm1[["Accession"]]); n_organisms = nlevels(accessions)
print("Number of organisms:"); print(n_organisms)
print("Number of probes:"); print(length(dm1$cutoff.ctrl))

#==========================================
#Perform Method TPPR on all accessions
#==========================================
total <- 0         #Count the number of probes whose p-values were below the set significance
summary <- NULL    #Data-frame that will be updated with results for each organism (accession)
pn_vector <- numeric(0)

#---------------------
# For each accession
#---------------------
for (acc in levels(accessions)){
   #---------------------
   #Extract the subset of data for that accession (all possible probes)
   #---------------------
   subset <- acc_extract(dm1, acc)
   aver_subj <- mean(subset$avg.subj)
   aver_ctrl <- mean(subset$avg.ctrl)
   probe_number <- length(subset[["Description"]])
   pn_vector = append(pn_vector, probe_number)
   #---------------------
   #Perform TPPR method
   #---------------------
   res <- use_TPPR(subset)
   p_result <- res[[1]]

   #--------------------------------------------------------------------
   #Print organism information to screen if binomial test is significant
   #--------------------------------------------------------------------
   if (p_result < significance){
     print("------------------------------------------------------------------------------")
     print("Found an organism with a significantly higher number of top-percentile-probes:")
     print("------------------------------------------------------------------------------")
     print(acc)
     print(subset[["Description"]][1])
     print(paste0("Subject probed-averaged intensity: ", format(aver_subj, digits = 2)))
     print(paste0("Controls probed-averaged intensity: ", format(aver_ctrl, digits = 2)))
     print(paste0("Total number of probes: ", format(probe_number, digits = 2)))
     print(paste0("# of probes above top-percentile (subject):", format(res[2], digits = 2)))
     print(paste0("# of probes above top-percentile (controls): ", format(res[3], digits = 2)))
     total = total + 1   #Count it towards positive organisms
   }
   
   #---------------------
   #Update results dataframe
   #---------------------
   summary <- rbind(summary, (data.frame(Accession = subset[["Accession"]][1], subset[["Description"]][1], probe_number, "pos_subj" = res[[2]], "pos_ctrl" = res[[3]], "p_value" = res[[1]], "aver_subj" = mean(subset$avg.subj), "aver_ctrl" = mean(subset$avg.ctrl), "cutoff_subj" = mean(subset$cutoff.subj), "cutoff_control" = mean(subset$cutoff.ctrl))))
}#For each accession

#========================================================
#Sort by p value of binomial test, write results to file
#========================================================
summary <- summary[order(summary$p_value), ]
write.csv(summary, "results_TPPR_analysis.csv")