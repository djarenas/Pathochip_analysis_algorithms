rm(list=ls())
#-----------------------------
#Libraries
#-----------------------------
library(FSA)

#==============================
#Set directory
#==============================
setwd("C:/Users/arenasd/Desktop/MPS_112821")

upper_c <- 0.05
significance = 0.05


#--------------------------------------------------------------------------------------------------
#Read PathoChip data
#--------------------------------------------------------------------------------------------------
dm1 = read.csv("1_Original_data_organized.csv", header=TRUE, stringsAsFactors = FALSE)

#======================================
# Log-transform the data before t-test?\
# Set percentile comparison
#======================================
logchoice = FALSE
upper_c <- 0.05   #Not used in the algorithm, only for reporting amount of probes in each organism above the percentile of all probes

#======================================
#Which data do you want to analyze?
#======================================
source("group_colnames.R")   #This file contains the colnames of
subj_colnames <- SerumiMCD_colnames
ctrl_colnames <- SerumControls_colnames
print(subj_colnames)
print(ctrl_colnames)

#----------------------------------------------------
#Source the separate file that contains the functions
#----------------------------------------------------
source("functions.R") #R likes "/"
 
#--------------------------------------------
#Create factors of all possible accessions
#---------------------------------------------
accessions <- as.factor(dm1[["Accession"]]); n_organisms = nlevels(accessions)
print("Number of organisms:"); print(n_organisms)

#------------------------------------------------------------------
#Calculate the percentile cutoff 
#Not used for t-test, only for reporting in the output summary
#------------------------------------------------------------------
c_subj <- calc_cutoff_from_data(rowMeans(dm1[, subj_colnames]), upper_c)
c_ctrl <- calc_cutoff_from_data(rowMeans(dm1[, ctrl_colnames]), upper_c)

#===========================================================================
#              PERFORM T-TEST ON PROBE-AVERAGE OF EVERY ORGANISM
#===========================================================================
summary <- NULL
zero_reads = 0
for (acc in levels(accessions)){
      #Extract the subset of data for the accession
      subset <- acc_extract(dm1, acc)
      #Calculate number of probes
      probe_number <- length(subset[["Description"]])
      #----------------------------------------------
      #Separate data-frames for subjects and controls
      #-----------------------------------------------
      subset_subj <- subset[, subj_colnames]
      subset_ctrl <- subset[, ctrl_colnames]
      #------------------------------------------------------------------
      #For each individual, average all the probes
      #------------------------------------------------------------------
      a = colMeans(subset_subj)
      b = colMeans(subset_ctrl)
      #------------------------------------------------------------------
      #Calculate how many probes above the percentile cutoff 
      #Not used for t-test, only for reporting in the output summary
      #------------------------------------------------------------------
      pac_subj = sum(rowMeans(subset_subj) > c_subj)
      pac_ctrl = sum(rowMeans(subset_ctrl) > c_ctrl)
      #------------------------------------------------------------------
      #                       Perform one-sided t-test
      #------------------------------------------------------------------
      #--------------------------
      #Check that there is signal
      #--------------------------
      no_signal = sum(a) == 0 & sum(b) == 0
      #------------------------------------
      #If there is no signal, p = 1
      #If there is signal, conduct t-test
      #------------------------------------
      if (no_signal){
         p_v = 1 #p value equals one
         zero_reads = zero_reads + 1
      }
      else{
         #------------------------
         #Log-transform, yes or no
         #------------------------
         if (logchoice == FALSE){
            p_v = t.test(a,b, alternative="greater")$p.value
         }
         if (logchoice == TRUE){
            p_v = t.test(log(a+0.000001),log(b+0.000001), alternative="greater")$p.value
         }
      }
      #------------------------------------------------------------------
      #Update overall results with this accession
      #------------------------------------------------------------------
      summary <- rbind(summary, (data.frame(Accession = subset[["Accession"]][1], subset[["Description"]][1], subset[["Superkingdom"]][1], "avg_subj" = mean(a), "avg_ctrl" = mean(b), "p_value" = p_v, "pac_subj" = pac_subj, "pac_ctrl" = pac_ctrl, "number of probes" = probe_number, "No signal" = no_signal ) )  )   
}
#---------------------------------------
#Sort by p value, adjust by FDR, write results to file
#----------------------------------------
# #Get rid of organisms with no signal
# summary <- summary[!(summary$No.signal == TRUE), ]
#Sort the data frame by p values
summary <- summary[order(summary$p_value), ]
#Now add a k column
summary$k = c(1:length(summary$p_value))
#Adjust p values by BH procedure p*m/k
summary$p_adj = summary$p_value*length(summary$p_value)/summary$k
write.csv(summary, "results_Probe-Averaged_per_Organism.csv")