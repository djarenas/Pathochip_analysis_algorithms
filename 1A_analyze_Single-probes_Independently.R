rm(list=ls())
#-----------------------------
#Libraries
#-----------------------------
library(FSA)
library(matrixStats)
library(data.table)

#==============================
#Set directory
#==============================
setwd("C:/Users/arenasd/Desktop/MPS_112621")

#======================================
#Read PathoChip data
#======================================
dm1 = read.csv("1_Original_data_organized.csv", header=TRUE, stringsAsFactors = FALSE)

#======================================
#Which data do you want to analyze?
#======================================
source("group_colnames.R")   #This file contains the colnames of
subset_subj <- dm1[, SerumiMCD_colnames]
subset_ctrl <- dm1[, SerumControls_colnames]

print(colnames(subset_subj))
print(colnames(subset_ctrl))

#----------------------------------------------------
#Source the separate file that contains the functions
#----------------------------------------------------
source("functions.R") #R likes "/"

#======================================
# Log-transform the data before t-test?\
# Set percentile comparison
#======================================
logchoice = TRUE
upper_c <- 0.05   #Not used in the algorithm, only for reporting amount of probes in each organism above the percentile of all probes

#======================================
#Calculate the percentile cutoff
#Not used for this algorithm only for reporting in the output summary
#======================================
c_subj <- calc_cutoff_from_data(rowMeans(subset_subj), upper_c)
c_ctrl <- calc_cutoff_from_data(rowMeans(subset_ctrl), upper_c)

#========================================
#Perform t-test Method on all probes
#========================================
#---------------------------------------------------------
#Get the data from the subject group, then controls group
#----------------------------------------------------------
np = length(dm1$ProbeID)
print("Number of total probes:"); print(np)
print(colnames(subset_subj))
print(colnames(subset_ctrl))
#---------------------------------------------------------
#For each probe, check there is some signal, then use t test
#----------------------------------------------------------
p_vector <- NULL; avg.subj <- NULL; avg.ctrl <- NULL    
zero_reads = 0     #A counter to keep track of how many organisms had zero intensity accross all probes in the SUBJECT group
#---------------
#For each probe
#---------------
for (i in c(1:np)){
  #--------------------------
  #Check that there is signal
  #--------------------------
  a = as.numeric(subset_subj[i, ])
  b = as.numeric(subset_ctrl[i, ])
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
  avg.subj <- append(avg.subj, mean(a))
  avg.ctrl <- append(avg.ctrl, mean(b))
  p_vector = append(p_vector, p_v)
}#For each probe
#======================================
#Sort by p value, adjust by FDR
#======================================
dm1$avg.subj <- avg.subj
dm1$avg.ctrl <- avg.ctrl
dm1$p_value = p_vector
#-------------------------------------
#Adjust p values
#-------------------------------------
dm1 <- dm1[order(dm1$p_value), ] #Sort the data frame by p values
k_v = c(1:length(dm1$p_value))     #Add a k column: 1,2,3,...
#BH procedure: p --> p*m/k
dm1$p_adj = dm1$p_value*length(dm1$p_value)/k_v
summary <- dm1[, c("ProbeID", "Sequence", "Accession", "Description", "avg.subj", "avg.ctrl", "p_value", "p_adj" )]
write.csv(summary, "results_single-probe_independently.csv")