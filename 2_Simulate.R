rm(list=ls())
library(matrixStats)
library(invgamma)
setwd("C:/users/arenasd/Desktop/MPS_112421")
source("f_Simulation.R") #R likes "/"
source("functions.R") #R likes "/"

#========================================================
#Parameters
#========================================================
#--------------------------------------------------------
#Significance, simulation method, and algorithm
#--------------------------------------------------------
significance <- 0.05
simul_method <- "analytical"   #"analytical" or "boot"
algorithm <- "TPPR"   #"PPR" or "mean_ttest"
#--------------------------------------------------------
#Log choice for t-test, upper_c for PPR
#--------------------------------------------------------
logchoice <- FALSE
upper_c = 1 - 0.95
#--------------------------------------------------------
#Number of simulations and group sizes
#--------------------------------------------------------
simul <- 10000
s_v <- c(4,8,6,13,8,11,30,30) #c(subjsize1, ctrlsize1, subjsize2, ctrlsize2,...)

#========================================================
#Effect sizes. Relevant for analytical-simulations choice
#========================================================
mu_subj <- 1 #Average (across individuals) for subjects
sigma_subj <- 0.1 #stdev (across individuals) for subjects
mu_ctrl <- 0 #Average (across individuals) for controls    
sigma_ctrl <- 0.1 #stdev (across individuals) for controls

#========================================================
#Experimental data for boot simulations
#========================================================
#----------------
#File information
#----------------
dm1 = read.csv("1_Original_data_organized.csv", header=TRUE, stringsAsFactors = FALSE)
#-------------------------------------------
#Choose column data and organism for subject 
#-------------------------------------------
dm1$avg.subj <- dm1$PTLD
dm1$avg.ctrl <- dm1$Not_PTLD
dm2 = dm1[dm1$Accession == "NC_007605.1", ]
#----------------
#Number of probes
#----------------
pr_exp = length(dm1$avg.subj)

#===========================
#Simulate various group sizes
#===========================
for (j in c(1:(length(s_v)/2))){
    n_s = s_v[2*j-1]
    n_c = s_v[2*j]
    print(paste0(n_s, ", ", n_c))
    #-------------------------------------
    #Calculate cutoffs  
    #-------------------------------------
    if (simul_method == "analytical"){
        c_subj_v <- rep(calc_cutoff(n_s, 0.05), simul/100)
        c_ctrl_v <- rep(calc_cutoff(n_c, 0.05), simul/100)
    }
    if (simul_method == "boot"){
        c_subj_v <- sim_cutoff_from_bootstrap(dm1$avg.subj, n_s, pr_exp, upper_c, simul/100)
        c_ctrl_v <- sim_cutoff_from_bootstrap(dm1$avg.ctrl, n_c, pr_exp, upper_c, simul/100)
    }
    #=======================================
    #Simulate various probe numbers
    #=======================================
    results_pr <- NULL
    results_posrate <- NULL
    #1,2,3, 4, 5, 6, 7,8,9,10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100
    for (pr in c(1,2,3, 4, 5, 6, 7,8,9,10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100)){
        #========================================================
        #Simulate s experiments
        #========================================================
        p_histo <- NULL #A histogram of the p values for each simulation
        for (s in c(1:simul)){
            #Simulate subject and control data after averaging probe numbers
            if (simul_method == "analytical"){
                subj <- simulate_analytical(pr, n_s, mu_subj, sigma_subj)
                ctrl <- simulate_analytical(pr, n_c, mu_ctrl, sigma_ctrl)
            }
            if (simul_method == "boot"){
                subj <- bootstrap_n_pr(dm2$avg.subj, n_s, pr)
                ctrl <- bootstrap_n_pr(dm2$avg.ctrl, n_c, pr)
            }
            #-----------------------------
            #Perform algorithm of choice
            #-----------------------------
            #t-test of 
            if (algorithm == "mean_ttest"){
                #Average over probes
                a = rowMeans(subj); b = rowMeans(ctrl)
                #Log transformation if desired
                if (logchoice == TRUE){a = log(a+0.000001); b = log(b+0.000001)}
                #t-test of the probe-averages for two groups of subjects
                pt = t.test(a, b, alternative="greater")    
                p_histo <- append(p_histo, pt$p.value)
            }
            if (algorithm == "TPPR"){
                #Average over subjects
                avg.subj = colMeans(subj); avg.ctrl = colMeans(ctrl)
                c.subj <- sample(c_subj_v, 1)   #Sample one of the cutoffs you calculated earlier 
                c.ctrl <- sample(c_ctrl_v, 1)   #Sample one of the cutoffs you calculated earlier 
                df <- data.frame(avg.subj = avg.subj, avg.ctrl = avg.ctrl, cutoff.subj = c.subj, cutoff.ctrl = c.ctrl)
                rv = use_TPPR(df)[1]    
                p_histo <- append(p_histo, rv)
            }
        }
        #How many simulations were positive?
        posrate = sum(p_histo<significance)/length(p_histo); print(paste0(pr, ": ", posrate))
        #Append to a vector that keeps track for each variable probe number
        results_pr = append(results_pr, pr)
        results_posrate = append(results_posrate, posrate)
    }#Loop for different number of probes
    #-------------------------------------
    #Write to output for each group size
    write.csv(data.frame(results_pr, results_posrate), paste0(n_s, "_", n_c, "_", simul_method,"_", algorithm,"_output.csv"))
}#Loop different group sizes
