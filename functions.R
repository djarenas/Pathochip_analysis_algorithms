#=============================================
#normalize_columns
#=============================================
normalize_columns <- function(df, string){
  #Check each column
  for (i in colnames(df)){  
    #For column names that have the string 
    if (grepl(string, i)){     
      soo <- mean(df[[i]]) #Find the sum of the column
      #If the sum is more than zero, normalize so that the mean is 1
      if (soo > 0){                    
        df[i] = df[i]/soo
      }
    }
    
  }
  return(df)  
}
#=============================================
#acc_extract
#=============================================
acc_extract <- function(dm, str_acc){
  #Input1: Data table with all data, as data frame. 
  #Input2: Desired accession, as string.
  #Return: Data frame containing only the desired accession
  subset_dm <- dm[dm$Accession %in% str_acc, ]
  return(subset_dm)
}
#=============================================
#calc_cutoff
#=============================================
calc_cutoff <- function(n, t_one_rate){
  #---------------------------------------------
  #Calculate cutoffs based on backgrounds
  #---------------------------------------------
  simul = 300000
  dat <- NULL
  for (j in 1:n){
    bac = (rnorm(simul, mean=0, sd=1))^2
    dat = rbind(dat, bac)
  } 
  #Average over subjects (Same as summing over k square normal distributions)
  dat = colMeans(dat)
  dat = sort(dat)
  hist(dat, breaks = 100, xlab = "Intensity")
  index <- as.integer(simul*(1-t_one_rate))
  return(dat[index])
}

#=============================================
#use_TPPR
#=============================================
use_TPPR <- function(subset){
  #Input1: Data table with only the desired accession's data
  #Must have four columns
  #   One column for the subject-averaged intensity for each probe ("avg.subj"), another one for controls ("avg.ctrl")
  #   One column for the top-percentile cutoff for the subject ("cutoff.subj"), another one for controls ("cutoff.ctrl")
  #Return: A vector [p-value, # of positive subject probes, # of positive control probes]
  
  #How many probes for this accession?
  probes <- length(subset[[1]])
  #How many probes positive for subjects and control?
  subj_pos <- subset["avg.subj"] > subset["cutoff.subj"]
  ctrl_pos <- subset["avg.ctrl"] > subset["cutoff.ctrl"]
  x1 <- sum(subj_pos)  #Number of probes above the top-percentile of all probes: Subjects
  x2 <- sum(ctrl_pos)  ##Number of probes above the top-percentile of all probes: Controls
  #===========================
  #Calculate proportion t-test 
  #===========================
  #Ignore if number of positive probes is the same for subjects and controls
  if(x1 == x2){ 
      result <- result <- c(1, x1, x2) #result vector (p-value, # of (+) probes in subject, # of (+) probes in control)
  }
  else {
    #Calculate only if there are at least two probes and if there is at least one positive probe in subjects
    if (x1 > 0 & probes > 1){  
      #Proportion test with correction
      re <- prop.test(c(x1, x2), c(probes, probes), p = NULL, alternative = "greater", correct = TRUE)
      result <- c(re$p.value, as.integer(x1), as.integer(x2))  #result vector (p-value, # of (+) probes in subject, # of (+) probes in control)
    }
    else {
      result <- c(1.0,as.integer(x1),as.integer(x2))   
    }
  }
  return(result) #result vector (p-value, # of top-percentile-probes in subject, # of top-percentile-probes in control)
}  

#=============================================
#calc_cutoff_from_data
#=============================================
calc_cutoff_from_data <- function(signal_vector, t_one_rate){
  #---------------------------------------------
  #Calculate cutoffs based on backgrounds
  #---------------------------------------------
  #hist(signal_vector, breaks = 100, xlab = "Intensity", xlim = c(0,2))
  signal_vector <- sort(signal_vector)
  index <- as.integer(length(signal_vector)*(1-t_one_rate))
  cutoff <- signal_vector[index]
  return(cutoff)
}

#=============================================
#probe_histogram
#=============================================
probe_histogram <- function(dm){
  probes_vector <- NULL #A vector or probe numbers so that we can do a histogram
  accessions <- as.factor(dm[["Accession"]])
  for (acc in levels(accessions)){
    #Extract the subset of data for that accession (all possible probes)
    subset <- acc_extract(dm, acc)
    probe_number <- length(subset[["Description"]])
    probes_vector <- append(probes_vector, probe_number)
  }
  print("total number of probes:"); print(sum(probes_vector))
  print("mean number of probes:"); print(mean(probes_vector))  
  
  print("number of organisms with only one probe") 
  print(sum(probes_vector == 1))
  print("number of organisms with probes between 1 and 10")
  b <- probes_vector < 11; print(sum(b)-sum(a))
  print("number of organisms with probes between 11 and 20")
  c<- probes_vector < 21; print(sum(c)-sum(b))
  print("number of organisms with probes between 21 and 50")
  d<- probes_vector <= 50; print(sum(d)-sum(c))
  print("number of organisms with probes between 51 and 99")
  e<- probes_vector < 99; print(sum(e)-sum(d))
  print("number of organisms with probes above 100")
  f<- probes_vector > 99; print(sum(f))
  
  r <- hist(probes_vector, breaks = 250)
  
  jpeg("C:/Users/arenasd/Desktop/histogram.jpg", width = 5, height=3.25, units = "in", res = 1200)
  plot(r$breaks[-1], r$counts, log='y', type='h', xlab = "Number of probes", ylab = "Frequency", xlim = c(0,250))  
  #hist(probes_vector, xlab = "Number of probes", main = "", xlim = c(0,250), ylim = c(0,50), breaks = 250)
  dev.off() 
}
#=============================================
#sim_cutoff_from_bootstrap
#=============================================
sim_cutoff_from_bootstrap <- function(exp_signal, n, pr, upper_c, simul_number){
  simul_number = as.integer(simul_number)
  res_v <- NULL
  for (s in c(1:simul_number)){
    #Simulate n subjects for pr probes, and average over subjects
    background <- colMeans(bootstrap_n_pr(exp_signal, n, pr))
    res_v <- append(res_v, calc_cutoff_from_data(background, upper_c))
  }  
  return(res_v)
}