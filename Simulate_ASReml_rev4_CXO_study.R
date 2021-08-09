#######################################################################
# Simulating and fitting longitudinal cluster randomized trial designs
# Created by Rhys Bowden
# Modified by J Kasza to simulate data with a discrete-time-decay correlation structure
# 2021-02-03
#######################################################################
# rev4: for rerunning settings with many errors
library("mvtnorm") #Required to simulate the discrete-time-decay random effects
# with thanks to Ren et al. 2019 "A Simulation Study of Statistical Approaches to Data" for the form of the calls to lme4.
library("lme4")
library("tictoc")
library("asreml")
#setwd("")
options(warn = 0) # careful, this changes the setting for other programs too. options(warn=1) to undo it if the previous setting was warn=1.
set.seed(12) # only changes time effects
design = "cxo" #Changes

# simulation parameters
num_repeats = 1000
starting_seed = 1
individual_level = TRUE # TRUE means fitting model to individual observations, rather than cluster-period cell means.
store_full_data = FALSE # stores all the data used to do the inference. Could use a lot of memory if set to TRUE.
store_each_replicate = TRUE
resuming = TRUE # if you set this to true, also set summary_fname and details_fname below.

# design and model parameters. You can specify any of these (except sigma_eps2) as a list.
T_list = c(2,4,2,6,10) # number of periods
c_list = c(10,25,50) # Number of clusters in total
N_list = c(100,50,10)  # number of individual observations per cluster-period
icc_list =0.05 #c(0.01,0.1,0.2) # intra-cluster correlation
cac_list = c(2-0.8,1,0.8) # cluster auto-correlation, r in Jess' misspecification paper. Values greater than 1 indicate BE model instead of DTD model and will be transformed by 2-value later.
sigma_eps2 = 1 # error variance. The other variance parameters are determined by this and the correlation values above

# pre-specifying the time effects -- time effects beta will be the first TT elements of beta_base
beta_base = 20*(sample.int(3,max(T_list),replace=TRUE)-1) # the exact value of these shouldn't matter

# loading which setting to rerun
#rerun_list = read.table("incomplete.txt",header=F)
#rerun_vec = rerun_list$V1

# output files
if(resuming){ # set these two file names if you are resuming from a previous partially completed run
  summary_fname = '2021-07-07_20-19-02_1.txt'
  if(store_each_replicate){
    details_fname = '2021-07-07_20-19-02_1_details.txt'
  }
  error_fname = '2021-07-07_20-19-02_1_err.txt'
  err_file = file(error_fname,open="at")
  progress = read.table(summary_fname,header=T,sep=",")
  starting_setting = max(progress[,'setting_num'])+1
}else{
  starting_setting = 1
  time_str = format(Sys.time(),"%Y-%m-%d_%H-%M-%S")
  summary_str = '1'
  summary_fname = sprintf('%s_%s.txt',time_str,summary_str)# file name for summary file
  if(store_each_replicate){
    details_fname = sprintf('%s_%s_details.txt',time_str,summary_str)
  }
  error_fname = sprintf('%s_%s_err.txt',time_str,summary_str)# file name for error file
  err_file = file(error_fname,open="at")
  save.image(sprintf('%s_%s_Parameters.RData',time_str,summary_str)) # records the parameters used in R data format
}

# generate each combination of parameters
param_grid = expand.grid(T_list,c_list,N_list,icc_list,cac_list)
names(param_grid) = c('T','C','ni','icc_long','cac_long')
num_param_tuples = dim(param_grid)[1]

param_grid$sigmaE2 = rep_len(sigma_eps2,num_param_tuples) # error (residual) variance
param_grid$dgm = rep('',num_param_tuples)
for(param_row in 1:nrow(param_grid)){
  cac_val = param_grid[param_row,'cac_long']
  if(cac_val==1){
    param_grid[param_row,'dgm'] = 'HH'
  }
  if(cac_val>1){
    param_grid[param_row,'dgm'] = 'BE'
    param_grid[param_row,'cac_long'] = 2-param_grid[param_row,'cac_long']
  }
  if(cac_val<1){
    param_grid[param_row,'dgm'] = 'DTD'
  }
}
# pre-allocate results variables----
AIC_pc1 = array(0,dim=c(num_param_tuples,1))
BIC_pc1 = array(0,dim=c(num_param_tuples,1))
AIC_pc2 = array(0,dim=c(num_param_tuples,1))
BIC_pc2 = array(0,dim=c(num_param_tuples,1))
AIC_pc3 = array(0,dim=c(num_param_tuples,1))
BIC_pc3 = array(0,dim=c(num_param_tuples,1))
warning_HH_m = array(0,dim=c(num_param_tuples,1))
warning_BE_m = array(0,dim=c(num_param_tuples,1))
warning_DTD_m = array(0,dim=c(num_param_tuples,1))
error_HH_m = array(0,dim=c(num_param_tuples,1))
error_BE_m = array(0,dim=c(num_param_tuples,1))
error_DTD_m = array(0,dim=c(num_param_tuples,1))
na_HH_m = array(0,dim=c(num_param_tuples,1))
na_BE_m = array(0,dim=c(num_param_tuples,1))
na_DTD_m = array(0,dim=c(num_param_tuples,1))

# ----
if(store_full_data){
  YY1ifull = vector("list",length=num_param_tuples)
  YY1full = vector("list",length=num_param_tuples)
}
if(resuming==FALSE){
  summary_line_df = data.frame(
    setting_num = numeric(),
    dgm = character(),
    TT = numeric(),
    C = numeric(), #  Number of clusters in total
    ni = numeric(),
    icc = numeric(),
    cac = numeric(),
    AIC_pc1 = numeric(),
    BIC_pc1 = numeric(),
    AIC_pc2 = numeric(),
    BIC_pc2 = numeric(),
    AIC_pc3 = numeric(),
    BIC_pc3 = numeric(),
    warning_HH_m = numeric(),
    warning_BE_m = numeric(),
    warning_DTD_m = numeric(),
    error_HH_m = numeric(),
    error_BE_m = numeric(),
    error_DTD_m = numeric(),
    na_HH_m = numeric(),
    na_BE_m = numeric(),
    na_DTD_m = numeric(),
    times_elapsed = numeric()
  )
  write.table( summary_line_df,  
               file=summary_fname, 
               append = T, 
               sep=',', 
               row.names=F, 
               col.names=T )
  if(store_each_replicate){
    setting_full_df = data.frame(
      rep = numeric(),
      setting = numeric(),
      AIC_bin1 = numeric(),
      BIC_bin1 = numeric(),
      AIC_bin2 = numeric(),
      BIC_bin2 = numeric(),
      AIC_bin3 = numeric(),
      BIC_bin3 = numeric(),
      status_HH = numeric(),
      status_BE = numeric(),
      status_DTD = numeric()
    )
    write.table( setting_full_df,  
                 file=details_fname, 
                 append = T, 
                 sep=',', 
                 row.names=F, 
                 col.names=T )
  }
}
# starting_setting=180
 times_elapsed = list()
for (param_ind in starting_setting:num_param_tuples) {
  if(!any(param_ind==rerun_vec)){
    next
  }
  set.seed(param_ind+starting_seed-1)
  writeLines(sprintf('Fitting setting %s',param_ind))
  writeLines('******************............')
  tic()
  # preallocate results variables for this loop
  status_HH = array(0,dim=c(1,num_repeats))
  status_BE = array(0,dim=c(1,num_repeats))
  status_DTD= array(0,dim=c(1,num_repeats))
  AIC_HH = array(0,dim=c(1,num_repeats))
  BIC_HH = array(0,dim=c(1,num_repeats))
  AIC_BE= array(0,dim=c(1,num_repeats))
  BIC_BE = array(0,dim=c(1,num_repeats))
  AIC_DTD= array(NA,dim=c(1,num_repeats))
  BIC_DTD= array(NA,dim=c(1,num_repeats))
  AIC_bin1 = array(0,dim=c(1,num_repeats))
  BIC_bin1= array(0,dim=c(1,num_repeats))
  AIC_bin2 = array(0,dim=c(1,num_repeats))
  BIC_bin2 = array(0,dim=c(1,num_repeats))
  AIC_bin3 = array(0,dim=c(1,num_repeats))
  BIC_bin3 = array(0,dim=c(1,num_repeats))
  
  # design parameters ----
  Generation_model = param_grid[param_ind,'dgm']
  TT = param_grid[param_ind,'T'] # number of time periods
  #cc = param_grid[param_ind,'c'] # number of clusters per treatment sequence
  cc = param_grid[param_ind,'C'] #  Number of clusters in total
  ni = param_grid[param_ind,'ni'] # number of individuals per cluster
  ###changes
  if(design=="cxo"){
    S = 2 # number of treatment sequences
    treat_mat_s = matrix(c(1,0,0,1),2,2)
    treatM = treat_mat_s %x% matrix(1,nrow=cc,ncol=(TT/2)) # rows of treatM are clusters, rather than treatment sequences.
  }
  C = S*cc # total number of clusters
  if(store_full_data){
    YY1ifull[[param_ind]] = matrix(nrow=num_repeats,ncol=ni*C*TT)
    YY1full[[param_ind]] = matrix(nrow=num_repeats,ncol=ni*C*TT)
  }
  
  # model parameters ----
  sigmaE = sqrt(param_grid[param_ind,'sigmaE2']) # sd of epsilon (errors)
  icc = param_grid[param_ind,'icc_long']
  cac = param_grid[param_ind,'cac_long']
  if(Generation_model == 'DTD'){
    sigmaA = sqrt(icc/(1-icc))*sigmaE
  }
  if(Generation_model =='BE' | Generation_model =='HH'){
    sigmaA = sqrt(cac*icc/(1-icc))*sigmaE # standard deviation of cluster random effect
    sigmaNU = sqrt((1-cac)*icc/(1-icc))*sigmaE # standard deviation of cluster-period random effect in block-exchangeable model
  }
  beta = beta_base[1:TT] # vector of fixed time effects
  theta = 0 # treatment effect
  
  # factors ----
  timeM = t(matrix(1:TT,TT,C))
  clusterM = matrix(1:C,C,TT)
  clustimeM = matrix(1:(C*TT),C,TT)
  timeV = as.vector(t(timeM)) # needs transpose for asreml order
  clusterV = as.vector(t(clusterM))
  clustimeV = as.vector(t(clustimeM))
  treatV = as.vector(t(treatM))
  treatVi = rep(treatV,each=ni)
  clusterVi= rep(clusterV,each=ni)
  timeVi = rep(timeV,each=ni)
  clustimeVi = rep(clustimeV,each=ni)
  
  if(!individual_level){
    treatf = factor(treatV)
    timef = factor(timeV)
    clusterf = factor(clusterV)
    clustimef = interaction(clusterf,timef)
  }else{
    treatf = factor(treatVi)
    timef = factor(timeVi)
    clusterf = factor(clusterVi)
    clustimef = interaction(clusterf,timef)
  }
  
  # iteration ---- 
  for(i in 1:num_repeats){
    # ---- sample ----
    epsi = rnorm(C*TT*ni,mean=0,sd=sigmaE) # epsilon (error term)
    eps = epsi # cluster-period cell mean version of epsilon
    dim(eps) = c(ni,C*TT)
    eps = colSums(eps)/ni
    
    betaM = t(matrix(beta,TT,C))
    betaV = as.vector(t(betaM))
    betaVi = rep(betaV,each=ni) # time effect
    thetaV = theta*treatV # treatment effect
    thetaVi = theta*treatVi
    
    if(Generation_model =='DTD'){
      #Need to generate the autocorrelated cluster-period random effects
      #First, construct the covariance matrix of the cluster-period random effects
      Vi <- (sigmaA^2)*(cac^abs(matrix(1:TT,nrow=TT, ncol=TT, byrow=FALSE) - matrix(1:TT,nrow=TT, ncol=TT, byrow=TRUE)))
      
      #Now, generate the correlated random effects for each cluster:
      #this generates a matrix, with one row for each cluster.
      randomeffs <- rmvnorm(C, mean = rep(0,TT), sigma = Vi)
      #cov(randomeffs)
      #cor(randomeffs)
      #Transform to a vector (transposing necessary to ensure the right order):
      randomeffsvec <- as.vector(t(randomeffs))
      
      randeffsi <- rep(randomeffsvec, each=ni) #one for each participant in each cluster-period
      
      YY1 =  eps +randomeffsvec +betaV +thetaV   # cluster-period cell mean level
      YY1i = epsi+randeffsi+betaVi+thetaVi # individual level version
      #YY1 =  eps +randomeffsvec +betaV    # cluster-period cell mean level
      #YY1i = epsi+randeffsi+betaVi # individual level version
    }
    if(Generation_model =='BE' | Generation_model =='HH'){
      za = rnorm(C) # cluster-level random effects (normalized)
      a = za*sigmaA # a (cluster random effects)
      aV = rep(a,each=TT)
      nu = sigmaNU*rnorm(C*TT)
      nuV = nu
      
      aVi = rep(aV,each=ni) # cluster effect
      nuVi = rep(nuV,each=ni) # cluster-period effect
      if(Generation_model =='BE'){
        YY1 =  eps +aV +betaV +thetaV+ nuV  # cluster-period cell mean level #changes
        YY1i = epsi+aVi+betaVi+thetaVi+ nuVi # individual level version #changes
      }
      if(Generation_model =='HH'){
        YY1 =  eps +aV +betaV+thetaV  # cluster-period cell mean level  #changes
        YY1i = epsi+aVi+betaVi+thetaVi # individual level version  #changes
      }
    }
    
    # ----- data allocation and factors -----
    if(individual_level){
      Y1 = YY1i
      Ydf1 = data.frame(YY1i,treatf, timef,clusterf,clustimef)  #changes
    }else{
      Y1 = YY1
      Ydf1 = data.frame(YY1,treatf, timef,clusterf,clustimef)  #changes
    }
    
    # make the data column name the same
    names(Ydf1)[1] = "Y1"
    # ----- store data -----
    if(store_full_data){
      YY1ifull[[param_ind]][i,] = YY1i
      YY1full[[param_ind]][i,] = YY1
    }
    # ---- inference ----
    # fit model
    # wrap in try-catch block in case of errors
    fitY1_DTD = NULL 
    output1 = tryCatch({
      fitY1_DTD = asreml(fixed= Y1 ~ timef+treatf, random= ~id(clusterf):ar1v(timef),
                         residual=~idv(units),data=Ydf1,maxit = 50,trace=F)
      list(fitY1_DTD,0) 
    },
    warning=function(war){
      print(paste("MY_WARNING1:  ",war," param_ind: ",param_ind," rep: ",i))
      writeLines(sprintf("Setting: %d, Rep: %d, Warning: %s",param_ind,i,war),err_file)
      return(list(fitY1_DTD,1))
    },
    error=function(err){
      print(paste("MY_ERROR1:  ",err," param_ind: ",param_ind," rep: ",i))
      writeLines(sprintf("Setting: %d, Rep: %d, Error: %s",param_ind,i,err),err_file)
      return(list(fitY1_DTD,2))
    })
    
    fitY1_DTD = output1[[1]]
    status_DTD[1,i]=output1[[2]]
    
    fitY1_BE = NULL 
    output2 = tryCatch({
      fitY1_BE = asreml(fixed= Y1 ~ timef+treatf, random= ~id(clusterf):corv(timef),
                        residual=~idv(units),data=Ydf1,maxit = 50, trace=F)
      list(fitY1_BE,0)
    },
    warning=function(war){
      print(paste("MY_WARNING2:  ",war," param_ind: ",param_ind," rep: ",i))
      writeLines(sprintf("Setting: %d, Rep: %d, Warning: %s",param_ind,i,war),err_file)
      return(list(fitY1_BE,1))
    },
    error=function(err){
      print(paste("MY_ERROR1:  ",err," param_ind: ",param_ind," rep: ",i))
      writeLines(sprintf("Setting: %d, Rep: %d, Error: %s",param_ind,i,err),err_file)
      return(list(fitY1_BE,2))
    })
    fitY1_BE = output2[[1]]
    status_BE[1,i]=output2[[2]]
    
    fitY1_HH =NULL
    output3 = tryCatch({
      fitY1_HH = asreml(fixed= Y1 ~ timef+treatf, random= ~id(clusterf),
                        residual=~idv(units),data=Ydf1,maxit = 50, trace=F)
      list(fitY1_HH,0)
    },
    warning=function(war){
      print(paste("MY_WARNING1:  ",war," param_ind: ",param_ind," rep: ",i))
      writeLines(sprintf("Setting: %d, Rep: %d, Warning: %s",param_ind,i,war),err_file)
      return(list(fitY1_HH,1))
    },
    error=function(err){
      print(paste("MY_ERROR1:  ",err," param_ind: ",param_ind," rep: ",i))
      writeLines(sprintf("Setting: %d, Rep: %d, Error: %s",param_ind,i,err),err_file)
      return(list(fitY1_HH,2))
    })
    fitY1_HH = output3[[1]]
    status_HH[1,i]=output3[[2]]
    rm(output1)
    rm(output2)
    rm(output3)
    
    if (is.null(fitY1_HH)) {
      AIC_HH[1,i] = NA
      BIC_HH[1,i] = NA
    } else {
      AIC_HH[1,i]=summary(fitY1_HH)$aic[[1]]
      BIC_HH[1,i]=summary(fitY1_HH)$bic[[1]]
    }
    if (is.null(fitY1_BE)) {
      AIC_BE[1,i] = NA
      BIC_BE[1,i] = NA
    } else {
      AIC_BE[1,i]=summary(fitY1_BE)$aic[[1]]
      BIC_BE[1,i]=summary(fitY1_BE)$bic[[1]]
    }
    if (is.null(fitY1_DTD)) {
      AIC_DTD[1,i] = NA
      BIC_DTD[1,i] = NA
    } else {
      AIC_DTD[1,i]=summary(fitY1_DTD)$aic[[1]]
      BIC_DTD[1,i]=summary(fitY1_DTD)$bic[[1]]
    }
    ###
    # treat nonconvergence as an indication that the model is inappropriate
    AICs = c(AIC_HH[1,i],AIC_BE[1,i],AIC_DTD[1,i])
    minAIC = min(AICs,na.rm=TRUE)
    minAICind = which(AICs==minAIC)
    if(length(minAICind)>1){
      minAICind=sample(minAICind,1)
    }
    AICbins = rep(0,3)
    AICbins[minAICind]=1
    AIC_bin1[1,i] = AICbins[1]
    AIC_bin2[1,i] = AICbins[2]
    AIC_bin3[1,i] = AICbins[3]
    
    BICs = c(BIC_HH[1,i],BIC_BE[1,i],BIC_DTD[1,i])
    minBIC = min(BICs,na.rm=TRUE)
    minBICind = which(BICs==minBIC)
    if(length(minBICind)>1){
      minBICind=sample(minBICind,1)
    }
    BICbins = rep(0,3)
    BICbins[minBICind]=1
    BIC_bin1[1,i] = BICbins[1]
    BIC_bin2[1,i] = BICbins[2]
    BIC_bin3[1,i] = BICbins[3]
    
    rm(fitY1_DTD)
    rm(fitY1_BE)
    rm(fitY1_HH)
  } # end of reps
  AIC_pc1 [param_ind,1]=  sum(AIC_bin1[1,],na.rm = TRUE)
  BIC_pc1 [param_ind,1] = sum(BIC_bin1[1,],na.rm = TRUE)
  AIC_pc2 [param_ind,1]=  sum(AIC_bin2[1,],na.rm = TRUE)
  BIC_pc2 [param_ind,1] = sum(BIC_bin2[1,],na.rm = TRUE)
  AIC_pc3 [param_ind,1]=  sum(AIC_bin3[1,],na.rm = TRUE)
  BIC_pc3 [param_ind,1] = sum(BIC_bin3[1,],na.rm = TRUE)
  
  warning_HH_m[param_ind] <- sum(status_HH[1,]==1)
  warning_BE_m[param_ind] <- sum(status_BE[1,]==1)
  warning_DTD_m[param_ind] <- sum(status_DTD[1,]==1)
  error_HH_m[param_ind] <- sum(status_HH[1,]==2)
  error_BE_m[param_ind] <- sum(status_BE[1,]==2)
  error_DTD_m[param_ind] <- sum(status_DTD[1,]==2)
  na_HH_m[param_ind] <- sum(is.na(AIC_bin1[1,]))
  na_BE_m[param_ind] <- sum(is.na(AIC_bin2[1,])) 
  na_DTD_m[param_ind] <-sum(is.na(AIC_bin3[1,]))
  
  times_elapsed[[param_ind]] = toc()
  if(store_each_replicate){
    setting_full_df = data.frame(
      rep = 1:num_repeats,
      setting = rep(param_ind,num_repeats),
      AIC_bin1 = AIC_bin1[1,],
      BIC_bin1 = BIC_bin1[1,],
      AIC_bin2 = AIC_bin2[1,],
      BIC_bin2 = BIC_bin2[1,],
      AIC_bin3 = AIC_bin3[1,],
      BIC_bin3 = BIC_bin3[1,],
      status_HH = status_HH[1,],
      status_BE = status_BE[1,],
      status_DTD = status_DTD[1,]
    )
    write.table( setting_full_df,  
                 file=details_fname, 
                 append = T, 
                 sep=',', 
                 row.names=F, 
                 col.names=F )
  }
  summary_line_df = data.frame(
    setting_num = param_ind,
    dgm = Generation_model,
    TT = TT,
    C = C, #  Number of clusters in total
    ni = ni,
    icc = icc,
    cac = cac,
    AIC_pc1 = AIC_pc1[param_ind,1],
    BIC_pc1 = BIC_pc1[param_ind,1],
    AIC_pc2 = AIC_pc2[param_ind,1],
    BIC_pc2 = BIC_pc2[param_ind,1],
    AIC_pc3 = AIC_pc3[param_ind,1],
    BIC_pc3 = BIC_pc3[param_ind,1],
    warning_HH_m = warning_HH_m[param_ind],
    warning_BE_m = warning_BE_m[param_ind],
    warning_DTD_m = warning_DTD_m[param_ind],
    error_HH_m = error_HH_m[param_ind],
    error_BE_m = error_BE_m[param_ind],
    error_DTD_m = error_DTD_m[param_ind],
    na_HH_m = na_HH_m[param_ind],
    na_BE_m = na_BE_m[param_ind],
    na_DTD_m = na_DTD_m[param_ind],
    times_elapsed = (times_elapsed[[param_ind]]$toc-times_elapsed[[param_ind]]$tic)[["elapsed"]]
  )
  write.table( summary_line_df,  
               file=summary_fname, 
               append = T, 
               sep=',', 
               row.names=F, 
               col.names=F )
}   
# end of parameter choice loop