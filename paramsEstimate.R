##############################################################################
#
# Stiefel_Bayes :  Conjugate priors and Posterior Inference for the 
#                  matrix Langivin distribution on Stiefel manifolds 
# Copyright (C) <2019>  <Subhajit Sengupta/Subhadip Pal>
# This file is part of Stiefel_Bayes

#  Stiefel_Bayes is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.

#  Stiefel_Bayes is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

#  You should have received a copy of the GNU General Public License
#  along with Stiefel_Bayes If not, see <http://www.gnu.org/licenses/>.

#############################################################################*/


### parametric mixture modeling of Matrix Langevin with fixed clusters
paramsEstimate <- function(data,max_iter=2,output_file,vague_prior=1, delta){
  
  nBin = ceiling(1/delta)
  
  print(date())
  init_run()
  
  N = dim(data)[[3]]
  n_row = dim(data)[[1]]
  
  MCMC_sample = vector("list", max_iter)
 
  data_with_init_with_MCMC_samples = NULL	
  data_with_init_with_MCMC_samples$data = data 
  #MCMC_output_file = sprintf("MCMC_sample_%s_%s_%s.RData",nc,run_id,true_nc)
  #init_param_output_file = sprintf("init_param_MLE_%s_%s_%s.RData",nc,run_id,true_nc)
  ########### hyper parameters #######
  ### need to select empirical prior
#   hyper=NULL
#   hyper$G = (matrix(c(0,0,0,0,0,0),ncol = 2))
#   hyper$H = (matrix(c(0,0,0,0),ncol = 2))
#   hyper$alpha = 1.0
#   hyper$beta = 0.0
#   hyper$dir_alpha = rep(0.0,nc)
#   for (j in 1:nc)
#     hyper$dir_alpha[j] = 1.0
#   
  ####################################
  
  #### initialization ####
  #Rprof(interval = 0.02)
  ###curr_param = init_param_true(nc)
  ### from MLE (Chikuse)
  init_param = initial_parameter_est_from_MLE(data,n_row)
  
  #init_param = perturb_init_param(init_param)
  
  data_with_init_with_MCMC_samples$init_param = init_param
  ### init_param contains initial value for M, D, V and F
  #save(init_param_MLE=init_param,file=init_param_output_file)
  ########################
  print(date())
  ########################
  #vague_prior = 1
  hyper = hyper_selection(n_row,vague_prior,init_param)
  ########################
  
  curr_param = init_param
# Z = curr_param$id_arr
  
#   prob_vec = rep(1,nc)
#   prob_vec = prob_vec/nc
#   for (i in 1:N){
#     r = rmultinom(1,1,prob_vec)
#     Z[i] = which(r==1)
#   }
  
  ### prior for mixture weights
  #pi_vec = rep(1,nc)/nc
  #pi_vec = rdirichlet(1,hyper$dir_alpha) 
  #cluster_id_cnt = sapply(1:nc, function(x) sum(curr_param$id_arr == x))
  #pi_vec = rdirichlet(1,hyper$dir_alpha+cluster_id_cnt)
  
  d_tmp_val = rep(0,N)
  
  for(iter in 1:max_iter){
    ######################
    if(iter%%10 == 0){
      print(paste0("MCMC main iter = ",iter))
    }
    ##### gibbs step #####
#     t1=Sys.time()
#     clust_assign_prob_mat_0 = matrix(rep(0,N*nc),nrow=N)
#     for (i in 1:N){
#       clust_assign_prob_i = rep(0,nc)
#       for (cluster_id in 1:nc){
#         d_tmp_val[i,cluster_id] = dMatrixLangevin(data[,,i],curr_param$M[,,cluster_id],curr_param$D[,,cluster_id],curr_param$V[,,cluster_id])
#         
#         clust_assign_prob_i[cluster_id] = pi_vec[cluster_id]*d_tmp_val[i,cluster_id]
#       }
#       clust_assign_prob_mat_0[i,] = clust_assign_prob_i
#       #clust_assign_prob_i = clust_assign_prob_i/sum(clust_assign_prob_i)
#       Z[i] = sample(1:nc,1,FALSE,clust_assign_prob_i)
#       #print(clust_assign_prob_i)
#     }
#     Sys.time()-t1
#     
    ### compact way to run the above commented code
    t2=Sys.time()
    #clust_assign_prob_mat = matrix(rep(0,N*nc),nrow=N)
    #d_tmp_val = apply(data,3,function(x) dMatrixLangevin(x,curr_param$M,curr_param$D,curr_param$V))
        
    #clust_assign_prob_mat = d_tmp_val*t(matrix(rep(pi_vec,N),ncol=N))
    #Z = apply(clust_assign_prob_mat,1,function(x) sample(1:nc,1,FALSE,x))
    #Sys.time()-t2
  
  
    #cluster_assign_vec = Z
    #write.table(d_tmp_val,file="d_tmp_val.txt")
    #tmp_table = as.matrix(table(Z))
    #cluster_id_cnt = cbind(as.integer(row.names(tmp_table)),tmp_table)
    #cluster_id_cnt = rep(0,length(hyper$dir_alpha))
    
    #for (cluster_id in 1:nc){
    #  cluster_id_cnt[cluster_id] = sum(cluster_assign_vec == cluster_id)
    #}
    ## as an alternative to previous 3 statements
    #cluster_id_cnt = sapply(1:nc, function(x) sum(cluster_assign_vec == x))
      
    #print(hyper$dir_alpha)
    #print(cluster_id_cnt)
    #pi_vec = rdirichlet(1,hyper$dir_alpha+cluster_id_cnt)
    
    #### for all the clusters ####
    #if(iter == 1){
    #  load("ML_dataset.Rdata")
    #  curr_param = L$curr_param
    #  cluster_assign_vec = L$clust
    #}
    curr_param = curr_param_update(data, curr_param, n_row, hyper, nBin) 
    
    print(paste0("curr_param_update"))
    ### need to update with updated cluster_assign_vec
    ##################################
    
    print(paste0("Parameter Estimation: MCMC iteration update done ",iter))
    print(Sys.time()-t2)
    MCMC_sample[[iter]]$curr_param = curr_param
    #MCMC_sample[[iter]]$pi_vec = pi_vec
    
    if(iter%%10==0){
  	  data_with_init_with_MCMC_samples$MCMC_sample = MCMC_sample 
      save(data_with_init_with_MCMC_samples, file = output_file)
    }
  }
  
  
  
  #G = NULL
  #G$curr_param = curr_param
  #G$clust = cluster_assign_vec
  #print(G$curr_param$D[,,1])
  
  #return(G)
  #save(D1_1,D2_1,D1_2,D2_2,D1_3,D2_3,file="D1_D2.Rdata")
  data_with_init_with_MCMC_samples$MCMC_sample = MCMC_sample 
  save(data_with_init_with_MCMC_samples, file = output_file)
  #save("MCMC_sample", file = MCMC_output_file)
  #Rprof(NULL)
  #out=summaryRprof()
  
  
  return(data_with_init_with_MCMC_samples)
  #return(out)
}


dMatrixLangevin <- function(X,M,D,V){
  
  #dyn.load("./C_code/functions_hyper_2by2_R.so")
  
  F = M %*% D %*% V
  D = diag(D)
  eigenValues=D^2/4;
  dRet = 0.0
  
  n_row = dim(M)[1]
  
  ### normalizing constant
  hyper0F1_val = .C("hyper_2by2_R",a=n_row/2,eigenValues,dRet)[[3]]
  
  ### compute density value
  ML_density =exp(sum(diag(t(F)%*%X))) / hyper0F1_val
  
  return(ML_density)
  
}

init_run <-function(){
  source("utility.R")
  load_src_libs()
  
}


initial_parameter_est_from_MLE <-function(data,n_row){
  
  library(expm)
  return_initial_est = NULL
  
  # sum_X = matrix(rep(0.0,3*2),c(3,2))
  X_bar = matrix(rep(0.0,n_row*2),c(n_row,2))
  #   L1 = length(idx)
  #   for(i in 1:L1){
  #     sum_X = sum_X + L$data[,,idx[i]]
  #   }
  #   X_bar = sum_X/L1
  #   
  X_bar = apply(data,c(1,2),mean)
  
  #   #%%%% following Mardia's notation
  #   R_bar = sqrtm(t(X_bar)%*%X_bar)
  #   polar_M = X_bar%*%solve(R_bar)
  #   eig = eigen(R_bar)  ##### R_bar = t(U)*D*U  
  #   U=t(eig$vectors)
  #   g = eig$values ## as a vector
  #   M = polar_M%*%t(U)
  #   
  ### following Chikuse book page 111 (for M and V as well)
  res = stiefel_SVD(X_bar) ### X_bar = M*D*V
  
  M = res$u
  return_initial_est$M = M
  
  V = res$v
  return_initial_est$V = V
  
  #### phi needs to be solved
  ### for small phi, phi = p*g_i ************************** TO BE DONE 
  p = dim(X_bar)[2]
  
  ##if(phi are small){
  #phi = p*g 
  #}
  phi = bivariate_NR_method_d1_d2(p,res$d[1],res$d[2],n_row,max_iter=25)
  
  ##if(phi are small){
  #phi = 2(1-g)
  #}
  #print(phi)
  D = diag(phi)
  return_initial_est$D = D
  
  return_initial_est$F = M%*%D%*%V
  
  return(return_initial_est) 
}


hyper_selection <-function(n_row,vague_prior,init_param){
  
  hyper=NULL
  
  if(vague_prior == 1){
    print(paste0("using vague hyper-parameters"))
    
    hyper$G = (matrix(rep(0,n_row*2),ncol = 2))
    hyper$H = (matrix(c(0,0,0,0),ncol = 2))
    hyper$alpha = 1.0
    hyper$beta = 0.15 ###### note that, min 5 data made a cluster 5*0.01 = 0.05, 
                      ###### for N==1 if S>0.99 then distribution of d1 hv a heavy tail 
  }
  if(vague_prior == 2)
  { ## emperical
    
    print(paste0("using empirical hyper-parameters"))
    
    hyper$G = (matrix(rep(0,n_row*2),ncol = 2))
    hyper$H = (matrix(c(0,0,0,0),ncol = 2))
    hyper$alpha = 1.0
    
    mean_D1_D2 = sum(init_param$D)/(2*nc)
    beta_mean = 1/mean_D1_D2
    
    hyper$beta = beta_mean
    
  }

  if(vague_prior == 3){
    print(paste0("using alpha and beta set to 1 and 0,respectively; uniform improper prior"))
    
    hyper$G = (matrix(rep(0,n_row*2),ncol = 2))
    hyper$H = (matrix(c(0,0,0,0),ncol = 2))
    hyper$alpha = 1.0
    hyper$beta = 0 ###### note that, min 5 data made a cluster 5*0.01 = 0.05, 
    ###### for N==1 if S>0.99 then distribution of d1 hv a heavy tail 
    
  }
  
  hyper$debug = 0
  
  return(hyper)
  
}

perturb_init_param <- function(init_param_0){
  init_param_1 = init_param_0
  n_row = dim(init_param_0$M)[1]
  n_col = dim(init_param_0$M)[2]
  
  for(i in 1:n_row){
    for(j in 1:n_col){
      init_param_1$M[i,j] = rnorm(1,0.4,0.1)
    }
  }
  for(i in 1:n_col){
    init_param_1$D[i,i] = rnorm(1,0.6,0.1)
    for(j in 1:n_col){
      init_param_1$V[i,j] = rnorm(1,0.4,0.1)
    }
  }
  
  init_param_1$F = init_param_1$M%*%init_param_1$D%*%init_param_1$V
  
  return(init_param_1)
}
