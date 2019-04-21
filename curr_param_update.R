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


curr_param_update<-function(data, curr_param, n_row, hyper, nBin){

    cat("*** inside curr_param_update ***\n")
  
    ###################
    ### Description:
    ### update of all the parameters one by one 
    ### (first D, then M , then V) for all the clusters
    ###################
  
    # setwd("~/Dropbox/ClusteringDTIonstiefel/DensitySVD_sep1_16/")
    # L = generate_simulated_data_ML(N=100) 
    # load("ML_dataset.Rdata")
    # new_param = cluster_param_update(L$data,L$curr_param,L$clust,c(2,3,1),hyper)
    # D=c(5.4,0.7);S=c(14,15);hyper=NULL;hyper$alpha=1;hyper$beta=0
  
    ### cluster_id_vec contains all the id for those clusters that needed to be udpated (e.g: c(3,1,2))
    # for(cluster_id in cluster_id_vec){
    
    #ID_cluster = which(cluster_assign_vec == cluster_id)
    #N = length(ID_cluster)
    #if(N == 0){
    #  sum_X = data[,,1]*0
    #}else{
    N = dim(data)[[3]]
    X = data
    sum_X = apply(X,c(1,2),sum)
    
  #       t1= Sys.time()
  #       if(N == 1){
  #         sum_X = X
  #       }else{
  #         sum_X = X[,,1]*0
  #         for (i in 1:N){
  #           sum_X = sum_X + X[,,i]
  #         }
  #       }
  #       Sys.time()-t1
    #}
    ####################################
    ############# inference ############
    ####################################
  
    M = curr_param$M
    D = curr_param$D
    V = curr_param$V
    
    #if(N > 5){ #### if cluster has at least one data point; how many points?? 
                #### for practical purpose may be we can take N>5 to avoid the 
                #### case S being close to N
    
    ####################################
    ### old M,D,V parameter update
    ####################################
      
    #### conditional of D|M,V for Gibbs
    S = t(M)%*%sum_X%*%t(V)
    ### for D[1] = d_1
    #print(D)
    #print(S)
  
    #####%%%%%%%%%#########**********************
    #### here according to cluster size parameter beta is changing
    hyper$beta = N*0.004
    #############################################
    print(paste0("inside curr_param_update 1."))
  
    ret_d1 = rdensity_d1_d2(D, S, dimIndex=1, N, nBin=nBin, sampleSize=1, n_row, hyper)
    D[1,1] = ret_d1

	  #print(paste0("nBin:",nBin,":D1:accept_ratio_count:",ret_d1$accept_ratio_count))
    print_debug(paste0("new d1 = ",d1, " n_row = ",n_row),hyper$debug)
    ### for D[2] = d_2
    ret_d2 = rdensity_d1_d2(D, S, dimIndex=2, N, nBin=nBin, sampleSize=1, n_row, hyper)
    D[2,2] = ret_d2
	  #print(paste0("nBin:",nBin,":D2:accept_ratio_count:",ret_d2$accept_ratio_count))
	  #str = sprintf("nbin=%d:D1_acpt_cnt=%d:D2_acpt_cnt=%d",nBin,ret_d1$accept_ratio_count,ret_d2$accept_ratio_count)
    #write(str,file=fName,append=TRUE)
	  print_debug(paste0("new d2 = ",d2, " n_row = ",n_row),hyper$debug)
      
    #### conditional of M|D,V for Gibbs
    M_param = hyper$G + sum_X%*%t(V)%*%(D)
    #print(M_param)
    M = rmf.matrix(M_param) 
    ### TBD: make the max absolute value for columnwise vector positive 
      
    #### this should automatically make V with proper sign
    M = change_sign_M(M)
     
    #### conditional of V|M,D for Gibbs
    V_param = hyper$H + (D)%*%t(M)%*%sum_X
    V = rmf.matrix(V_param) 
    ### TBD: make the max absolute value for columnwise vector positive 
    #} ### if N > 0
    
    #### update the current parameter set
    
    curr_param$M = M
    curr_param$D = D
    curr_param$V = V
    
    curr_param$F = M%*%D%*%V
    print(paste0("update done for the data"))
    
  #}
  
    return(curr_param)
  
}
