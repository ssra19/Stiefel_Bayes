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
#### preprocessing data and getting intial clusters
####### hierarchical clustering for initialization of parameters
create_hclust<-function(data,K=3){
  
  N = dim(data)[3]

  
  #dist_mat0 = matrix(rep(0.0,N*N),nrow=N,ncol=N)
  #for(i in 1:N){
  #  for(j in 1:N){
  #    d = 2.0-sum(diag(t(data[,,i])%*%data[,,j]))
  #    dist_mat0[i,j] = d
  #  }
  #}
  ### generating pairwise distance matrix
  date()
  vec_len = dim(data)[1]*dim(data)[2]
  tmp_data = matrix(rep(0.0,N*vec_len),ncol=vec_len)
  tmp_data = t(matrix(data,c(vec_len,N)))
  dist_mat = 2.0-(tmp_data%*%t(tmp_data))
  date()
  
  input_hclust_dist = as.dist(dist_mat)
  h = hclust(input_hclust_dist)
  id_arr = cutree(h,K)
  #plot(h)
  #heatmap(dist_mat,Rowv=NA,Colv=NA)
  
  
  ##### initial parameter finding by Mardia book
  cluster_initial_est = NULL
  cluster_initial_est$M = array(rep(0,3*2*K), c(3, 2, K))
  cluster_initial_est$D = array(rep(0,2*2*K), c(2, 2, K))
  cluster_initial_est$V = array(rep(0,2*2*K), c(2, 2, K))
  cluster_initial_est$id_arr = id_arr
  
  for(clust_id in 1:K){
    idx = which(id_arr == clust_id)
    return_all_cluster_initial_est = initial_parameter_est_each_cluster(data,idx)
    cluster_initial_est$M[,,clust_id] = return_all_cluster_initial_est$M
    cluster_initial_est$D[,,clust_id] = return_all_cluster_initial_est$D
    cluster_initial_est$V[,,clust_id] = return_all_cluster_initial_est$V
  }
  
  return(cluster_initial_est)
  
  ## 2-sum(diag(t(L$curr_param$M[,,1])%*%init_est$M[,,1]))
  
}

init_param_true <- function(nClust=3){
  
  M_hyper_arr = array(c(1,0,0,0,1,0,1,0,0,0,-1,0,0,-1,0,0,0,-1), c(3, 2, 3))
  
  ### mean matrix generation nClust matrices
  M = array(rep(0,3*2*nClust), c(3, 2, nClust))
  for(i in 1:nClust){
    M[,,i] = rmf.matrix(M_hyper_arr[,,i]) 
  }
  
  D = array(rep(0,2*2*nClust), c(2, 2, nClust))
  for(i in 1:nClust){
    D[,,i] = matrix(c(5+3*i,0,0,3+3*i))
  }
  
  V = array(rep(0,2*2*nClust), c(2, 2, nClust))
  for(i in 1:nClust){
    V[,,i] = matrix(c(1,0,0,1))
  }
  
  curr_param = NULL
  curr_param$M = M
  curr_param$D = D
  curr_param$V = V
  
  return(curr_param)
}
#####################

initial_parameter_est_each_cluster <-function(data,idx){
  library(expm)
  return_all_cluster_initial_est = NULL
  
  # sum_X = matrix(rep(0.0,3*2),c(3,2))
  X_bar = matrix(rep(0.0,3*2),c(3,2))
  #   L1 = length(idx)
  #   for(i in 1:L1){
  #     sum_X = sum_X + L$data[,,idx[i]]
  #   }
  #   X_bar = sum_X/L1
  #   
  X_bar = apply(data[,,idx],c(1,2),mean)
  
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
  return_all_cluster_initial_est$M = M
  
  V = res$v
  return_all_cluster_initial_est$V = V
  
  #### phi needs to be solved
  ### for small phi, phi = p*g_i ************************** TO BE DONE 
  p = dim(X_bar)[2]
  
  ##if(phi are small){
  #phi = p*g 
  #}
  phi = bivariate_NR_method_d1_d2(p,res$d[1],res$d[2],n_row = 3, max_iter=20)
  
  ##if(phi are small){
  #phi = 2(1-g)
  #}
  print(phi)
  D = diag(phi)
  return_all_cluster_initial_est$D = D
  
  
  return(return_all_cluster_initial_est) ### ***** NOT WORKING WELL
}


###############

init_param_from_MLE <-function(data,K){
  
  cat("*** inside init_param_from_MLE ***\n")
  cluster_initial_est = create_hclust(data,K)
  
  init_param = NULL
  init_param$M = cluster_initial_est$M
  init_param$D = cluster_initial_est$D
  init_param$V = cluster_initial_est$V
  init_param$F = array(rep(0,3*2*K), c(3, 2, K))
  for(i in 1:K){
    init_param$F[,,i] = init_param$M[,,i]%*%init_param$D[,,i]%*%init_param$V[,,i]
  }
  init_param$id_arr = cluster_initial_est$id_arr
  return(init_param)
  
}

################################
initial_parameter_est_from_MLE <-function(data,n_row){
  
  cat("*** inside initial_parameter_est_from_MLE ***\n")
  library(expm)
  return_all_cluster_initial_est = NULL
  
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
  return_all_cluster_initial_est$M = M
  
  V = res$v
  return_all_cluster_initial_est$V = V
  
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
  print(phi)
  D = diag(phi)
  return_all_cluster_initial_est$D = D
  
  return_all_cluster_initial_est$F = M%*%D%*%V
  
  return(return_all_cluster_initial_est) ### ***** NOT WORKING WELL
}


compare_initial_est <-function(L){
  
  K = dim(L$curr_param$M)[3]
  cluster_initial_est = create_hclust(L$data,K)
  for(i in 1:K){
    print(paste0("****************************"))
    print(cbind(L$curr_param$M[,,i],cluster_initial_est$M[,,i]))
    print(paste0("distance of M ",i," = ", 2-sum(diag(t(L$curr_param$M[,,i])%*%cluster_initial_est$M[,,i]))))
    print(paste0("============================="))
    print(cbind(L$curr_param$D[,,i],cluster_initial_est$D[,,i]))
    print(paste0("distance of D ",i," = ", diag((L$curr_param$D[,,i]-cluster_initial_est$D[,,i])) ))
    print(paste0("============================="))
    print(cbind(L$curr_param$V[,,i],cluster_initial_est$V[,,i]))
    print(paste0("distance of V ",i," = ", 2-sum(diag(t(L$curr_param$V[,,i])%*%cluster_initial_est$V[,,i]))))
    
    print(paste0("****************************"))
  }
  for(i in 1:K){
    print(paste0("============= F MATRIX ================"))
    F_true = L$curr_param$M[,,i]%*%L$curr_param$D[,,i]%*%L$curr_param$V[,,i]
    F_est  = cluster_initial_est$M[,,i]%*%cluster_initial_est$D[,,i]%*%cluster_initial_est$V[,,i]
    print(cbind(F_true,F_est))
    print(paste0("distance of F ",i," = ", norm(F_est-F_true)))
  }
  return(cluster_initial_est)
}
