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


### for p=2
generate_single_cluster_simulated_data_ML<-function(n_row = 3, N=1000, data_dir="."){
    
    #source("bivariate_NR_method_d1_d2.R")
    #source("stiefel_SVD.R")
    ###################
    ### Description:
    
    ### Simulated data generation for fixed number of clusters 
    ### Visualize in some random direction on 2D
    ### Clustering heatmap visualization
    
    ###################
    
    ##################
    library(rstiefel)
    ##################
    
    
    ### hyperparameter generation; first col = <1,0,0,...> ; second col = <0,1,0,...>
    M_hyper = matrix(0, nrow = n_row, ncol = 2)
    M_hyper[1,1] = 1
    mid_row = ceiling(n_row/2)
    M_hyper[2,2] = 1
    #M_hyper[n_row,2] = 1
    #M_hyper[mid_row,2] = 1
    
    
    ### mean matrix generation nClust matrices
    M = matrix(rep(0,n_row*2), nrow = n_row)
    #M = rmf.matrix(M_hyper) 
    M = M_hyper
    
    ### disgonal entries
    D = matrix(rep(0,2*2), nrow = 2)
    s = 0.5*sort(rgamma(c(1,1),c(160,100),c(20,20)),decreasing=TRUE)
    D = diag(s)
    
    ### V matrix generation
    V = matrix(rep(0,2*2), nrow = 2)
    V = matrix(c(1,0,0,1),nrow = 2)
    
    ###### main parameter
    F = matrix(rep(0,n_row*2), nrow = n_row)
    
    ##########################
    ### to keep the sign thing consistent
    ### generate the main parameter F
    M_V_list = change_sign_M_V(M,V)
    M = M_V_list$M
    V = M_V_list$V
    F = M%*%D%*%V
    ##########################
    
    
    ##########################
    data = array(rep(0,n_row*2*N), c(n_row, 2, N)); 
    for(i in 1:N){
        data[,,i] = rmf.matrix(F) 
    }
    
    #print(data)
    L = NULL
    L = list(data=data,curr_param=list(M=M,D=D,V=V,F=F))
    #print(L)
    
    file_name = sprintf("%s/ML_dataset_n_%d_p_2.RData",data_dir,n_row)
    save(L,file=file_name)
    #return(L) 
}

###### some random projection of data on 2D
visualize_data<-function(L){
    
    data = L$data
    r1 = runif(dim(data)[1],-1,1)
    random_direction_vec = r1/sqrt(sum(r1*r1))
    
    N = dim(data)[3]
    proj = matrix(rep(0,2*N), c(2, N)); 
    for(i in 1:N){
        proj[,i] = random_direction_vec%*%data[,,i]   
    }
    arr = c('red','blue','green')
    plot(proj[1,],proj[2,],col=arr[L$clust])
    
}


####### clustering visualization
calc_dist_data<-function(L){
    data = L$data
    N = dim(data)[3]
    
    data1 = array(rep(0,3*2*N), c(3, 2, N));
    nClust = length(unique(L$clust))
    
    #   ID1 = NULL
    #   ID2 = NULL
    #   ID3 = NULL
    #   cnt1 = 0
    #   cnt2 = 0
    #   cnt3 = 0
    #   
    #   
    #   
    #   for(i in 1:N){
    #     if(L$clust[i] == 1){
    #       cnt1 = cnt1+1
    #       ID1[cnt1] = i
    #     }
    #     if(L$clust[i] == 2){
    #       cnt2 = cnt2+1
    #       ID2[cnt2] = i
    #     }
    #     if(L$clust[i] == 3){
    #       cnt3 = cnt3+1
    #       ID3[cnt3] = i
    #     }
    #     
    #   }
    
    id_start = 1
    for(i in 1:nClust){
        cnt1 = sum(L$clust == i)
        id_end =  id_start + cnt1-1
        data1[,,id_start:id_end] = data[,,which(L$clust == i)]
        id_start = id_end+1
    }
    
    #data1[,,(cnt1+1):(cnt1+cnt2)] = data[,,ID2]
    #data1[,,(cnt1+cnt2+1):N] = data[,,ID3]
    
    
    dist_mat = matrix(rep(0.0,N*N),nrow=N,ncol=N)
    
    for(i in 1:N){
        for(j in 1:N){
            
            d = 2.0-sum(diag(t(data1[,,i])%*%data1[,,j]))
            dist_mat[i,j] = d
        }
    }
    
    #save(dist_mat,file="dist_mat.Rdata")
    heatmap(dist_mat,Rowv=NA,Colv=NA)
}



