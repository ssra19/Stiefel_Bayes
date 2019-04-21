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
rdensity_d1_d2<-function(D,S,dimIndex,N,nBin=10,sampleSize=1,n_row, hyper){
  
  ###################
  ### Description:
  ### this function is used for sampling d1 and d2 (depends on dimIndex = 1 or 2)
  ### for d1: currently approximate mode of the distribution is searched and a 
  ### nearby value is returned by putting a grid between d2 to d2+3*|d2-mode|
  ### for d2: value is returned by putting a grid in between 0 and d1
  ####################
  
  #print(paste0("inside rdensity_d1_d2"))
  ##out=density1(y=(1:8000)/200,D=c(0.7,5.4),S=c(14,15),1,10,1)
  ##plot(exp(out-max(out)))
  #dyn.load("./C_code/functions_hyper_2by2_R.so")
  
  #### y is the supprt of the distribution discretized in a grid
  #### if dimIndex = 1 then we have to sample for d_1 (> d_2)
  #### if dimIndex = 2 then we have to sample for d_2 (< d_1)
  
  if(length(dim(D)) == 2 ){
    D = diag(D)
  }
  if(length(dim(S)) == 2 ){
    S = diag(S)
  }
  if(nBin > 100)
  {
    nBin = 100
  }

  alpha=hyper$alpha;
  beta=hyper$beta;
  
  #print(paste0("D[1] D[2] : ",D[1]," ",D[2]))
  
  if((dimIndex == 1) & ((S[1]-beta) > 0)) #### d_1 should be in d_2 to \infty
  {
      
    print_debug(paste0("before NR ",D),hyper$debug)
    #### S[1] - N < \epsilon how large \epsilon shoud be ?
    bFlag = 1
    d1_near_mode_0 = D[1]
    #if(abs(S[1] - N) < (N*1.0/100)){
    #if(abs(S[1] - N) < 5){
        
      beta_approx_0 = N/2
      beta_approx_1 = N/2+alpha-1
      beta_approx_2 = -(S[1]-beta-N)
      b_approx = beta_approx_2*D[2] - beta_approx_1 - beta_approx_0
      
      d1_near_mode_0 = (-b_approx+sqrt(b_approx*b_approx+4*beta_approx_1*beta_approx_2*D[2]))/(2*beta_approx_2)
      print_debug(paste0("d1_near_mode_0 ",d1_near_mode_0),hyper$debug)
    #}
    if(bFlag == 1){
      d1_near_mode_1 = d1_near_mode_0
      if(d1_near_mode_0 < 60){
        d1_near_mode_1 = NR_method_derivative_d1(D,S,dimIndex,N,max_iter=10,hyper,n_row) ##### replaced by NR method
        print_debug(paste0("d1_near_mode_1 ",d1_near_mode_1," n_row = ",n_row),hyper$debug)
      }
    }
    print(d1_near_mode_0)
    if(!((d1_near_mode_1 >= 100) || ((d1_near_mode_1+D[2]) > 150))){
        
      d_tmp0 = d_density_d1_log(d1_near_mode_0,D,S,dimIndex,N,hyper,bLog=1,n_row)
      d_tmp1 = d_density_d1_log(d1_near_mode_1,D,S,dimIndex,N,hyper,bLog=1,n_row)
      
      print_debug(paste0("d_tmp0 = ",d_tmp0," d_tmp1 = ",d_tmp1),hyper$debug)
      d1_near_mode = d1_near_mode_1
      
      if(d_tmp0 > d_tmp1)
        d1_near_mode = d1_near_mode_0 
    }else{
      d1_near_mode = d1_near_mode_0
    }
    print_debug(paste0("d1_near_mode =  ",d1_near_mode, " n_row = ",n_row),hyper$debug)
    ################
    
    if(d1_near_mode <= D[2]){
      d1_near_mode = D[2]  ##### because support starts from D[2] and it is maximum at D[2]
    }
    
    #diff_from_d2 = abs(D[2]-d1_near_mode) #### how far should we search; need a lower bound
    #pseudoInfty = D[2] + max(5,20*diff_from_d2)   #### totally adhoc
    pseudoInfty = return_pseudoInfty(D,S,dimIndex,N,hyper,d1_near_mode,n_row) ### to be checked
    #print(paste0("1. pseudoInfty = ",pseudoInfty))
    y = seq(D[2],pseudoInfty,1/nBin) 
  }
  
  
  
  #### in the following case when S[1] <= 0 then we don't need NR method
  if((dimIndex == 1) & ((S[1]-beta) <= 0)) #### d_1 should be in d_2 to \infty but as S[1] is negative so there is no mode so just sample near d2 
  {
    pseudoInfty =  return_pseudoInfty(D,S,dimIndex,N,hyper,D[2],n_row)
    #print(paste0("2. pseudoInfty = ",pseudoInfty)) 
    
    #pseudoInfty = D[2] + 5   #### totally adhoc
    y = seq(D[2],pseudoInfty,1/nBin) 
  }
  
  if(dimIndex == 2) #### d_2 should be 0 to d_1
  {
    if(D[1] <= 1/nBin)
    {
      nBin = 100*nBin
    }
    y = seq(1/nBin,D[1],1/nBin) 
  }
  
  val=0*seq(1:length(y));
  
  for (i in  1:length(y)){
    
    x=y[i];
    
    D[dimIndex]=x;
    eigenValues=D^2/4;
    d = 0.0
    bLog = 1 ### hard-coded to avoid problem
    if(bLog==1){
      
      val[i] = d_density_d1_log(x,D,S,dimIndex,N,hyper,bLog,n_row)
      
      #hyper0F1_val = .C("hyper_2by2_R",a=n_row/2,eigenValues,d)[[3]]
      #val[i]=(alpha-1)*log(x)+(x*(S[dimIndex]-beta))- N * log(hyper0F1_val)
    }else{
      #hyper0F1_val = .C("hyper_2by2_R",a=n_row/2,eigenValues,d)[[3]]
      #val[i]=x^(alpha-1)*exp(x*(S[dimIndex]-beta))/(hyper0F1_val)^N
      val[i] = exp(d_density_d1_log(x,D,S,dimIndex,N,hyper,bLog,n_row))
    } 
    
  }
  #print_debug(val,hyper$debug)
  #write.table(val,file="./rdensity_d1_d2.txt")
  prob=exp(val-max(val))
  #if(dimIndex >= 1){
  #  plot(y,prob,type='l',main=paste0(' d',dimIndex)) ### ploting (optional)
  #}
  
  S1=y[which(prob>0)]
  prob1=prob[which(prob>0)]
  #S1
  x=sample(S1, size=sampleSize, replace = TRUE, prob = prob1)
  
  return(x)  
  
}

return_pseudoInfty <- function(D,S,dimIndex,N,hyper,d1_near_mode,n_row){
  
  epsilon = 0.001
  dist_1 = 5 ### make adaptive ?
  MAX_ITER = 6
  
  #print(paste0("************* IN THE FUNCTION ******"))
  #print(paste0("near mode = ",d1_near_mode))
  #diff_from_d2 = abs(D[2]-d1_near_mode)
  #tmp1 = D[2] + max(5,20*diff_from_d2)
  
  
  d_tmp1 = d_density_d1_log(d1_near_mode,D,S,dimIndex,N,hyper,bLog=1,n_row)
  prop_d1_end = d1_near_mode + dist_1
  d_tmp2 = d_density_d1_log(prop_d1_end,D,S,dimIndex,N,hyper,bLog=1,n_row)
  diff_ratio = exp(d_tmp2-d_tmp1)
  
  print_debug(paste0("d1_near_mode = ",d1_near_mode," prop_d1_end = ",prop_d1_end),hyper$debug)
  print_debug(paste0("d_density_d1_log (near_mode) = ", d_tmp1, " d_density_d1_log (right_tail) = ",d_tmp2),hyper$debug)
  cnt_inr = 1
  
  print_debug(paste0("cnt", cnt_inr, " diff_ratio = ",diff_ratio),hyper$debug)
  while((diff_ratio > epsilon) && (cnt_inr < MAX_ITER)){
    prop_d1_end = d1_near_mode + (2*cnt_inr)*dist_1
    d_tmp2 = d_density_d1_log(prop_d1_end,D,S,dimIndex,N,hyper,bLog=1,n_row)
    diff_ratio = exp(d_tmp2-d_tmp1)
    cnt_inr = cnt_inr + 1
    print_debug(paste0("cnt", cnt_inr, " diff_ratio = ",diff_ratio),hyper$debug)
  }
  return(prop_d1_end)
}


d_density_d1_log<-function(y,D,S,dimIndex,N,hyper,bLog,n_row,approx=0){
  
  ##out=d_density_d1_log(y=(1:4000)/2000,D=c(0.7,5.4),S=c(14,15),1,10,1,n_row)
  ##plot(exp(out-max(out)))

  #dyn.load("./C_code/functions_hyper_2by2_R.so")
  
  #print(paste0("d_density_d1_log 1"))
  
  val=0*seq(1:length(y))
  
  for (i in  1:length(y)){
    
    x=y[i]
    if(x > 100){
      approx = 1
    }
    alpha=hyper$alpha
    beta=hyper$beta
    D[dimIndex]=x
    eigenValues=D^2/4
    d = 0.0
    if(bLog==1){
        
      if(approx==1){
        log_hyper0F1_val = approx_hyper_log(eigenValues,n_row)  #### defined in utility.R 
        #print(paste0("******######******** ",hyper0F1_val))
      }else{
        
        hyper0F1_val = .C("hyper_2by2_R",a=n_row/2,eigenValues,d)[[3]]
        log_hyper0F1_val = log(hyper0F1_val)
        
      }
      val[i]=(alpha-1)*log(x)+(x*(S[dimIndex]-beta))- N * log_hyper0F1_val
      
    }else{
      hyper0F1_val = .C("hyper_2by2_R",a=n_row/2,eigenValues,d)[[3]]
      val[i]=x^(alpha-1)*exp(x*(S[dimIndex]-beta))/(hyper0F1_val)^N
    } 
    
  }
  #write.table(val,file="./density1.txt")
    
  return(val)  
  
}



