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

NR_method_derivative_d1<-function(D,S,dimIndex,N,max_iter,hyper,n_row){

  #print(paste0("D = ",D))
  #print(paste0("S = ",S))
  #print(paste0("N = ",N))
  ###################
  ### Description:
  
  ### this method implements NR method to find out the approximate
  ### mode location from density of d1 by trying to find the zero 
  ### of the gradient of density function by iterative method
  
  ###################
  
  #print(paste0("inside NR method for updating d1"))
  
  #g_d1 = f_prime_x/f_x
  #g_prime_d1 = f_prime2_x/f_x - (f_prime_x/f_x)^2
  #x_new = n_old - g_d1/g_prime_d1

  #dyn.load("./C_code/functions_hyper_2by2_R.so")
 
  ###########################  
  d2 = D[2]
  d1_old = D[1]
  d1_new = d1_old

  alpha = hyper$alpha
  beta = hyper$beta
  	
  ############################
  ### check the derivative of log(f(d1;...)) at D[2]
  D_chk = c(D[2],D[2])
  eigenValues=D_chk^2/4
  dRet = 0.0
  hyper0F1_val = .C("hyper_2by2_R",cc=n_row/2,eigenValues,dRet)[[3]]
  dRet = 0.0
  Y1 = .C("partial_d1_hyper_2by2_R",n_row/2,eigenValues,dRet)[[3]]
  LY1 = (S[dimIndex]-beta) - N*Y1/hyper0F1_val ### Oct26
  if(LY1 < 0){
    print_debug(paste0("error due to mode is less than D[2] !!"),hyper$debug)
    return(D[2])
  }
  
  #x = (1:4000)/10
  #out=d_density_d1_log(y=x,D,S,1,N,hyper,1,n_row)
  #pdf("NR_debug_plot_new.pdf")
  #plot(x,exp(out-max(out)),type='l')
  #dev.off()
  
	############################	
  for (i in  1:max_iter){
  
	  d1_old = d1_new
		
    D = c(d1_old,d2)
    eigenValues=D^2/4

    #print(paste0("iter: ",i," d1_new = ",d1_new)) 
    
  	dRet = 0.0
  	print_debug(paste0("eigenvalue = ",eigenValues," N = ",N," S = ",S," n_row", n_row),hyper$debug)
  	
  	#print(paste0("D[1] = ",D[1]))
  	
  	#print(paste0("D[2] = ",D[2]))
  	if((D[1] >= 100) || (D[2] >= 100) || ((D[1]+D[2]) > 150)){
  	  err_msg(paste0("D's element is too large ",eigenValues))
  	  print(paste0("eigenvalue = ",eigenValues," N = ",N," S = ",S," D[1] = ",D[1]," D[2] = ",D[2]))
  	  x = (1:4000)/10
  	  out=d_density_d1_log(y=x,D,S,1,N,hyper,1,n_row)
  	  pdf("NR_debug_plot_new.pdf")
  	  plot(x,exp(out-max(out)),type='l')
  	  dev.off()
  	  warning("warnig NR_method!!!")
  	  hyper0F1_val = approx_hyper(eigenValues,n_row)
  	}else{
  	  hyper0F1_val = .C("hyper_2by2_R",cc=n_row/2,eigenValues,dRet)[[3]]
  	  #log_f_d1_old = (alpha-1.0)*log(d1_old)+(d1_old*(S[dimIndex]-beta))- N * log(hyper0F1_val)
	    print_debug(paste0("0F1 val = ",hyper0F1_val, " n_row = ", n_row),hyper$debug)
  	}
	  dRet = 0.0
	  #hyper0F1_val_1 = .C("hyper_2by2_R",n_row/2,eigenValues,d)[[3]]
	  #D[dimIndex]=d1_old+0.001;
	  #eigenValues2=D^2/4;
	  
	  #hyper0F1_val_2 = .C("hyper_2by2_R",n_row/2,eigenValues2,dRet)[[3]]
	  #log_f_d1_old_1 = (alpha-1.0)*log(d1_old+0.001)+((d1_old+0.001)*(S[dimIndex]-beta))- N * log(hyper0F1_val_2)
	  
	  #LY1 = (log_f_d1_old_1 - log_f_d1_old)/0.001
	  Y1 = .C("partial_d1_hyper_2by2_R",n_row/2,eigenValues,dRet)[[3]]
	  YY1 = .C("partial_d1_partial_d1_hyper_2by2_R",n_row/2,eigenValues,dRet)[[3]]
	  
	  LY1 = (S[dimIndex]-beta) - N*Y1/hyper0F1_val ### Oct26
	  LYY1 = - N*(YY1/hyper0F1_val - (Y1/hyper0F1_val)*(Y1/hyper0F1_val))
	  
	  #Y1 = (hyper0F1_val_2-hyper0F1_val) / 0.001
	  
	  #D[dimIndex]=d1_old-0.001;
	  #eigenValues2=D^2/4;
	  #hyper0F1_val_3 = .C("hyper_2by2_R",n_row/2,eigenValues2,dRet)[[3]]
	  #log_f_d1_old_3 = (alpha-1.0)*log(d1_old-0.001)+((d1_old-0.001)*(S[dimIndex]-beta))- N * log(hyper0F1_val_3)
	  
	  #LY2 = (log_f_d1_old - log_f_d1_old_3)/0.001
	  
	  #Y2 = -(hyper0F1_val_3-hyper0F1_val) / 0.001
	  
	  #YY1 = (Y1-Y2)/0.002
	  #LYY1 = (LY1-LY2)/0.002
	  
	  #del_d1_log_f_d1_old = (S[dimIndex]-beta) - N*Y1/hyper0F1_val 
    
	  #d1_new = d1_old - (log_f_d1_old/del_d1_log_f_d1_old)

	  d1_new = d1_old - LY1/LYY1
	  
	  if(d1_new < D[2]){
	    d1_new = D[2]
	 }
	  
	  #print(LYY1)
	  print_debug(paste0("d1 new (from NR) ",d1_new),hyper$debug)
	  
  }
  #write.table(val,file="./density1.txt")
  
  return(d1_new)  
  
}


NR_method_function_d1<-function(D,S,dimIndex,N,max_iter,hyper,epsilon = 0.1,n_row){
  
  ###################
  ### Description:
  
  ### this method implements NR method to find out the approximate
  ### mode location from density of d1 by trying to find the zero 
  ### of the gradient of density function by iterative method
  
  ###################
  
  #print(paste0("inside NR function method"))
  
  
  
  #g_d1 = f_prime_x/f_x
  #g_prime_d1 = f_prime2_x/f_x - (f_prime_x/f_x)^2
  #x_new = n_old - g_d1/g_prime_d1
  
  #dyn.load("./C_code/functions_hyper_2by2_R.so")
  
  ###########################  
  d2 = D[2]
  d1_old = D[1]
  d1_new = d1_old
  
  alpha = hyper$alpha
  beta = hyper$beta
  
  ############################	
  for (i in  1:max_iter){
    
    d1_old = d1_new
    
    D = c(d1_old,d2)
    eigenValues=D^2/4;
    
    
    dRet = 0.0
    hyper0F1_val = .C("hyper_2by2_R",cc=n_row/2,eigenValues,dRet)[[3]]
    log_f_d1_old = (alpha-1.0)*log(d1_old)+(d1_old*(S[dimIndex]-beta))- N * log(hyper0F1_val)
    
    #print(log_f_d1_old)
    dRet = 0.0
    
    Y1 = .C("partial_d1_hyper_2by2_R",n_row/2,eigenValues,dRet)[[3]]
    #YY1 = .C("partial_d1_partial_d1_hyper_2by2_R",n_row/2,eigenValues,dRet)[[3]]
    
    LY1 = (S[dimIndex]-beta) - N*Y1/hyper0F1_val
    #LYY1 = - N*(YY1/hyper0F1_val - (Y1/hyper0F1_val)*(Y1/hyper0F1_val))
    
    
    d1_new = d1_old - (log_f_d1_old-log(epsilon))/LY1
    #print(d1_new)
    #print(density1(d1_new,D,S,1,N,1))
  }
  #write.table(val,file="./density1.txt")
  
  return(d1_new)  
  
}

