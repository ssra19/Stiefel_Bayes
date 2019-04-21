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


bivariate_NR_method_d1_d2 <- function(p,g1,g2,n_row,max_iter){

 ### input D is a initial value from Mardia's book for small value estimate
 ### d1 = p*g_1; d2 = p*g2  
 
 print(paste0("inside bivariate NR method for updating d1 and d2"))
 
 #dyn.load("./C_code/functions_hyper_2by2_R.so") 

 nrow_D = p 
 ncol_D = p
 ### compute H
 H = matrix(rep(0.0,nrow_D*ncol_D),ncol=ncol_D)
 
 G = rep(0.0,2)
 G[1] = g1
 G[2] = g2
 
 ### D = p*G
 D = rep(0.0,2)
 D_old = rep(0.0,2)
 D_new = rep(0.0,2)
 
 D[1] = p*g1
 D[2] = p*g2
 
 D_old = D
 D_new = D_old
 
 for (i in  1:max_iter){
   
   D_old = D_new
  
   eigenValues=D_old^2/4
   dRet = 0.0
   
   hyper0F1_val = .C("hyper_2by2_R",cc=n_row/2,eigenValues,dRet)[[3]]
   
   Y1 = .C("partial_d1_hyper_2by2_R",n_row/2,eigenValues,dRet)[[3]]
   YY1 = .C("partial_d1_partial_d1_hyper_2by2_R",n_row/2,eigenValues,dRet)[[3]]
   
   YY12 = .C("partial_d1_partial_d2_hyper_2by2_R",n_row/2,eigenValues,dRet)[[3]]
   YY21 = YY12
   
   Y2 = .C("partial_d2_hyper_2by2_R",n_row/2,eigenValues,dRet)[[3]]
   YY2 = .C("partial_d2_partial_d2_hyper_2by2_R",n_row/2,eigenValues,dRet)[[3]]
   
   H[1,1] = YY1/hyper0F1_val - (Y1/hyper0F1_val)*(Y1/hyper0F1_val)
   H[1,2] = YY12/hyper0F1_val - (Y1/hyper0F1_val)*(Y2/hyper0F1_val)
   H[2,1] = H[1,2]
   H[2,2] = YY2/hyper0F1_val - (Y2/hyper0F1_val)*(Y2/hyper0F1_val)
   
   F = rep(0.0,2)
   F[1] = Y1/hyper0F1_val
   F[2] = Y2/hyper0F1_val
   
   H_inv = solve(H)
   
   D_new = D_old - H_inv%*%(F-G)
   
   print(D_new)
 }
 
 return(as.vector(D_new)) 
}