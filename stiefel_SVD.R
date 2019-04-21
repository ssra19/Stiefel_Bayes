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
stiefel_SVD<- function(F){
### to make the first row positive for
  res = svd(F)
  res = adjust_SVD_sign(res)
  
  return(res)
}

adjust_SVD_sign <-function(res){
  
  U = res$u
  D = diag(res$d)
  V = res$v
  
  U1 = apply(U,2,'make_max_postive')
  sign_vec = U1[1,]/U[1,]
  V1 = V%*%diag(sign_vec)
  D1 = D
  
  res = NULL
  res$d = diag(D1)
  res$u = U1
  res$v = t(V1)
  return(res)
  
}

make_max_postive <- function(x){
  s = sign(x[which(abs(x)==max(abs(x)))])
  #s = sign(x[1]) ### alternative Chikuse's definition
  x = x*s
}

change_sign_M_V<- function(M,V){
  ### to keep the sign thing consistent
  M1 = apply(M,2,'make_max_postive')
  tmp_sign = dim(M)[1] - apply(M1==M,2,sum)
  id1 = which(tmp_sign==0)
  sign_vec = rep(-1,length(tmp_sign))
  sign_vec[id1] = 1
  V1 = V%*%diag(sign_vec)
 
  M_V_list = list(M=M1,V=V1)
  return(M_V_list)
  
}

change_sign_M <- function(M){
  M1 = apply(M,2,'make_max_postive')
  return(M1)
}