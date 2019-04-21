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

### initialization
rm(list=ls(all=TRUE))
set.seed(72357)
source("utility.R")
load_src_libs()

### generate simulate data
generate_single_cluster_simulated_data_ML(n_row = 3, N=1000, data_dir="./data")
### data loading 
load("./data/ML_dataset_n_3_p_2.RData")


### set other parameters
max_iter = 100
output_file = sprintf("output.RData")
vague_prior = 1  
delta = 0.01

### run parameter estimation part to generate the output
data_with_init_with_MCMC_samples = paramsEstimate(L$data, max_iter, output_file, vague_prior, delta)


