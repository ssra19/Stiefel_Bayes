# Stiefel_Bayes
An example of how to run the code starting from simulated data generation to the posterior inference is given in the file "runAllCode.R"

In order to run the example, do the following. 

1. Compile the C code "functions_hyper_2by2_R.c" with the following command. 
   R CMD SHLIB functions_hyper_2by2_R.c -I/usr/local/include -L/usr/local/lib -lgsl -lgslcblas -lm
   You need to install GSL library in your computer. Here we assume it is installed in /usr/local
2. Generate data with the help of the function written in R
   generate_single_cluster_simulated_data_ML(n_row = 3, N=1000, data_dir="./data")
   This example is for p=2 and 1000 data points are generated each of them is a orthogonal matrix with dim 3x2.
   The data is generated and saved in a file "ML_dataset_n_3_p_2.RData" inside ./data directory.
3. Load the data using 
   load("./data/ML_dataset_n_3_p_2.RData")
   This would populate a R variable L where L$data is the datset and L$curr_param is the true parameter values (M,D,V) with which the data are generated.
4. Using this data, run the posterior inference using the following R command
   data_with_init_with_MCMC_samples = paramsEstimate(L$data, max_iter, output_file, vague_prior, delta)
   One set of values for vague_prior and delta is given in runAllCode.R file. Output is written in a RData file named "output_file" after running "max_iter" MCMC iterations. 
   
  
