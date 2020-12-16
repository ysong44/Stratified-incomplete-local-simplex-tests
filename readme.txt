There are two folders. One is for single bandwidth and the other for multiple. For each algorithm, it involves two steps.



-------------------------
1. Single_bandwidth

Step 1: specify data and parameters in the GENERATE_input_single.R file, and run the file to generate the input file "INPUT_TO_SILS_Single.RData"

Step 2: run the RUN_SILS_single.R file, and make sure "INPUT_TO_SILS_Single.RData" has been generated in Step 1. The output is the p-value of the test.

------------------------
2. Multiple_bandwidths

Step 1: specify data and parameters in the GENERATE_input_multiple.R file, and run the file to generate the input file "INPUT_TO_SILS_Multiple.RData"

Step 2: run the RUN_SILS_multiple.R file, and make sure "INPUT_TO_SILS_Multiple.RData" has been generated in Step 1. The output is the p-value of the test.

---------------
Note that for d=2,n=1000,b_n=0.5, it takes around 5 mins with 40 computing cores. 

Of course, the smaller the p-value, the stronger evidence against the null (H_0: the regression function is concave)

@misc{song2020stratified,
      title={Stratified incomplete local simplex tests for curvature of nonparametric multiple regression}, 
      author={Yanglei Song and Xiaohui Chen and Kengo Kato},
      year={2020},
      eprint={2003.09091},
      archivePrefix={arXiv},
      primaryClass={math.ST}
}


