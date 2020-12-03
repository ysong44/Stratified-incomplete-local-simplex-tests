****inputGen.R:

1. create a folder inputConfig
2. run inputGen.R to generate inputConfig_** files, which will appear under the folder inputConfig. Here, ** is a number that indexes one parameter configuration.


****nonsym_cluster_para_qc_test_incom_jk.R

1. Create a folder tmpFiles/, where the generated output will be stored

2. Take as input the inputConfig/inputConfig_** files generated before

3. Line 77 describes how the error is generated (currently it is a mixture of normals). Change Line 77 to the desired noise distribution.


****Note:
1. The bernoulliSampling.R and util_funs.R are used in the file "nonsym_cluster_para_qc_test_incom_jk.R"

2. We will clean the code later to make it more user friendly. 
