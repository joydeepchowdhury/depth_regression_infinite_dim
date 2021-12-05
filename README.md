The function **conditionaldepth** returns the conditional depth of a number of response value given a
fixed covariate value based on one of the halfspace depth, the spatial depth and the projection depth.
Its arguments are the following:

  (1) `Response`: A n-by-p matrix, where n is the sample size and p (> 1) is the dimension of the response.
                For projection depth, dimnesion of the response must be either 2 or 3.
                
  (2) `Covariate`: A n-by-q matrix, where n is the sample size and q is either the length of the equispaced
                grid where the functional covariate is recorded, or the dimension of the finite dimensional
                covariate. If it is a vector, it is treated as a n-by-1 matrix.
                
  (3) `fixed.covariate`: A vector whose length is equal to q, the number of columns of Covariate.
  
  (4) `response.value`: A m-by-p matrix, whose each row is a response value. If it a vector of length p,
                m is taken to be 1.
                
  (5) `depth.type`: 1 for halfspace depth, 2 for spatial depth and 3 for projection depth. Default value
                is 1.


The function **conditionaltrimmedmean** returns the conditional 100r% trimmed mean given a fixed covariate
value based on one of the halfspace depth, the spatial depth and the projection depth. Its arguments are
the following:
  (1) Response: A n-by-p matrix, where n is the sample size and p (> 1) is the dimension of the response.
                For projection depth, dimnesion of the response must be either 2 or 3.
  (2) Covariate: A n-by-q matrix, where n is the sample size and q is either the length of the equispaced
                grid where the functional covariate is recorded, or the dimension of the finite dimensional
                covariate. If it is a vector, it is treated as a n-by-1 matrix.
  (3) fixed.covariate: A vector whose length is equal to q, the number of columns of Covariate.
  (4) trimming.proportion: A number between 0 and 1, equal to r in the definition of conditional 100r%
                trimmed mean. Default value is 0.5.
  (5) depth.type: 1 for halfspace depth, 2 for spatial depth and 3 for projection depth. Default value
                is 1.


The function conditionalmedian returns the conditional median given a fixed covariate value based on
one of the halfspace depth, the spatial depth and the projection depth. Its arguments are the following:
  (1) Response: A n-by-p matrix, where n is the sample size and p (> 1) is the dimension of the response.
                For projection depth, dimnesion of the response must be either 2 or 3.
  (2) Covariate: A n-by-q matrix, where n is the sample size and q is either the length of the equispaced
                grid where the functional covariate is recorded, or the dimension of the finite dimensional
                covariate. If it is a vector, it is treated as a n-by-1 matrix.
  (3) fixed.covariate: A vector whose length is equal to q, the number of columns of Covariate.
  (4) depth.type: 1 for halfspace depth, 2 for spatial depth and 3 for projection depth. Default value
                is 1.


The function conditionalspread returns the value of the conditional spread measure Delta( r | x ) given
a fixed covariate value 'x' and a number 'r' with 0 < r < 1, based on one of the halfspace depth, the
spatial depth and the projection depth. Its arguments are the following:
  (1) Response: A n-by-p matrix, where n is the sample size and p (> 1) is the dimension of the response.
                For projection depth, dimnesion of the response must be either 2 or 3.
  (2) Covariate: A n-by-q matrix, where n is the sample size and q is either the length of the equispaced
                grid where the functional covariate is recorded, or the dimension of the finite dimensional
                covariate. If it is a vector, it is treated as a n-by-1 matrix.
  (3) fixed.covariate: A vector whose length is equal to q, the number of columns of Covariate.
  (4) r: A number between 0 and 1, equal to r in the definition of Delta( r | x ). Default value is 0.5.
  (5) depth.type: 1 for halfspace depth, 2 for spatial depth and 3 for projection depth. Default value
                is 1.


The function conditionalskewness returns the value of one of the conditional skewness measures
Psi_1( r_1, r_2 | x ) and Psi_2( r_2 | x ) given a fixed covariate value 'x' and numbers 'r_1' and 'r_2'
with 0 < r_1, r_2 < 1, based on one of the halfspace depth, the spatial depth and the projection depth.
Its arguments are the following:
  (1) Response: A n-by-p matrix, where n is the sample size and p (> 1) is the dimension of the response.
                For projection depth, dimnesion of the response must be either 2 or 3.
  (2) Covariate: A n-by-q matrix, where n is the sample size and q is either the length of the equispaced
                grid where the functional covariate is recorded, or the dimension of the finite dimensional
                covariate. If it is a vector, it is treated as a n-by-1 matrix.
  (3) fixed.covariate: A vector whose length is equal to q, the number of columns of Covariate.
  (4) measure.type: 1 for Psi_1( r_1, r_2 | x ) and 2 for Psi_2( r_1 | x ). Default is 1.
  (5) r_1: A number between 0 and 1, equal to r_1 in the definition of Psi_1( r_1, r_2 | x ).
                Default value is 0.1.
  (6) r_2: A number between 0 and 1, equal to r_2 in the definitions of Psi_1( r_1, r_2 | x ) and
                Psi_2( r_2 | x ). Default value is 0.5.
  (7) depth.type: 1 for halfspace depth, 2 for spatial depth and 3 for projection depth. Default value
                is 1.


The file projmed.rcpp contains some internal functions required in the functions hetertestpvalue
and skewtestpvalue. One has to make a R package from projmed.rcpp to use the functions hetertestpvalue
and skewtestpvalue. The package can be made executing the R file makepackageprojmed.


The function hetertestpvalue returns the p value for the test of heteroscedasticity based on random
permutations and one of the halfspace depth, the spatial depth and the projection depth. The test is
based on Delta( 0.5 | x ). See the chapter (section 4 and section 6) for the description of
Delta( 0.5 | x ). Its arguments are the following:
  (1) Response: A n-by-p matrix, where n is the sample size and p (> 1) is the dimension of the response.
                For projection depth, dimnesion of the response must be either 2 or 3.
  (2) Covariate: A n-by-q matrix, where n is the sample size and q is either the length of the equispaced
                grid where the functional covariate is recorded, or the dimension of the finite dimensional
                covariate. If it is a vector, it is treated as a n-by-1 matrix.
  (3) depth.type: 1 for halfspace depth, 2 for spatial depth and 3 for projection depth. Default value
                is 1.
  (4) randpermnum: Number of random permutations on which the test of heteroscedasticity is based on.
                Default value is 1000.


The function skewtestpvalue returns the p value for the test of conditional skewness based on random
samples and one of the halfspace depth, the spatial depth and the projection depth. Its arguments
are the following:
  (1) Response: A n-by-p matrix, where n is the sample size and p (> 1) is the dimension of the response.
                For projection depth, dimnesion of the response must be either 2 or 3.
  (2) Covariate: A n-by-q matrix, where n is the sample size and q is either the length of the equispaced
                grid where the functional covariate is recorded, or the dimension of the finite dimensional
                covariate. If it is a vector, it is treated as a n-by-1 matrix.
  (3) fixed.covariate: A vector whose length is equal to q, the number of columns of Covariate.
  (4) depth.type: 1 for halfspace depth, 2 for spatial depth and 3 for projection depth. Default value
                is 1.
  (5) pvalue.type: 1 for the p value of the test based on only Psi_1(0.1, 0.9 | x),
                   2 for the p value of the test based on only Psi_2(0.5 | x),
                   3 for both the above p values in a bivariate vector in the above order. See the chapter
                (section 5 and section 6) for the descriptions of Psi_1(0.1, 0.9 | x) and Psi_2(0.5 | x).
  (6) bootstrap.repl.num: Number of random samples on which the test of conditional skewness is based on.
                Default value is 1000.

