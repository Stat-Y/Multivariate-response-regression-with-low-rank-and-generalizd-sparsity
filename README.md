# Multivariate response regression with low-rank and generalized sparsity

## Overview

Multiple response regression with low-rank and generalized sparsity, also can be noted as Multi-task learning with low-rank and generalized sparsity is a multi-task learning model which can be used for models assuming various structures for the coefficient matrix, e.g., one can assume low-rank and sparse coeffcient matrix or low-rank and generalized sparse coefficient matrix. Optimization for this model can be done with ADMM algorithm. See the main paper for more details.

## Data

One can analyze Cancer Cell Line Encyclopedia (CCLE) data by the proposed model. CCLE data is made up of 482 cancer cell lines and for each line, there are drug resistance responses for 24 drugs and 18988 gene information. More detailed description for the CCLE data is in the main paper. The source for the CCLE data can be found in https://depmap.org/portal/.

In the following data analysis example, we use screened version of CCLE data. 

1. [Example_Data_Z_eq_I.Rdata]<https://github.com/Stat-Y/CCLE/blob/main/Example_Code_Z_eq_I.R> : screened CCLE data for proposed model with Z=I case.
2. Example_Data_Z_neq_I.Rdata : screened CCLE data for proposed model with Z!=I case.

## Functions


## Data Analysis Examples

## Authors

## Contact
