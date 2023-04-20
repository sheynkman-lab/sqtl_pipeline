Output 1:
Output 1 determines the junction set for each SNP through. It determines, which junctional events are associated with eeach SNP above a given p-value threshold.

Input: 
In order to prepare Output 1, three inputs are required:
(1) The sQTL file generated from leafcutter.
(2) The user-specified set of SNPs included in the leafcutter results.
(3) A p-value threshold. The default is 0.05 however, can be adjusted to a value specified by the user.

Steps: 
The first step is to read in the leafcutter results. 
The scond step is to read in the user-specified SNPs.

Output:
After running, Output 1 will generate a tabular output in CSV format and a graph illustrating the relationship between each SNP and its associated junction set.
