Input
To determine the junction set for each SNP, this process requires the following inputs:
1) Leafcutter sQTL results file: This file provides the sQTLs of interest, which contains the necessary junctions of interest. 
Within the Leafcutter output, there are junction files generated from bam files. These junctions specify the chromosome where 
the junction is located as well as the start and stop position on the chromosome.
2) set of colocalized user-specified SNPs: This file contains the SNP used to run leafcutter and generate the first input file 
above. Included in this is the position of the SNP and which chromosome it is located on. Additionally the original nucleotide 
base as well as the swapped version is included. 
3) alpha (p-value threshold): The default threshold is 0.05 but can be adjusted by the user to a custom value.

Steps
1) For each SNP LeafCutter is used to find the junctions containing the SNP
2) For each SNP the junctions obtained from the prevous step are filtered based on the significance level. This forms the junction set for an individual SNP.

Output
This process produces a dictionary that defines the junction set for each SNP entered by the user. The key for the dictionary is the user entered SNPs while the values contain a list of dictionaries. Within each dictionary is the infromation from the leafcutter generated sQTL file.
