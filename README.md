# INK-SQUID
Interpreting novel and known sQTLs underlying isoforms in diseases


## Documentation

The steps for each tasks performed in the tool are broken down into the following.
1. [`A_SNP_junction_set`](https://github.com/sheynkman-lab/sqtl_pipeline/blob/master/A_SNP_junction_set/README.md)
2. [`B_SNP_donor_acceptor_set`](https://github.com/sheynkman-lab/sqtl_pipeline/blob/master/B_SNP_donor_acceptor_set/donor_acceptor_set.md)
3. [`C_overlap_set`](https://github.com/sheynkman-lab/sqtl_pipeline/blob/master/C_overlap_set/overlap_set.md)



## Input Data

- Sample MTCL1 datasets are provided in the `/data` directory. 
- Toy GTF annotation file can be downloaded from Biosurfer Zenodo ([link](https://zenodo.org/record/7004071)) or by clicking [here](https://zenodo.org/record/7004071/files/gencode.v38.toy.gtf?download=1). 

The included MTCL1 example can be run with the following base command:

```
python3 main.py data/MTCL1/MTCL1_QTL_results.tsv data/MTCL1/MTCL1_moloc_results.tsv data/gencode_toy/gencode.v38.toy.gtf
```

