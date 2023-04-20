import argparse
import logging
import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
from helpers import read_leafcutter_sQTL_file, read_gwas_snp_file, read_moloc_sqtl_file
from analysis import get_filtered_phenotypes, get_filtered_sqtl_junctions, print_sqtl_data
from output1 import get_junction_set

def main():
    parser = argparse.ArgumentParser(description='Analyze splice junctions and GWAS SNPs')
    parser.add_argument('junction_file', help='Path to the Leafcutter sQTL junction file')
    parser.add_argument('gwas_snp_file', help='Path to the GWAS SNP file')
    parser.add_argument('moloc_snp_file', help='Path to the moloc SNP file')
    parser.add_argument('--plot', action='store_true', help='Generate a plot of the splice graph')

    args = parser.parse_args()

    # set up logging configuration
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s', filename='output.log', filemode='w')

    #read LeafCutter sQTL file
    leafcutter_sqtl = read_leafcutter_sQTL_file(args.junction_file)
    logging.info(f'Read {len(leafcutter_sqtl)} junctions from {args.junction_file}')

    #read GWAS SNP file
    gwas_snp_data = read_gwas_snp_file(args.gwas_snp_file)
    logging.info(f'Read {len(gwas_snp_data)} SNPs from {args.gwas_snp_file}')

    #read moloc file
    moloc_snp_data = read_moloc_sqtl_file(args.moloc_snp_file)
    logging.info(f'Read {len(moloc_snp_data)} SNPs from {args.moloc_snp_file}')



    
    phenotype_cluster = 'clu_2300'
    marker = 'chr18:8809447:G:C'
    max_pval_nominal = 0.05


if __name__ == '__main__':
    main()

# To run this script:
# python3 main.py data/MTCL1/MTCL1_QTL_results.tsv data/MTCL1/MTCL1_GWAS_results.tsv data/MTCL1/MTCL1_moloc_results.tsv




