import argparse
import logging
import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
from A_SNP_junction_set.output1 import get_snp_junctions

def main():
    parser = argparse.ArgumentParser(description='Analyze splice junctions and GWAS SNPs')
    parser.add_argument('junction_file', help='Path to the Leafcutter sQTL junction file')
    parser.add_argument('moloc_snp_file', help='Path to the moloc SNP file')
    parser.add_argument('--plot', action='store_true', help='Generate a plot of the splice graph')

    args = parser.parse_args()

    # set up logging configuration
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s', filename='output.log', filemode='w')

    # Function call for Output 1 tasks
    snp_junctions = get_snp_junctions(args.moloc_snp_file, args.junction_file)

if __name__ == '__main__':
    main()

# To run this script:
# python3 main.py data/MTCL1/MTCL1_QTL_results.tsv data/MTCL1/MTCL1_moloc_results.tsv




