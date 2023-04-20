#!/usr/bin/env python3

import argparse
import logging
import pandas as pd
import matplotlib.pyplot as plt
from A_SNP_junction_set.output1 import get_snp_junctions
from B_SNP_donor_acceptor_set.output2 import extract_splice_sites_gtf

def main():
    parser = argparse.ArgumentParser(description='Analyze splice junctions and GWAS SNPs')
    parser.add_argument('junction_file', help='Path to the Leafcutter sQTL junction file')
    parser.add_argument('moloc_snp_file', help='Path to the moloc SNP file')
    parser.add_argument('gtf_file', nargs='?', help='Path to GTF annotation file')

    args = parser.parse_args()

    # set up logging configuration
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s', filename='output.log', filemode='w')

    # Function call for Output 1 tasks
    snp_junctions = get_snp_junctions(args.moloc_snp_file, args.junction_file)

    #Function call for Output 2 tasks
    # 1. Identify spplice sites
    extract_splice_sites_gtf(args.gtf_file)

if __name__ == '__main__':
    main()

# To run this script:
# python3 main.py data/MTCL1/MTCL1_QTL_results.tsv data/MTCL1/MTCL1_moloc_results.tsv




