#!/usr/bin/env python3

import argparse
from argparse import ArgumentParser, FileType
import logging
import pandas as pd
import matplotlib.pyplot as plt
from A_SNP_junction_set.output1 import get_snp_junctions
from B_SNP_donor_acceptor_set.output2 import extract_splice_sites_gtf, extract_splice_sites_lc

def main():
    parser = argparse.ArgumentParser(description='Analyze splice junctions and GWAS SNPs')
    parser.add_argument('junction_file', help='Path to the Leafcutter sQTL junction file')
    parser.add_argument('moloc_snp_file', help='Path to the moloc SNP file')
    parser.add_argument('gtf_file', nargs='?', type=FileType('r'), help='input GTF file (use "-" for stdin)') 

    args = parser.parse_args()

    # set up logging configuration
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s', filename='output.log', filemode='w')

    # Function call for Output 1 tasks
    leafcutter_list, snp_junctions = get_snp_junctions(args.moloc_snp_file, args.junction_file)

    #Function call for Output 2 tasks
    # 1. Identify splice sites from GTF annotation file
    extract_splice_sites_gtf(args.gtf_file)
    # 2. Identify splice sites from Leafcutter sQTL junction file
    extract_splice_sites_lc(leafcutter_list)
    # 3. Identify if the colocalized SNP disrupts the splice site from GTF and LC
    # snp_disrupts_splice_site()
if __name__ == '__main__':
    main()

# To run this script:
# python3 main.py data/MTCL1/MTCL1_QTL_results.tsv data/MTCL1/MTCL1_moloc_results.tsv 




