#!/usr/bin/env python3

import argparse
from argparse import ArgumentParser, FileType
import logging
import pandas as pd
from helpers import read_leafcutter_sQTL_file, read_moloc_sqtl_file
import matplotlib.pyplot as plt
from A_SNP_junction_set.output1 import get_snp_junctions
from B_SNP_donor_acceptor_set.output2 import extract_splice_sites_gtf, extract_splice_sites_lc, snp_disrupts_splice_sites
from C_overlap_set.output3 import find_set_overlaps

def main():
    parser = argparse.ArgumentParser(description='Analyze splice junctions and GWAS SNPs')
    parser.add_argument('junction_file', help='Path to the Leafcutter sQTL junction file')
    parser.add_argument('moloc_snp_file', help='Path to the moloc SNP file')
    parser.add_argument('gtf_file', nargs='?', type=FileType('r'), help='input GTF file (use "-" for stdin)') 

    args = parser.parse_args()

    # set up logging configuration
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s', filename='output.log', filemode='w')

    # Read the colocalized SNP file
    snp_list = read_moloc_sqtl_file(args.moloc_snp_file)
    logging.info(f'Read {len(snp_list)} SNPs from {args.moloc_snp_file}')

    # Read the LeafCutter sQTL file
    leafcutter_list = read_leafcutter_sQTL_file(args.junction_file)
    logging.info(f'Read {len(leafcutter_list)} junctions from {args.junction_file}')
    
    # from IPython import embed; embed() 
    # Function call for Output 1 tasks
    jx_set = get_snp_junctions(snp_list, leafcutter_list)
    #Function call for Output 2 tasks
    # 1. Identify splice sites from GTF annotation file
    ss_gtf = extract_splice_sites_gtf(args.gtf_file)
    # 2. Identify splice sites from Leafcutter sQTL junction file
    ss_lc = extract_splice_sites_lc(leafcutter_list)
    # 3. Identify if the colocalized SNP disrupts the splice site from GTF and LC
    snp_ss = snp_disrupts_splice_sites(snp_list, leafcutter_list, ss_lc, ss_gtf)

    # Function call for Output 3 task
    overlap_set = find_set_overlaps(jx_set, snp_ss)

if __name__ == '__main__':
    main()

# To run this script:
# python3 main.py data/MTCL1/MTCL1_QTL_results.tsv data/MTCL1/MTCL1_moloc_results.tsv data/gencode_toy/gencode.v38.toy.gtf




