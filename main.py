#!/usr/bin/env python3

import argparse
from argparse import ArgumentParser, FileType
from collections import Counter, defaultdict
from sys import stderr, exit
import pandas as pd
import csv
import matplotlib.pyplot as plt

#%% Parsing arguments from user
parser = argparse.ArgumentParser(description='Analyze splice junctions and GWAS SNPs')
parser.add_argument('junction_file', help='Path to the Leafcutter sQTL junction file')
parser.add_argument('moloc_snp_file', help='Path to the moloc SNP file')
parser.add_argument('gtf_file', nargs='?', type=FileType('r'), help='input GTF file (use "-" for stdin)') 
args = parser.parse_args()

#%% Reads a MOLoc sQTL file and returns a list of dictionaries with the data.
sqtl_list = []
with open(args.moloc_snp_file, 'r') as f:
    header = f.readline().strip().split('\t')
    for line in f:
        data = line.strip().split('\t')
        sqtl_dict = {}
        for i, value in enumerate(data):
            sqtl_dict[header[i]] = value
        sqtl_list.append(sqtl_dict)
snp_list = sqtl_list
#%% Reads a LeafCutter sQTL file and returns a list of dictionaries with the data.
with open(args.junction_file, 'r') as f:
    header = f.readline().strip().split('\t')  # Read and parse the header line
    sqtls = []  # Initialize an empty list to store the data
    for line in f:
        values = line.strip().split('\t')
        row_dict = {header[i]: values[i] for i in range(len(header))}
        sqtls.append(row_dict)
leafcutter_list = sqtls
#%% Function call for Output 1 task:  Retrieves junctions information for a list of SNPs from a LeafCutter file.
# Output: A tuple containing two dictionaries. The first dictionary contains SNP junction information for all SNPs in snp_list.
# The second dictionary contains junction information for SNPs with a p-value < 0.05.
# from IPython import embed; embed() 
  
# Create a dictionary to store the junctions containing each SNP
snp_junctions_tested = {}
return_snp_junctions_tested = {}
junction_set = {}
return_junction_set={}
# Iterate over the SNPs
for snp in snp_list:
    snp_coord = snp['best.snp.coloc']
    snp_coord = '18:8809447:G:C'
    # Dictionary to store the SNP tested junctions
    snp_junctions_tested[snp_coord] = []
    junction_set[snp_coord] = []
    # Iterate over the LeafCutter sQTLs
    for sqtl in leafcutter_list:
        variant_id = sqtl['variant_id'].replace('chr', '')
        if variant_id == snp_coord:
            snp_junctions_tested[snp_coord].append({
                'phenotype_id': sqtl['phenotype_id'],
                'pval_nominal': sqtl['pval_nominal'],
                'slope': sqtl['slope']
            })
            # Extracting SNP-junction set where the SNP is associated with the junction
            if float(sqtl['pval_nominal']) < 0.05:
                junction_set[snp_coord].append({
                'phenotype_id': sqtl['phenotype_id'],
                'pval_nominal': sqtl['pval_nominal'],
                'slope': sqtl['slope']
            })

return_snp_junctions_tested = snp_junctions_tested
return_junction_set = junction_set

## create the intermediate TSV file
with open('A_SNP_junction_set/snp_junctions_tested_only.tsv', 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\t')
    # write the header row
    writer.writerow(['SNP (variant_id)', 'Filtered junctions (phenotype_id)', 'strand', 'pval_nominal', 'slope'])
    # initialize a flag for first row
    jx_set = []
    # iterate through the dictionary and write rows
    # TODO: @Will to add sanity checks to this for loop
    for snp_coord, junctions_info in snp_junctions_tested.items():
        # split the SNP coordinate to get chromosome and position
        chrom, pos, ref, alt = snp_coord.split(':')
        variant_id = f"{chrom}:{pos}:{ref}:{alt}"
        first_row = True
        for i, junction_info in enumerate(junctions_info):
            phenotype_id = junction_info['phenotype_id']
            _, _, _, clu_strand = phenotype_id.split(':')
            strand = str(clu_strand.split('_')[2])
            pval_nominal = junction_info['pval_nominal']
            slope = junction_info['slope']
            
            # write the variant_id only for the first row for a given variant_id
            if first_row:
                writer.writerow([variant_id, phenotype_id, strand, pval_nominal, slope])
                jx_set.append([variant_id, phenotype_id, strand, pval_nominal, slope])
                first_row = False
            else:
                writer.writerow(['', phenotype_id, strand, pval_nominal, slope])
                jx_set.append(['', phenotype_id, strand, pval_nominal, slope])


## create the intermediate TSV file
with open('A_SNP_junction_set/junction_set.tsv', 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\t')
    # write the header row
    writer.writerow(['SNP (variant_id)', 'Filtered junctions (phenotype_id)', 'strand', 'pval_nominal', 'slope'])
    # initialize a flag for first row
    jx_set = []
    # iterate through the dictionary and write rows
    # TODO: @Will to add sanity checks to this for loop
    for snp_coord, junctions_info in junction_set.items():
        # split the SNP coordinate to get chromosome and position
        chrom, pos, ref, alt = snp_coord.split(':')
        variant_id = f"{chrom}:{pos}:{ref}:{alt}"
        first_row = True
        for i, junction_info in enumerate(junctions_info):
            phenotype_id = junction_info['phenotype_id']
            _, _, _, clu_strand = phenotype_id.split(':')
            strand = str(clu_strand.split('_')[2])
            pval_nominal = junction_info['pval_nominal']
            slope = junction_info['slope']
            
            # write the variant_id only for the first row for a given variant_id
            if first_row:
                writer.writerow([variant_id, phenotype_id, strand, pval_nominal, slope])
                jx_set.append([variant_id, phenotype_id, strand, pval_nominal, slope])
                first_row = False
            else:
                writer.writerow(['', phenotype_id, strand, pval_nominal, slope])
                jx_set.append(['', phenotype_id, strand, pval_nominal, slope])


#%% Function call for Output 2 tasks
# Task 1. Identify splice sites from GTF annotation file (rHISAT code modification)
# Output: ss_gtf (pd.DataFrame): DataFrame containing annotated splice sites.

genes = defaultdict(list)
trans = {}
splicesites = []
# Parse valid exon lines from the GTF file into a dict by transcript_id
for line in args.gtf_file:
    line = line.strip()
    if not line or line.startswith('#'):
        continue
    if '#' in line:
        line = line.split('#')[0].strip()

    try:
        chrom, source, feature, left, right, score, \
            strand, frame, values = line.split('\t')
    except ValueError:
        continue
    left, right = int(left), int(right)

    if feature != 'exon' or left >= right:
        continue

    values_dict = {}
    for attr in values.split(';'):
        if attr:
            attr, _, val = attr.strip().partition(' ')
            values_dict[attr] = val.strip('"')

    if 'gene_id' not in values_dict or \
            'transcript_id' not in values_dict:
        continue

    transcript_id = values_dict['transcript_id']
    if transcript_id not in trans:
        trans[transcript_id] = [chrom, strand, [[left, right]]]
        genes[values_dict['gene_id']].append(transcript_id)
    else:
        trans[transcript_id][2].append([left, right])

# Sort exons and merge where separating introns are <=5 bps
for tran, [chrom, strand, exons] in trans.items():
    exons.sort()
    tmp_exons = [exons[0]]
    for i in range(1, len(exons)):
        if exons[i][0] - tmp_exons[-1][1] <= 5:
            tmp_exons[-1][1] = exons[i][1]
        else:
            tmp_exons.append(exons[i])
    trans[tran] = [chrom, strand, tmp_exons]

# Calculate and print the unique junctions and associated transcript IDs
#TODO: Add sanity check for junctions/splicesites etc. etc.
ss = []
junctions = {}
for chrom, strand, exons in trans.values():
    for i in range(1, len(exons)):
        junction = (chrom, exons[i-1][1], exons[i][0], strand)
        transcript_id = [k for k, v in trans.items() if v == [chrom, strand, exons]][0]
        junctions.setdefault(junction, set()).add(transcript_id)

junctions = sorted(junctions.items())      
# Write each junction to the TSV file
for junction, transcript_ids in junctions:
    chrom, left, right, strand = junction
    # Coverting to 0-based coordinates
    left, right = left -1, right -1
    jx_coord = str(left) + '_' + str(right)
    if strand == '-':
        donor2 = str(chrom) + '_' + str(left-2) 
        donor1 = str(chrom) + '_' + str(left-1)
        acceptor1 = str(chrom) + '_' + str(right+1)
        acceptor2 = str(chrom) + '_' + str(right+2)
    else:
        donor1 = str(chrom) + '_' + str(left+1)
        donor2 = str(chrom) + '_' + str(left+2)
        acceptor1 = str(chrom) + '_' + str(right-2) 
        acceptor2 = str(chrom) + '_' + str(right+1)
    ss.append([jx_coord, donor2, strand, 'donor2', ','.join(transcript_ids)])
    ss.append([None, donor1, strand, 'donor1', None])
    ss.append([None, acceptor1, strand, 'acceptor1', None])
    ss.append([None, acceptor2, strand, 'acceptor2', None])

# Create a dataframe from the splice sites list
ss_gtf = pd.DataFrame(ss, columns=['junction_coordinate', 'splicesite_coord', 'strand', 'splicesite_category', 'matched_transcripts'])

# Write the annotated splice sites to a TSV file
ss_gtf.to_csv(f'B_SNP_donor_acceptor_set/splice_sites_gtf.tsv', sep='\t', index=False))

#%% 2. Extracts splice site information from a list of dictionaries containing sQTL data from LeafCutter.
# Output: ss_lc (pd.DataFrame): DataFrame containing annotated splice sites.

# Create an empty list to store splice site information
splice_sites = []

# Iterate over each row in the sqtl list
for row_dict in sqtl_list:
    phenotype_id = row_dict['phenotype_id']
    tss_distance = row_dict['tss_distance']
    chrom, acceptor, donor, clu_strand = phenotype_id.split(':')
    acceptor, donor = int(acceptor), int(donor)

    #TODO: Verify if coloc SNPs coords are 0-based or 1-based
    acceptor, donor = acceptor +1, donor +1
    jx_coord = str(acceptor) + '_' + str(donor)
    cluster = clu_strand.split('_')[1]
    strand = clu_strand.split('_')[2]
    if strand == '-':
        donor2 = 'chr' + str(chrom) + '_' + str(acceptor-2) 
        donor1 = 'chr' + str(chrom) + '_' + str(acceptor-1)
        acceptor1 = 'chr' + str(chrom) + '_' + str(donor+1)
        acceptor2 = 'chr' + str(chrom) + '_' + str(donor+2)
    else:
        donor1 = 'chr' + str(chrom) + '_' + str(acceptor+1)
        donor2 = 'chr' + str(chrom) + '_' + str(acceptor+2)
        acceptor1 = 'chr' + str(chrom) + '_' + str(donor-2) 
        acceptor2 = 'chr' + str(chrom) + '_' + str(donor-1)
    splice_sites.append([phenotype_id, jx_coord, donor2, strand, 'donor2'])
    splice_sites.append([None, None, donor1, strand, 'donor1'])
    splice_sites.append([None, None, acceptor1, strand, 'acceptor1'])
    splice_sites.append([None, None, acceptor2, strand, 'acceptor2'])
    
# Create a dataframe from the splice sites list
ss_lc = pd.DataFrame(splice_sites, columns=['phenotype_id', 'junction_coordinate', 'splicesite_coord', 'strand', 'splicesite_category'])

# Write the annotated splice sites to a TSV file
ss_lc.to_csv(f'B_SNP_donor_acceptor_set/splice_sites_lc.tsv', sep='\t', index=False)

#%% 3. Extracts splice site information from a list of dictionaries containing sQTL data from LeafCutter and GTF 
# that are disrupted by the SNP.
# Output:  snp_disrupts_ss (list): List of dictionaries containing sQTL data from LeafCutter

# Iterate over each row in the sqtl list
snp_disrupts_ss = []
for snp in snp_list:
    snp_coord = snp['best.snp.coloc'].split(":")[1]
    # For mock data simulation
    snp_coord = 8809447

    # Check if the coloc SNPs disrupt the splice sites annotated from LeafCutter sQTL
    for ss in ss_lc.itertuples():
        ss_coord = ss.splicesite_coord.split("_")[1]
        if int(ss_coord) == int(snp_coord):
            snp_disrupts_ss.append([snp['best.snp.coloc'], snp['phenotype'], 'LeafCutter', 'yes', ss.splicesite_coord, ss.splicesite_category, ''])
        else:
            snp_disrupts_ss.append([snp['best.snp.coloc'], snp['phenotype'], 'LeafCutter', '', ss.splicesite_coord, ss.splicesite_category, ''])

    # Check if the coloc SNPs disrupt the splice sites annotated from GTF 
    for ss in ss_gtf.itertuples():
        ss_coord = ss.splicesite_coord.split("_")[1]
        if int(ss_coord) == int(snp_coord):
            snp_disrupts_ss.append([snp['best.snp.coloc'], snp['phenotype'], 'GTF', 'yes', ss.splicesite_coord, ss.splicesite_category, ss.matched_transcripts])
        else:
            snp_disrupts_ss.append([snp['best.snp.coloc'], snp['phenotype'], 'GTF', '', ss.splicesite_coord, ss.splicesite_category, ss.matched_transcripts])

# Create a dataframe from the splice sites list
df = pd.DataFrame(snp_disrupts_ss, columns=['SNP(variant_id)', 'Filtered junctions(phenotype_id)', 'source', 'disrupts splice site', 'splicesite_coord', 'splicesite_category', 'matched_transcripts'])

# Write the annotated splice sites to a TSV file
df.to_csv(f'B_SNP_donor_acceptor_set/snp_disrupts_splice_sites.tsv', sep='\t', index=False))

#%% Function call for Output 3 task

output_file = 'output.tsv'
jx_items = set(item['phenotype_id'] for sublist in jx_set.values() for item in sublist)
snp_items = set(snp_ss[i][1] for i in range(0, len(snp_ss), 4))
overlap_items = jx_items.intersection(snp_items)

# TODO: @Will- Write all attributes to the table + Make the table more readable.
table_rows = []
for key, value in jx_set.items():
    if any(item['phenotype_id'] in overlap_items for item in value):
        snp_matches = []
        for i in range(0, len(snp_ss), 4):
            if (snp_ss[i][1] in overlap_items and key == snp_ss[i][0]):
                snp_matches.append({
                    'splicesite_coord': snp_ss[i+1], 
                    'splicesite_category': snp_ss[i+2],
                    'matched_transcripts': snp_ss[i+3] if snp_ss[i+3] is not None else ''
                })

        for jx in value:
            if jx['phenotype_id'] in overlap_items:
                for snp_match in snp_matches:
                    row = {
                        'SNP(variant_id)': key,
                        'Filtered junctions(phenotype_id)': jx['phenotype_id'],
                
                        'disrupts splice site': '',
                        'splicesite_coord': snp_match['splicesite_coord'],
                        'splicesite_category': snp_match['splicesite_category'],
                        'matched_transcripts': snp_match['matched_transcripts']
                    }
                    table_rows.append(row)

# Write the table to TSV file
with open(output_file, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['SNP(variant_id)', 'Filtered junctions(phenotype_id)', 'disrupts splice site', 'splicesite_coord', 'splicesite_category', 'matched_transcripts'], delimiter='\t')
    writer.writeheader()
    writer.writerows(table_rows)

# To run this script:
# python3 main.py data/MTCL1/MTCL1_QTL_results.tsv data/MTCL1/MTCL1_moloc_results.tsv data/gencode_toy/gencode.v38.toy.gtf




