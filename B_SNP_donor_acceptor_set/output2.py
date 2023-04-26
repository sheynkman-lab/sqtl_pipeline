#!/usr/bin/env python3

from collections import Counter, defaultdict
from sys import stderr, exit
import csv
import pandas as pd


def extract_splice_sites_gtf(gtf_file):
    """
    Extracts splice site information GTF annotation file.
    From rHISAT

    Args:
        gtf_file (str): Path to the GTF annotation file.

    Returns:
        ss_gtf (pd.DataFrame): DataFrame containing annotated splice sites.
    """

    genes = defaultdict(list)
    trans = {}
    splicesites = []
    # Parse valid exon lines from the GTF file into a dict by transcript_id
    for line in gtf_file:
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
    # Write each junction to the CSV file
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
    
    # Write the annotated splice sites to a CSV file
    ss_gtf.to_csv(f'B_SNP_donor_acceptor_set/splice_sites_gtf.csv', index=False)
    
    # Return the dataframe
    return ss_gtf



def extract_splice_sites_lc(sqtl_list):
    """
    Extracts splice site information from a list of dictionaries containing sQTL data from LeafCutter.

    Args:
        sqtl_list (list): List of dictionaries containing sQTL data.

    Returns:
        ss_lc (pd.DataFrame): DataFrame containing annotated splice sites.
    """
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
    
    # Write the annotated splice sites to a CSV file
    ss_lc.to_csv(f'B_SNP_donor_acceptor_set/splice_sites_lc.csv', index=False)
    
    # Return the dataframe
    return ss_lc

def snp_disrupts_splice_sites(snp_list, leafcutter_list, ss_lc, ss_gtf):
    """
    Extracts splice site information from a list of dictionaries containing sQTL data from LeafCutter and GTF 
    that are disrupted by the SNP.

    Args:
        snp_list (list): List of dictionaries containing sQTL data.
        leafcutter_list (list): List of dictionaries containing sQTL data from LeafCutter.
        ss_lc (pd.DataFrame): DataFrame containing annotated splice sites.
        ss_gtf (pd.DataFrame): DataFrame containing annotated splice sites.

    Returns:
        snp_disrupts_ss (list): List of dictionaries containing sQTL data from LeafCutter.
    """
    
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

    # Write the annotated splice sites to a CSV file
    df.to_csv(f'B_SNP_donor_acceptor_set/snp_disrupts_splice_sites.csv', index=False)
    
    # Return the dataframe
    return snp_disrupts_ss 

