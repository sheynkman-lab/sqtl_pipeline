#!/usr/bin/env python
import csv

def get_snp_junctions(snp_list, leafcutter_list):
    """
    Retrieves junctions information for a list of SNPs from a LeafCutter file.

    Args:
        snp_list (list): A list of dictionaries containing SNP information with a 'best.snp.coloc' key.
        leafcutter_list (list): A list of dictionaries containing junction information with a 'variant_id' key.

    Returns:
        tuple: A tuple containing two dictionaries. The first dictionary contains SNP junction information for all SNPs in snp_list.
        The second dictionary contains junction information for SNPs with a p-value < 0.05.

    Example:
        snp_list = [{'best.snp.coloc': '18:8809447:G:C'}, {'best.snp.coloc': '18:8809448:G:T'}]
        leafcutter_list = [{'variant_id': '18:8809447:G:C:1'}, {'variant_id': '18:8809447:G:C:2'}, {'variant_id': '18:8809448:G:T:1'}]
        snp_junctions, junction_set = get_snp_junctions(snp_list, leafcutter_list)
    """
    
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

    ## create the intermediate CSV file
    with open('A_SNP_junction_set/snp_junctions_tested_only.csv', 'w', newline='') as f:
        writer = csv.writer(f)
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


    ## create the intermediate CSV file
    with open('A_SNP_junction_set/junction_set.csv', 'w', newline='') as f:
        writer = csv.writer(f)
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

    return return_snp_junctions_tested, return_junction_set

					