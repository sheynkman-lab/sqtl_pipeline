<<<<<<< HEAD
#!/usr/bin/env python3
=======
>>>>>>> 86ea07e4b4849f56dd891effd0e6a264801344c3

from helpers import read_leafcutter_sQTL_file, read_moloc_sqtl_file
import logging
import csv

def get_snp_junctions(snp_file_path, leafcutter_file_path, ):
    """
    Given a file with colocalized SNPs and a file with LeafCutter sQTL data, returns a dictionary where the key is the SNP
    coordinate and the value is a list of dictionaries containing junction info from LeafCutter sQTL file.

    Args:
        snp_file_path (str): Path to the file with colocalized SNPs (from args).
        leafcutter_file_path (str): Path to the file with LeafCutter sQTL data (from).

    Returns:
        Filtered SNP junctions dictionary.

    Outputs:
        An intermediate CSV file containing the filtered junctions (snp_junctions.csv).
    """
    
    # Read the colocalized SNP file
    snp_list = read_moloc_sqtl_file(snp_file_path)
    logging.info(f'Read {len(snp_list)} SNPs from {snp_file_path}')

    # Read the LeafCutter sQTL file
    leafcutter_list = read_leafcutter_sQTL_file(leafcutter_file_path)
    logging.info(f'Read {len(leafcutter_list)} junctions from {leafcutter_file_path}')

    # Create a dictionary to store the junctions containing each SNP
    snp_junctions = {}

    # Iterate over the SNPs
    for snp in snp_list:
        snp_coord = snp['best.snp.coloc']
        snp_junctions[snp_coord] = []

        # Iterate over the LeafCutter sQTLs
        for sqtl in leafcutter_list:
            variant_id = sqtl['variant_id'].replace('chr', '')
            if variant_id == snp_coord and float(sqtl['pval_nominal']) < 0.05:
                snp_junctions[snp_coord].append({
                    'phenotype_id': sqtl['phenotype_id'],
                    'pval_nominal': sqtl['pval_nominal'],
                    'slope': sqtl['slope']
                })

    
    ## create the intermediate CSV file
    with open('A_SNP_junction_set/snp_junctions.csv', 'w', newline='') as f:
        writer = csv.writer(f)

        # write the header row
        writer.writerow(['SNP (variant_id)', 'Filtered junctions (phenotype_id)', 'strand', 'pval_nominal', 'slope'])
        
        # initialize a flag for first row
        first_row = True

        # iterate through the dictionary and write rows
        for snp_coord, junctions_info in snp_junctions.items():
            # split the SNP coordinate to get chromosome and position
            chrom, pos, ref, alt = snp_coord.split(':')
            variant_id = f"{chrom}:{pos}:{ref}:{alt}"
            
            for i, junction_info in enumerate(junctions_info):
                phenotype_id = junction_info['phenotype_id']
                _, _, _, clu_strand = phenotype_id.split(':')
                strand = str(clu_strand.split('_')[2])
                pval_nominal = junction_info['pval_nominal']
                slope = junction_info['slope']
                
                # write the variant_id only for the first row for a given variant_id
                if first_row:
                    writer.writerow([variant_id, phenotype_id, strand, pval_nominal, slope])
                    first_row = False
                else:
                    writer.writerow(['', phenotype_id, strand, pval_nominal, slope])


    return snp_junctions

					