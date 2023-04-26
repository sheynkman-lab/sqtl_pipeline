import csv

def find_set_overlaps(jx_set, snp_ss):
    output_file = 'output.csv'
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
    
    # Write the table to CSV file
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['SNP(variant_id)', 'Filtered junctions(phenotype_id)', 'disrupts splice site', 'splicesite_coord', 'splicesite_category', 'matched_transcripts'])
        writer.writeheader()
        writer.writerows(table_rows)
