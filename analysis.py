import pandas as pd
import networkx as nx

   
def parse_junction_id(junction_id):
    """Extract chromosome, start, end, and cluster, strand information from a junction ID."""
    # print(junction_id)
    fields = junction_id.split(":")
    chrom = fields[0]
    start = int(fields[1])
    end = int(fields[2])
    junction_coords = fields[3].split("_")
    cluster = junction_coords[0] + '_' + junction_coords[1]
    strand = junction_coords[2]
    # intron_start = junction_coords[3]
    # intron_end = junction_coords[4]
    # type_d_a = junction_coords[5]
    return chrom, start, end, cluster, strand


def get_variant_ids(junction_ids, splice_graph):
    """Return a list of variant IDs corresponding to the given junction IDs."""
    variant_ids = []
    for junction_id in junction_ids:
        chrom, start, end, cluster, strand = parse_junction_id(junction_id)
        
        jx_id = str(chrom) + ':' + str(start) + ':' + str(end) + ':' + str(cluster) + '_' + strand
        variant_ids.append(splice_graph.edges[jx_id]['sqtl_data']['variant_id'])
    return variant_ids

def get_filtered_phenotypes(sg, phenotype_cluster):
    """
    This function takes in the moloc_snp_data, gwas_snp_data, and splice_graph dataframes as arguments, and returns
    a list of all junctions in the splice graph that have the phenotype cluster and GWAS SNP specified in the prompt.

    Note that the function assumes that the phenotype_id and snp_id fields in the moloc_snp_data and gwas_snp_data
    dataframes respectively are unique. If this is not the case, the function may not return the desired result.
    """
    sqtl_junctions = []
    for junction_id in sg.edges:
        chrom, start, end, cluster, strand = parse_junction_id(junction_id)
        parsed_id = str(chrom) + ':' + str(start) + ':' + str(end) + ':' + str(cluster) + '_' + strand
        if 'sqtl_data' in sg.edges[junction_id] and cluster == phenotype_cluster:
            sqtl_junctions.append(junction_id)
    return sqtl_junctions


def get_filtered_sqtl_junctions(sg, filtered_junctions, variant_id_str, max_pval_nominal):
    filtered_sqtl_junctions = []
    for junction_id in filtered_junctions:
        if 'sqtl_data' in sg.edges[junction_id]:
            sqtl_data = sg.edges[junction_id]['sqtl_data']
            for phenotype_id in sqtl_data:
                for data_dict in sqtl_data[phenotype_id]:
                    # print(f"variant_id: {data_dict['variant_id']}")
                    if data_dict['variant_id'] == variant_id_str and float(data_dict['pval_nominal']) <= float(max_pval_nominal):
                        filtered_sqtl_junctions.append(data_dict)
    return filtered_sqtl_junctions


def print_sqtl_data(sg, junction_id):
    if 'sqtl_data' in sg.edges[junction_id]:
        sqtl_data = sg.edges[junction_id]['sqtl_data']
        for phenotype_id in sqtl_data:
            print(f"SQTL data for phenotype ID {phenotype_id}:")
            for data_dict in sqtl_data[phenotype_id]:
                for key, value in data_dict.items():
                    print(f"{key}: {value}")
                print("")
    else:
        print("No SQTL data for this junction ID.")





