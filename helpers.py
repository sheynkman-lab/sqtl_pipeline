# Ways to import the below functions
# from my_tool.helpers import read_gtf, get_exon_boundaries, is_exon, get_gene_name
# from my_tool.splice_graph import SpliceGraph, read_sqtl_junctions, check_snp_in_junctions
# from my_tool.gwas_sno import GWASSNO, read_gwas_snps


def read_gtf(filename):
    ...


def get_exon_boundaries(gtf):
    ...


def is_exon(chrom, pos, exon_boundaries):
    ...


def get_gene_name(chrom, pos, gtf):
    ...


def check_snp_in_junctions(sg, snp_chrom, snp_pos):
    ...


def read_gwas_snp_file(file_path):
    """
    Reads a GWAS SNP file and returns a list of dictionaries with the data.

    Args:
        file_path (str): Path to the GWAS SNP file.

    Returns:
        A list of dictionaries with the data from the file.
    """
    with open(file_path, 'r') as f:
        header = f.readline().strip().split('\t')  # Read and parse the header line
        rows = []  # Initialize an empty list to store the data
        for line in f:
            values = line.strip().split('\t')
            row_dict = {header[i]: values[i] for i in range(len(header))}
            rows.append(row_dict)
    return rows


def read_leafcutter_sQTL_file(file_path):
    """
    Reads a LeafCutter sQTL file and returns a list of dictionaries with the data.

    Args:
        file_path (str): Path to the LeafCutter sQTL file.

    Returns:
        A list of dictionaries with the data from the file.
    """
    with open(file_path, 'r') as f:
        header = f.readline().strip().split('\t')  # Read and parse the header line
        sqtls = []  # Initialize an empty list to store the data
        for line in f:
            values = line.strip().split('\t')
            row_dict = {header[i]: values[i] for i in range(len(header))}
            sqtls.append(row_dict)
    return sqtls


def read_moloc_sqtl_file(file_path):
    """
    Reads a MOLoc sQTL file and returns a list of dictionaries with the data.

    Args:
        file_path (str): Path to the MOLoc sQTL file.

    Returns:
        A list of dictionaries with the data from the file.
    """
    sqtl_list = []
    with open(file_path, 'r') as f:
        header = f.readline().strip().split('\t')
        for line in f:
            data = line.strip().split('\t')
            sqtl_dict = {}
            for i, value in enumerate(data):
                sqtl_dict[header[i]] = value
            sqtl_list.append(sqtl_dict)
    return sqtl_list
