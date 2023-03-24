# from my_tool.helpers import read_gtf, get_exon_boundaries, is_exon, get_gene_name
# from my_tool.splice_graph import SpliceGraph, read_sqtl_junctions, check_snp_in_junctions
# from my_tool.gwas_sno import GWASSNO, read_gwas_snps

# __all__ = [
#     'read_gtf',
#     'get_exon_boundaries',
#     'is_exon',
#     'get_gene_name',
#     'SpliceGraph',
#     'read_sqtl_junctions',
#     'check_snp_in_junctions',
#     'GWASSNO',
#     'read_gwas_snps'
# ]


# # use the SpliceGraph class to create a splice graph
# sg = SpliceGraph()
# sg.add_junctions(read_sqtl_junctions('junctions.txt'))

# # use the read_gwas_snps function to read in GWAS SNPs
# snps = read_gwas_snps('gwas_snps.txt')
# for snp in snps:
#     sg.annotate_snp(snp['chrom'], snp['pos'], snp['rsid'])
