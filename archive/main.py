import argparse
import logging
import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
from splice_graph import SpliceGraph
from helpers import read_leafcutter_sQTL_file, read_gwas_snp_file, read_moloc_sqtl_file
from analysis import get_filtered_phenotypes, get_filtered_sqtl_junctions, print_sqtl_data
from output1 import get_junction_set

def main():
    parser = argparse.ArgumentParser(description='Analyze splice junctions and GWAS SNPs')
    parser.add_argument('junction_file', help='Path to the Leafcutter sQTL junction file')
    parser.add_argument('gwas_snp_file', help='Path to the GWAS SNP file')
    parser.add_argument('moloc_snp_file', help='Path to the moloc SNP file')
    parser.add_argument('--plot', action='store_true', help='Generate a plot of the splice graph')

    args = parser.parse_args()

    # set up logging configuration
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s', filename='output.log', filemode='w')

    #read LeafCutter sQTL file
    leafcutter_sqtl = read_leafcutter_sQTL_file(args.junction_file)
    logging.info(f'Read {len(leafcutter_sqtl)} junctions from {args.junction_file}')

    # create a SpliceGraph object and add the junctions
    sg = SpliceGraph()
    for junction in leafcutter_sqtl:
        sg.add_junction(junction)
    sg.generate_splice_graph()

    #read GWAS SNP file
    gwas_snp_data = read_gwas_snp_file(args.gwas_snp_file)
    logging.info(f'Read {len(gwas_snp_data)} SNPs from {args.gwas_snp_file}')

    #read moloc file
    moloc_snp_data = read_moloc_sqtl_file(args.moloc_snp_file)
    logging.info(f'Read {len(moloc_snp_data)} SNPs from {args.moloc_snp_file}')

    # create a NetworkX graph object from the SpliceGraph
    nx_graph = nx.Graph(sg.edges)

    # Print all the vertices
    for vertex in sg.vertices:
        logging.info(f'Vertex of splicegraph {vertex}')

    # Print all the edges
    for vertex, edges in sg.edges.items():
        for edge in edges:
            logging.info(f'Edge of splicegraph {vertex} -> {edge} ')

    # Visualize the graph
    # nx.draw(nx_graph, with_labels=True)
    # plt.show()

    
    phenotype_cluster = 'clu_2300'
    marker = 'chr18:8809447:G:C'
    max_pval_nominal = 0.05


    # Updated code for Output 1(https://docs.google.com/document/d/1OFUz3p-P1ZITwIBa6-zXuqV--Hc7dhjx7anv4kB1PBE/edit#heading=h.hi4ej0bwwneq)

    jx_set = get_junction_set(marker, sg)


    # filtered_phenotypes = get_filtered_phenotypes(sg, phenotype_cluster)
    # logging.info(f'Phenotype cluster provided for sQTL filteration : {phenotype_cluster}')
    # for phenotype in filtered_phenotypes:
    #     logging.info(f'Filtered phenotype : {phenotype}')

    # logging.info(f'Marker and p-val provided for sQTL filteration : {marker}, {max_pval_nominal}')
    # filtered_sqtl_junctions = get_filtered_sqtl_junctions (sg, filtered_phenotypes, marker, max_pval_nominal)
    # for sqtl_row in filtered_sqtl_junctions:
    #     logging.info(f'Filtered LeafCutter sQTL data : {sqtl_row}')

    # print_sqtl_data(sg, '18:8720496:8762078:clu_2297_+')

    # df = pd.DataFrame(filtered_junctions)
    # df.to_csv(f'trial.csv', sep='\t', index=False)

if __name__ == '__main__':
    main()

# To run this script:
# python3 main.py data/MTCL1/MTCL1_QTL_results.tsv data/MTCL1/MTCL1_GWAS_results.tsv data/MTCL1/MTCL1_moloc_results.tsv


# The vertices of this splice graph represent the different junctions and SNP positions, 
# while the edges represent the connections between them. Specifically, each junction is represented 
# by two vertices - one for the donor site and one for the acceptor site - and edges connect the donor 
# and acceptor vertices of the same junction as well as the SNP vertices to the corresponding donor and 
# acceptor vertices. This results in a complex network of edges connecting different junctions and SNPs, 
# representing the potential paths of splicing events in the gene being analyzed.


