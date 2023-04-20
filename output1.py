import pandas as pd
import networkx as nx
from splice_graph import SpliceGraph
import csv

def get_junction_set(snp, splicegraph, p_value_threshold=0.05):
    # Create a subgraph containing only the junctions containing the SNP
    subgraph = SpliceGraph()
    for node in splicegraph.vertices:
        if snp in node:
            subgraph.vertices[node] = splicegraph.vertices[node]
            subgraph.edges[node] = splicegraph.edges[node]

    # Filter junctions based on p-value threshold
    junction_set = []
    for node in subgraph.vertices:
        if 'sqtl_data' in subgraph.edges[node]:
            for phenotype_id in subgraph.edges[node]['sqtl_data']:
                for sqtl in subgraph.edges[node]['sqtl_data'][phenotype_id]:
                    if sqtl['pval_nominal'] <= p_value_threshold:
                        junction_set.append(sqtl['junction_id'])
    
    # Write to a csv file
    with open(f"{snp}_junction_set.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Junction ID"])
        for junction_id in junction_set:
            writer.writerow([junction_id])
    
    return junction_set

