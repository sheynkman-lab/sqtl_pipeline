class SQTLJunctionMapper:
    """
    A class to map colocalized SQTLs to LeafCutter junctions.
    """
    
    def __init__(self, sqtls, junction_clusters):
        """
        Initializes the SQTLJunctionMapper object.
        
        Args:
        sqtls (list): A list of GWASSQTL objects.
        junction_clusters (list): A list of LeafCutterJunctionCluster objects.
        """
        self.sqtls = sqtls
        self.junction_clusters = junction_clusters
        self.mappings = {}
        
    def map_sqtls_to_junctions(self):
        """
        Maps the colocalized SQTLs to LeafCutter junction clusters.
        """
        for sqtl in self.sqtls:
            for cluster in self.junction_clusters:
                if sqtl.chromosome == cluster.chromosome and cluster.contains(sqtl.position):
                    if cluster.cluster_id in self.mappings:
                        self.mappings[cluster.cluster_id].append(sqtl)
                    else:
                        self.mappings[cluster.cluster_id] = [sqtl]
                        
    def get_junction_sqtls(self, cluster_id):
        """
        Returns the colocalized SQTLs for a given junction cluster ID.
        
        Args:
        cluster_id (str): The ID of the junction cluster.
        
        Returns:
        A list of GWASSQTL objects that are colocalized with the junction cluster.
        """
        if cluster_id in self.mappings:
            return self.mappings[cluster_id]
        else:
            return []


    def map_sqtls_to_junctions_2(self):
        sqtl_junction_map = {}

        # Loop through the SQTLs and junctions and map each SQTL to its corresponding junction(s)
        for sqtl in sqtls:
            sqtl_chr, sqtl_pos = sqtl
            sqtl_key = (sqtl_chr, sqtl_pos)
            sqtl_junction_map[sqtl_key] = []
            for junction in junctions:
                junction_chr, junction_pos = junction
                if sqtl_chr == junction_chr and abs(sqtl_pos - junction_pos) <= 500:
                    sqtl_junction_map[sqtl_key].append(junction)






class MaserOutput:
    def __init__(self, filename):
        self.filename = filename
        self.data = []
        with open(filename, 'r') as f:
            for line in f:
                if line.startswith('GeneID'):
                    continue
                fields = line.strip().split('\t')
                gene_id = fields[0]
                gene_name = fields[1]
                gene_type = fields[2]
                expression = float(fields[3])
                tp = int(fields[4])
                tn = int(fields[5])
                fp = int(fields[6])
                fn = int(fields[7])
                auc = float(fields[8])
                f1score = float(fields[9])
                pvalue = float(fields[10])
                self.data.append({
                    'gene_id': gene_id,
                    'gene_name': gene_name,
                    'gene_type': gene_type,
                    'expression': expression,
                    'tp': tp,
                    'tn': tn,
                    'fp': fp,
                    'fn': fn,
                    'auc': auc,
                    'f1score': f1score,
                    'pvalue': pvalue
                })
    
    def get_significant_genes(self, pvalue_cutoff=0.05):
        """
        Get significant genes based on the P-value cutoff.
        """
        significant_genes = []
        for gene in self.data:
            if gene['pvalue'] < pvalue_cutoff:
                significant_genes.append(gene)
        return significant_genes
