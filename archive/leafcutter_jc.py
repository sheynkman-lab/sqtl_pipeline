class LeafCutterJunctionCluster:
    """
    A class to hold information about a single LeafCutter junction cluster.
    """
    
    def __init__(self, chrom, start, end, strand, name, junctions):
        """
        Initializes the LeafCutterJunctionCluster object.
        
        Args:
        chrom (str): The chromosome where the cluster is located.
        start (int): The start position of the cluster.
        end (int): The end position of the cluster.
        strand (str): The strand of the cluster ('+' or '-').
        name (str): The name of the cluster.
        junctions (list): A list of junctions that belong to the cluster.
        """
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.name = name
        self.junctions = junctions
    
    def __str__(self):
        """
        Returns a string representation of the object.
        """
        return f"{self.name} ({self.chrom}:{self.start}-{self.end})"
    

class LeafCutterJunctionClusterReader:
    """
    A class to read LeafCutter junction cluster files and return
    a list of LeafCutterJunctionCluster objects.
    """
    
    def __init__(self, cluster_file):
        """
        Initializes the LeafCutterJunctionClusterReader object.
        
        Args:
        cluster_file (str): The path to the cluster file.
        """
        self.cluster_file = cluster_file
    
    def read_clusters(self):
        """
        Reads the cluster file and returns a list of LeafCutterJunctionCluster objects.
        """
        clusters = []
        with open(self.cluster_file, 'r') as f:
            for line in f:
                if line.startswith('track'):
                    continue
                fields = line.strip().split('\t')
                chrom, start, end, strand, name = fields[:5]
                junctions = [int(x) for x in fields[5:]]
                cluster = LeafCutterJunctionCluster(chrom, int(start), int(end), strand, name, junctions)
                clusters.append(cluster)
        return clusters
