class RMATSOutput:
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
                chrom = fields[2]
                strand = fields[3]
                exon_id = fields[4]
                upstream_exon = fields[5]
                downstream_exon = fields[6]
                pvalue = float(fields[7])
                fdr = float(fields[8])
                inc_diff = float(fields[9])
                inc1 = float(fields[10])
                inc2 = float(fields[11])
                inc1_sd = float(fields[12])
                inc2_sd = float(fields[13])
                inc_diff_sd = float(fields[14])
                self.data.append({
                    'gene_id': gene_id,
                    'gene_name': gene_name,
                    'chrom': chrom,
                    'strand': strand,
                    'exon_id': exon_id,
                    'upstream_exon': upstream_exon,
                    'downstream_exon': downstream_exon,
                    'pvalue': pvalue,
                    'fdr': fdr,
                    'inc_diff': inc_diff,
                    'inc1': inc1,
                    'inc2': inc2,
                    'inc1_sd': inc1_sd,
                    'inc2_sd': inc2_sd,
                    'inc_diff_sd': inc_diff_sd
                })

    def get_significant_splicing_events(self, fdr_cutoff=0.05, inc_diff_cutoff=0.1):
        """
        Get significant splicing events based on the FDR and inclusion level difference cutoffs.
        """
        significant_events = []
        for event in self.data:
            if event['fdr'] < fdr_cutoff and abs(event['inc_diff']) > inc_diff_cutoff:
                significant_events.append(event)
        return significant_events
