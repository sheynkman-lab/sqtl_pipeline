class GWASSNO:
    def __init__(self, chrom, pos, rsid, pval):
        self.chrom = chrom
        self.pos = pos
        self.rsid = rsid
        self.pval = pval


def read_gwas_snps(filename):
    snps = []
    with open(filename, 'r') as f:
        for line in f:
            fields = line.strip().split()
            chrom = fields[0]
            pos = int(fields[1])
            rsid = fields[2]
            pval = float(fields[3])
            snps.append(GWASSNO(chrom, pos, rsid, pval))
    return snps
