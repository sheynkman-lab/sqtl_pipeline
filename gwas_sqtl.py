class GWASData:
    """
    A class to hold GWAS summary statistics and corresponding SQTL associations.
    """
    
    def __init__(self, gwas_file):
        """
        Initializes the GWASData object.
        
        Args:
        gwas_file (str): The path to the GWAS summary statistics file.
        """
        self.gwas_file = gwas_file
        self.snps = []
        self.sqtls = []
        
    def add_sqtl(self, chromosome, position, ref_allele, alt_allele, gene_id, tss_distance, beta, se, pval):
        """
        Adds a new SQTL association to the GWASData object. Assuming the user inputs the colocalized sQTLS

        Args:
        chromosome (str): The chromosome where the SNP is located.
        position (int): The genomic position of the SNP.
        ref_allele (str): The reference allele of the SNP.
        alt_allele (str): The alternative allele of the SNP.
        gene_id (str): The gene ID of the SQTL.
        tss_distance (int): The distance between the SNP and the transcription start site of the gene.
        beta (float): The effect size (beta) of the SNP on the expression level of the gene.
        se (float): The standard error of the effect size (beta).
        pval (float): The p-value of the association.
        """
        sqtl = GWASSQTL(chromosome, position, ref_allele, alt_allele, gene_id, tss_distance, beta, se, pval)
        self.sqtls.append(sqtl)
        
    def read_gwas_file(self):
        """
        Reads the GWAS summary statistics file and stores the information in the GWASData object.
        """
        with open(self.gwas_file) as f:
            header = f.readline().strip().split()
            for line in f:
                fields = line.strip().split()
                chromosome = fields[0]
                position = int(fields[1])
                ref_allele = fields[2]
                alt_allele = fields[3]
                beta = float(fields[4])
                se = float(fields[5])
                pval = float(fields[6])
                snp = GWASSNP(chromosome, position, ref_allele, alt_allele, beta, se, pval)
                self.snps.append(snp)
