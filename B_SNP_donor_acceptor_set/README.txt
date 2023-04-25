Output 2: 

Output 2 deteremines whether the SNPs disrupt slice donor or acceptor sites. #Will expand more

Input
In order to produce the second output the following is needed:
1) GTF file with transcript annotations. Each row of the file contains the gene name, start and end position of the gene, feature, reading frame, and strand (+ for forward and - for reverse).

Steps
1) The first step begins with reading in the GTF file
2) The GTF is then filtered for rows that are defined as an exon under the feature column. Additionally the acceptor and donor sites for each gene is located.
2) Once read, the SNP junctions from Output 1 are compared to the splice sites from the GTF. This is used to determine whether the SNPs of interest fall within the splice site.
3) Additionally, SNPs are also compared to intron regions. For SNPs that fall within introns, the location can be determined.

Output
This process produces a table that contains the junction coordinates, splice coordinates, strand (forward or reverse), acceptor/ donor classification, and associated transcripts with each splice site.
