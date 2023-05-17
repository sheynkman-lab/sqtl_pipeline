# Donor-Acceptor Set

## Input
To deteremines whether the SNPs disrupt slice donor or acceptor sites. 
1) GTF file with transcript annotations. Each row of the file contains the gene name, start and end position of the gene, feature, reading frame, and strand (+ for forward and - for reverse).

## Steps
1) Identify the splice sites from GTF annotation file (to bring intron-exon transcript information for splice disruption sites)
2) Identify the splice sites from LeafCutter input file 
3) Check if SNP falls within the splice site regions (GTF + LC)

## Output
This process produces a table that contains the junction coordinates, splice coordinates, strand (forward or reverse), acceptor/ donor classification, and associated transcripts with each splice site.
