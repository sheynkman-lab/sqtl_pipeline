#!/usr/bin/env python3

from collections import Counter, defaultdict
from sys import stderr, exit
import csv


def extract_splice_sites_gtf(gtf_file):
    """
    Extracts splice site information GTF annotation file.
    From rHISAT

    Args:
        gtf_file (str): Path to the GTF annotation file.

    Returns:
        None. Writes the annotated splice sites to a CSV file.
    """

    genes = defaultdict(list)
    trans = {}
    splicesites = []
    # Parse valid exon lines from the GTF file into a dict by transcript_id
    for line in gtf_file:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        if '#' in line:
            line = line.split('#')[0].strip()

        try:
            chrom, source, feature, left, right, score, \
                strand, frame, values = line.split('\t')
        except ValueError:
            continue
        left, right = int(left), int(right)

        if feature != 'exon' or left >= right:
            continue

        values_dict = {}
        for attr in values.split(';'):
            if attr:
                attr, _, val = attr.strip().partition(' ')
                values_dict[attr] = val.strip('"')

        if 'gene_id' not in values_dict or \
                'transcript_id' not in values_dict:
            continue

        transcript_id = values_dict['transcript_id']
        if transcript_id not in trans:
            trans[transcript_id] = [chrom, strand, [[left, right]]]
            genes[values_dict['gene_id']].append(transcript_id)
        else:
            trans[transcript_id][2].append([left, right])

    # Sort exons and merge where separating introns are <=5 bps
    for tran, [chrom, strand, exons] in trans.items():
        exons.sort()
        tmp_exons = [exons[0]]
        for i in range(1, len(exons)):
            if exons[i][0] - tmp_exons[-1][1] <= 5:
                tmp_exons[-1][1] = exons[i][1]
            else:
                tmp_exons.append(exons[i])
        trans[tran] = [chrom, strand, tmp_exons]

    # Calculate and print the unique junctions and associated transcript IDs
    #TODO: Add sanity check for junctions/splicesites etc. etc.
    junctions = {}
    for chrom, strand, exons in trans.values():
        for i in range(1, len(exons)):
            junction = (chrom, exons[i-1][1], exons[i][0], strand)
            transcript_id = [k for k, v in trans.items() if v == [chrom, strand, exons]][0]
            junctions.setdefault(junction, set()).add(transcript_id)

    junctions = sorted(junctions.items())
    with open(f'B_SNP_donor_acceptor_set/splice_sites_gtf.csv', mode='w', newline='') as csv_file:
        writer = csv.writer(csv_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        
        # Write the header row to the CSV file
        writer.writerow(['splicesite_coord', 'strand', 'splicesite_category', 'matched_transcripts'])
        
        # Write each junction to the CSV file
        for junction, transcript_ids in junctions:
            chrom, left, right, strand = junction
            left, right = left -1, right -1
            if strand == '-':
                donor2 = str(chrom) + '_' + str(left-2) 
                donor1 = str(chrom) + '_' + str(left-1)
                acceptor1 = str(chrom) + '_' + str(right+1)
                acceptor2 = str(chrom) + '_' + str(right+2)
            else:
                donor1 = str(chrom) + '_' + str(left+1)
                donor2 = str(chrom) + '_' + str(left+2)
                acceptor1 = str(chrom) + '_' + str(right-2) 
                acceptor2 = str(chrom) + '_' + str(right+1)

            writer.writerow([donor2, strand, 'donor2', ','.join(transcript_ids)])
            writer.writerow([donor1, strand, 'donor1', None])
            writer.writerow([acceptor1, strand, 'acceptor1', None])
            writer.writerow([acceptor2, strand, 'acceptor2', None])
        
    # exon_lengths, intron_lengths, trans_lengths = \
    #     Counter(), Counter(), Counter()
    # for chrom, strand, exons in trans.values():
    #     tran_len = 0
    #     for i, exon in enumerate(exons):
    #         exon_len = exon[1]-exon[0]+1
    #         exon_lengths[exon_len] += 1
    #         tran_len += exon_len
    #         if i == 0:
    #             continue
    #         intron_lengths[exon[0] - exons[i-1][1]] += 1
    #     trans_lengths[tran_len] += 1


    # print('genes: {}, genes with multiple isoforms: {}'.format(
    #         len(genes), sum(len(v) > 1 for v in genes.values())),
    #         file=stderr)
    # print('transcripts: {}, transcript avg. length: {:.0f}'.format(
    #         len(trans), sum(trans_lengths.elements())//len(trans)),
    #         file=stderr)
    # print('exons: {}, exon avg. length: {:.0f}'.format(
    #         sum(exon_lengths.values()),
    #         sum(exon_lengths.elements())//sum(exon_lengths.values())),
    #         file=stderr)
    # print('introns: {}, intron avg. length: {:.0f}'.format(
    #         sum(intron_lengths.values()),
    #         sum(intron_lengths.elements())//sum(intron_lengths.values())),
    #         file=stderr)
    # print('average number of exons per transcript: {:.0f}'.format(
    #         sum(exon_lengths.values())//len(trans)),
    #         file=stderr)


def extract_splice_sites_lc(sqtl_list):
    """
    Extracts splice site information from a list of dictionaries containing sQTL data from LeafCutter.

    Args:
        sqtl_list (list): List of dictionaries containing sQTL data.

    Returns:
        None. Writes the annotated splice sites to a CSV file.
    """
    with open(f'B_SNP_donor_acceptor_set/splice_sites_lc.csv', mode='w', newline='') as csv_file:
        writer = csv.writer(csv_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        
        # Write the header row to the CSV file
        writer.writerow(['splicesite_coord', 'strand', 'splicesite_category'])
        
        # Iterate over each row in the sqtl list
        for row_dict in sqtl_list:
            phenotype_id = row_dict['phenotype_id']
            tss_distance = row_dict['tss_distance']
            chrom, acceptor, donor, clu_strand = phenotype_id.split(':')
            acceptor, donor = int(acceptor), int(donor)
            cluster = clu_strand.split('_')[1]
            strand = clu_strand.split('_')[2]
            if strand == '-':
                donor2 = 'chr' + str(chrom) + '_' + str(acceptor-2) 
                donor1 = 'chr' + str(chrom) + '_' + str(acceptor-1)
                acceptor1 = 'chr' + str(chrom) + '_' + str(donor+1)
                acceptor2 = 'chr' + str(chrom) + '_' + str(donor+2)
            else:
                donor1 = 'chr' + str(chrom) + '_' + str(acceptor+1)
                donor2 = 'chr' + str(chrom) + '_' + str(acceptor+2)
                acceptor1 = 'chr' + str(chrom) + '_' + str(donor-2) 
                acceptor2 = 'chr' + str(chrom) + '_' + str(donor+1)
            writer.writerow([donor2, strand, 'donor2'])
            writer.writerow([donor1, strand, 'donor1'])
            writer.writerow([acceptor1, strand, 'acceptor1'])
            writer.writerow([acceptor2, strand, 'acceptor2'])

# def snp_disrupts_splice_sites():
#     sdf