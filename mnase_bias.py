import coord_ops as co
import csv
import housekeeping as hk
import os
import random
import read_and_write as rw

def map_kmers_to_positions(fasta, k = 2, focal_pos = 1):
    """
    Given a FASTA with headers formatted as "transcript_id(strand)",
    make a dictionary where for each k-mer, you indicate all the positions
    in the FASTA that end in that k-mer. NB! The FASTA should include however much of the pre-transcript sequence
    should be included in the first k-mer.
    :param fasta: FASTA file (sequences cannot break across lines)
    :param k: order of the k-mers
    :param focal_pos: the position of the read start
    :return: dictionary with dinculeotide strings as keys and list of tuples as values,
    each tuple consisting of (transcript_id, 0-based position counted from position -1)
    """
    kmer_dict = {}
    counter = 0
    with open(fasta) as f_file:
        for line in f_file:
            if line[0] == ">":
                # removing the greater than sign, the strand and the newline
                prev = line[1:-4]
            else:
                counter = hk.update_counter(counter, 1000, "{0} transcripts processed.")
                # - k because you don't want the last k-mer to contain the newline
                for pos in range(len(line) - k):
                    curr_dint = line[pos: pos + k]
                    hk.add_key(curr_dint, [], kmer_dict)
                    # the positions are counted from the start of the interval, not the start of the transcript
                    kmer_dict[curr_dint].append((prev, pos + focal_pos))
    return kmer_dict

def pick_random_positions(polII_file, dints_fasta, outfile, kmer_dict, transcripts_dict, chrom_sizes = None):
    """
    Replace NET-seq reads with simulated reads that start with the same 4-mer.
    :param polII_file: NET-seq reads
    :param dints_fasta: FASTA file with only the -2:+2 4-mer of the NET-seq reads,
    has to be in the same order as the polII_file
    :param outfile: output file, same format as polII_file
    :param kmer_dict: dictionary of kmer positions as output by map_kmers_to_positions()
    :param transcripts_dict: dictionary with transcript IDs as values and their BED-format
    :param chrom_sizes: if specified, should be a dictionary with chromosome IDs as keys and
    their lengths as values. Will be used to make sure control reads "fit" on the
    chromosome.
    coordinates as values
    :return: None
    """
    with open(polII_file) as read_file, open(dints_fasta) as dints_file, open(outfile, "w") as write_file:
        reader = csv.reader(read_file, delimiter = "\t")
        writer = csv.writer(write_file, delimiter = "\t")
        counter = 0
        for line in dints_file:
            if line[0] != ">":
                counter = hk.update_counter(counter, 100000, "{0} lines processed.")
                # -1 because you don't want the newline
                curr_dint = line[:-1]
                curr_bed_line = next(reader)
                length = int(curr_bed_line[2]) - int(curr_bed_line[1])
                # sometimes, a read will appear super long because it has been mapped across
                # an intron
                # so you don't pick a control position if the new read doesn't "fit"
                # on the chromosome
                new_start = -5
                while new_start < 0:
                    random_pos = random.choice(kmer_dict[curr_dint])
                    new_chrom, new_start, new_end, new_strand = transpose_to_new_coords(transcripts_dict, random_pos, length, chrom_sizes = chrom_sizes)
                if new_start < 0:
                    print(curr_bed_line)
                    print(random_pos)
                    print("\n")
                new_line = curr_bed_line.copy()
                new_line[0] = new_chrom
                new_line[1] = str(new_start)
                new_line[2] = str(new_end)
                new_line[5] = new_strand
                writer.writerow(new_line)

def transpose_to_new_coords(transcripts_dict, random_pos, length, chrom_sizes = None):
    picked_trans = transcripts_dict[random_pos[0]]
    strand = picked_trans[5]
    new_chrom = picked_trans[0]
    if strand == "+":
        new_start = int(int(picked_trans[1]) + random_pos[1])
        new_end = int(new_start + length)
    else:
        new_end = int(int(picked_trans[2]) - random_pos[1])
        new_start = int(new_end - length)
    if chrom_sizes:
        if new_chrom not in chrom_sizes:
            new_start = -1
        elif new_end > chrom_sizes[new_chrom]:
            # this way, the generated random value will not be considered and a new one will
            # be picked instead
            new_start = -1
    return new_chrom, new_start, new_end, strand


def main():
    description = "Generate a NET-seq control set that would have the same distribution of -2:2 nucleotides" \
                  "as the true set."
    args = hk.parse_arguments(description, ["active_genes_file", "gtf", "PolII_file", "fasta", "outfile", "chrom_sizes"])
    active_genes_file, gtf, PolII_file, fasta, outfile, chrom_sizes = args.active_genes_file, args.gtf, args.PolII_file, args.fasta, args.outfile, args.chrom_sizes

    chrom_sizes = rw.read_many_fields(chrom_sizes, delimiter = "\t")
    chrom_sizes = hk.list_to_dict(chrom_sizes, 0, 1, intify=True)

    # get transcriptionally active genes and make a BED file with their coordinates
    print("Getting the coordinates of transcriptionally active genes...")
    trans_active_genes = rw.read_many_fields(active_genes_file, "\t")[1:]
    trans_active_genes = [i[3] for i in trans_active_genes]
    transcripts_file = "{0}_transcripts_all.bed".format(gtf[:-4])
    co.get_transcripts(gtf, transcripts_file)

    transcripts_dict = {}
    # this will be used for getting the k-mers in the transcripts
    filtered_transcripts_file_plus2 = "{0}_trans_act_only_plus3.bed".format(transcripts_file[:-4])
    # this will be used for filtering the reads
    filtered_transcripts_file = "{0}_trans_act_only.bed".format(transcripts_file[:-4])
    with open(filtered_transcripts_file, "w") as ft_file, open(transcripts_file) as t_file, open(filtered_transcripts_file_plus2, "w") as ft_file2:
        reader = csv.reader(t_file, delimiter="\t")
        writer = csv.writer(ft_file, delimiter="\t")
        writer2 = csv.writer(ft_file2, delimiter="\t")
        for line in reader:
            if line[3] in trans_active_genes:
                # if line[0][0] not in ["G", "K"]:
                #     line[0] = "chr{0}".format(line[0])
                writer.writerow(line)
                # this is because if a read falls at the first position, you will need to know the
                # preceding two bases. Same if it falls at the last position.
                line[1] = str((int(line[1])) - 3)
                line[2] = str((int(line[2])) + 3)
                writer2.writerow(line)
                transcripts_dict[line[3]] = line

    print("Filtering reads to the transcripts...")
    # filter reads to only ones that overlap these transcripts
    transcripts_PolII = "{0}_transcripts.bed".format(PolII_file[:-4])
    co.intersect_bed(PolII_file, filtered_transcripts_file, force_strand=True, output_file=transcripts_PolII)

    print("Extracting FASTA from the transcript coordinates...")
    # the genome FASTA is formatted as N rather than chrN
    filtered_transcripts_file_no_chr = "{0}_trans_act_only_plus3_no_chr.bed".format(transcripts_file[:-4])
    hk.run_process(["sed", "s/^chr//", filtered_transcripts_file_plus2], file_for_output=filtered_transcripts_file_no_chr)
    filtered_transcripts_fasta_no_chr = "{0}_trans_act_only_plus3.fasta".format(transcripts_file[:-4])
    hk.run_process(["bedtools", "getfasta", "-fi", fasta, "-bed", filtered_transcripts_file_no_chr, "-fo", filtered_transcripts_fasta_no_chr, "-s", "-name"])

    print("Mapping kmers to transcript positions...")
    kmer_dict = map_kmers_to_positions(filtered_transcripts_fasta_no_chr, k = 6, focal_pos = 3)

    print("Extracting the starting dinucleotide for each read...")
    starting_dints_PolII = "{0}_transcripts_starting_6mers.bed".format(PolII_file[:-4])
    starting_dints_PolII_fasta = "{0}_transcripts_starting_6mers.fasta".format(PolII_file[:-4])
    co.extend_intervals(transcripts_PolII, starting_dints_PolII, 3, 3, remove_chr=True)
    hk.run_process(["bedtools", "getfasta", "-fi", fasta, "-bed", starting_dints_PolII, "-fo", starting_dints_PolII_fasta, "-s"])

    print("Picking random control positions...")
    pick_random_positions(transcripts_PolII, starting_dints_PolII_fasta, outfile, kmer_dict, transcripts_dict, chrom_sizes = chrom_sizes)

    print("Making single nucleotide resolution file...")
    snr_file = "{0}_snr.bed".format(outfile[:-4])
    co.snr_bed(outfile, snr_file)

    print("Removing reads that overlap potential splice intermediate positions...")
    no_si_snr_file = "{0}_snr_no_si.bed".format(outfile[:-4])
    co.intersect_bed(snr_file, "data/Genomes/GTFs/dm6/dmel-all-r6.18_exon_ends_chr.gtf", force_strand=True, exclude = True, no_dups=False)

if __name__ == "__main__":
    main()