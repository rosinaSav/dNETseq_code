'''Make a file with various gene architecture statistics for a set of introns.'''

import coord_ops as co
import housekeeping as hk
import NGS
import nucleotide_comp as nc
import numpy as np
from pyfaidx import Fasta
import read_and_write as rw

def add_to_array(out_array, new_dict, header):
    """
    Append another column to output array.
    :param out_array: Output array.
    :param new_dict: Dictionary with the values that are going to make up the column.
    :param header: Column name.
    :return: Updated output array.
    """
    new_list = [header]
    for junction in out_array[1:,0]:
        new_list.append(new_dict[junction])
    new_list = np.array(new_list)
    new_list.shape = (len(new_list), 1)
    out_array = np.hstack((out_array, new_list))
    return out_array

def get_dens_per_trans(truncated_exons_file, polII_bed, dens_per_trans_file, valid_junctions):
    """
    For each transcript that overlaps one of a set of pre-defined splice junctions,
    calculate Pol II read density per transcript.
    :param truncated_exons_file: BED file of exon coordinates, with the last nucleotide
    of the exon truncated so you wouldn't consider splicing intermediates.
    :param polII_bed: BED file of Pol II coordinates
    :param dens_per_trans_file: output file
    :param valid_junctions: list of valid 3' splice site IDs
    :return: dictionary with junction IDs as keys and transcript PolII densities as values
    """
    dens_per_trans = NGS.density_per_transcript(truncated_exons_file, polII_bed, dens_per_trans_file)
    dens_per_trans = hk.list_to_dict(dens_per_trans, 0, 1, floatify=True)
    dens_per_trans_junctions = {}
    for junction in valid_junctions:
        trans = junction.split(".")[0]
        dens_per_trans_junctions[junction] = dens_per_trans[trans]
    return(dens_per_trans_junctions)

def get_exon_GC(exons, exon_rank_start, genome):
    """
    For a set of 3' splice sites, calculate the GC content of the downstream exon.
    :param exons: dictionary of exon GTF lines, with transcript IDs as keys
    :param exon_rank_start: dictionary of downstream exon ranks, with junction IDs as keys
    :param genome: a pyfaidx Fasta object corresponding to the correct genome
    :return: dictionary of exon GC contents
    """
    out_dict = {}
    for junction in exon_rank_start:
        trans = junction.split(".")[0]
        curr_exon = exons[trans][exon_rank_start[junction]]
        seq = co.get_sequence(curr_exon, genome, impose_strand = True, strip_chr = True)
        GC = nc.get_GC(seq)
        out_dict[junction] = GC
    return(out_dict)

def get_exon_GC4(CDSs, exons, exon_rank_start, genome):
    """
    For a set of 3' splice sites, calculate the GC4 content of the downstream exon.
    :param CDSs: dictionary of exon CDS lines
    :param exons: dictionary of exon GTF lines, with transcript IDs as keys
    :param exon_rank_start: dictionary of downstream exon ranks, with junction IDs as keys
    :param genome: a pyfaidx Fasta object corresponding to the correct genome
    :return: dictionary of exon GC4 contents
    """
    out_dict = {}
    for junction in exon_rank_start:
        trans = junction.split(".")[0]
        curr_exon = exons[trans][exon_rank_start[junction]]
        curr_CDSs = CDSs[trans]
        # if it's on the plus strand, the starts have to match (the end might not match if the end
        # of the exon is non-coding
        # on the minus strand, the ends have to match
        if curr_exon[6] == "+":
            curr_CDS = [i for i in curr_CDSs if i[3] == curr_exon[3]][0]
        else:
            curr_CDS = [i for i in curr_CDSs if i[4] == curr_exon[4]][0]
        seq = co.get_sequence(curr_CDS, genome, impose_strand = True, strip_chr = True)
        GC4 = nc.get_GC4(seq, int(curr_CDS[7]))
        out_dict[junction] = GC4
    return(out_dict)

def get_upstream_intron_GC(exons, exon_rank_start, genome):
    """
    For a set of 3' splice sites, calculate the GC content of the upstream intron.
    :param exons: dictionary of exon GTF lines, with transcript IDs as keys
    :param exon_rank_start: dictionary of downstream exon ranks, with junction IDs as keys
    :param genome: a pyfaidx Fasta object corresponding to the correct genome
    :return: dictionary of upstream intron GC contents
    """
    introns = co.get_introns(exons)
    out_dict = {}
    for junction in exon_rank_start:
        trans = junction.split(".")[0]
        # -1 because the rank of the intron is 1 less than the rank of the downstream exon
        curr_intron = introns[trans][exon_rank_start[junction] - 1]
        seq = co.get_sequence(curr_intron, genome, impose_strand = True, strip_chr = True)
        GC = nc.get_GC(seq)
        out_dict[junction] = GC
    return(out_dict)

def main():
    description = "Aggregate various statistics on the splicing events you're studying."
    args = hk.parse_arguments(description, ["gtf", "polII_bed", "exon_start_coords", "truncated_exons_file", "genome_file", "output_file"])
    gtf, polII_bed, exon_start_coords, truncated_exons_file, genome_file, output_file = args.gtf, args.polII_bed, args.exon_start_coords, args.truncated_exons_file, args.genome_file, args.output_file

    CDSs = rw.read_gtf(gtf, "CDS", gene=False)
    exons = rw.read_gtf(gtf, "exon", gene=False)
    exon_starts = rw.read_many_fields(exon_start_coords, skip_header = False, delimiter = "\t")
    exon_starts = {i[3]: i for i in exon_starts}
    out_array = np.array(sorted(exon_starts.keys()), dtype="str")
    out_array.shape = (len(exon_starts.keys()), 1)
    out_array = np.vstack((["junction"], out_array))

    #1. exon size
    curr_dict = co.get_lengths(CDSs, exon_starts.keys())
    out_array = add_to_array(out_array, curr_dict, "exon_size")
    print("Exon size done.")

    #2. exon number
    curr_dict = co.get_exon_number(exons, exon_starts.keys())
    out_array = add_to_array(out_array, curr_dict, "exon_number")
    print("Exon number done.")

    #3. exon rank (from start and end)
    exon_rank_start, exon_rank_end = co.get_exon_rank(exons, exon_starts)
    out_array = add_to_array(out_array, exon_rank_start, "exon_rank_from_start")
    out_array = add_to_array(out_array, exon_rank_end, "exon_rank_from_end")
    print("Exon rank done.")

    #4. upstream intron size
    curr_dict = co.get_upstream_intron_size(exons, exon_rank_start)
    out_array = add_to_array(out_array, curr_dict, "upstream_intron_size")
    curr_dict = co.get_upstream_intron_size(exons, exon_rank_start, downstream=True)
    out_array = add_to_array(out_array, curr_dict, "downstream_intron_size")
    print("Intron size done.")

    if truncated_exons_file != "None":

        #5. Pol II density per transcript
        dens_per_trans_file = "{0}_dens_per_trans.txt".format(polII_bed[:-4])
        dens_per_trans_junctions = get_dens_per_trans(truncated_exons_file, polII_bed, dens_per_trans_file, out_array[1:,0])
        out_array = add_to_array(out_array, dens_per_trans_junctions, "polII_dens_per_trans")
        print("Pol II density done.")

    #6. exon GC4 and GC content
    genome = Fasta(genome_file)
    curr_dict = get_exon_GC4(CDSs, exons, exon_rank_start, genome)
    out_array = add_to_array(out_array, curr_dict, "exon_GC4")
    curr_dict = get_exon_GC(exons, exon_rank_start, genome)
    out_array = add_to_array(out_array, curr_dict, "exon_GC")
    print("Exon GC done.")

    #7. upstream intron GC content
    curr_dict = get_upstream_intron_GC(exons, exon_rank_start, genome)
    out_array = add_to_array(out_array, curr_dict, "upstream_intron_GC")
    print("Intron GC done.")

    #8. splice site strength
    curr_dict = nc.get_ss_strength(exons, genome_file, upstream = True, five = True, exonic = 3, intronic = 6)
    out_array = add_to_array(out_array, curr_dict, "upstream_5ss_strength")
    curr_dict = nc.get_ss_strength(exons, genome_file, upstream = True, five = False, exonic = 3, intronic = 20)
    out_array = add_to_array(out_array, curr_dict, "upstream_3ss_strength")
    curr_dict = nc.get_ss_strength(exons, genome_file, upstream = False, five = True, exonic = 3, intronic = 6)
    out_array = add_to_array(out_array, curr_dict, "downstream_5ss_strength")
    print("Splice site strength done.")

    with open(output_file, "w") as file:
        for line in range(0, out_array.shape[0]):
            line = out_array[line,:]
            line = "\t".join([str(i) for i in line])
            file.write(line)
            file.write("\n")

if __name__ == "__main__":
    main()