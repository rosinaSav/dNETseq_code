'''Identify spliced.unspliced reads and build a meta-profile for either.'''

import coord_ops as co
import housekeeping as hk
import NGS
import numpy as np
import read_and_write as rw

def make_dist_mat(distances, max_dist, lengths, scale_matrix = False):
    '''
    Make a matrix of spliced and unspliced read splicing distances,
    with introns in rows and positions in columns,
    and the number of spliced/unspliced reads ending at each position
    in each cell.
    :param distances: a dictionary with intron IDs as keys and lists of
    splicing distances as values
    :param max_dist: maximum distance to consider
    :param lengths: a dictionary with intron IDs as keys and downstream
    exon lengths as values
    :param scale_matrix: scale all values based on attached dictionary, where the keys
    are intron names and the values are the scaling factors
    :return: Sorted list of intron identifiers, distance matrix, and
    position of first read.
    '''
    mat = np.empty((len(distances), max_dist))
    mat[:] = np.nan
    first_read = []
    for pos, key in enumerate(sorted(distances.keys())):
        if lengths:
            curr_length = lengths[key]
        else:
            curr_length = max_dist
        mat[pos,:curr_length] = 0
        if scale_matrix:
            trans_name = key.split(".")[0]
            scaling_value = scale_matrix[trans_name]
        else:
            scaling_value = 1
        if len(distances[key]) > 0:
            first_read.append(min(distances[key]))
        for distance in distances[key]:
            if distance < max_dist and distance >= 0:
                mat[pos, distance] = mat[pos, distance] + (1/scaling_value)
    return(sorted(distances.keys()), mat, first_read)

def update_dist_dict(intron_name, dist_dict, splice_dist):
    '''
    Add a splicing distance to the correct intron in the splicing distance
    dictionary.
    :param intron_name: intron ID (dictionary key)
    :param dist_dict: distance dictionary
    :param splice_dist: distance between 3'ss and 3' read end
    :return: splicing distance dictionary (updated)
    '''
    if intron_name not in dist_dict:
        dist_dict[intron_name] = []
    dist_dict[intron_name].append(splice_dist)
    return(dist_dict)

def write_dist_mat(distances, max_dist, file_name, lengths, ID_file_name = None, first_spliced_file_name = None, scale_matrix = False):
    """
    Given a dictionary of splicing distances for different introns,
    turn them into a binary matrix and write to file.
    :param distances: Distance dictionary.
    :param max_dist: Maximum distance to consider.
    :param file_name: File in which to write the matrix.
    :param lengths: a dictionary with intron IDs as keys and downstream
    exon lengths as values
    :param ID_file_name: File in which to write the intron IDs.
    :param first_spliced_file_name: File in which to write the first position at
    which splicing is observed for each intron.
    :param scale_matrix: scale all values based on attached dictionary, where the keys
    are intron names and the values are the scaling factors
    :param filter_length: if specified, then instead of filtering transcripts by length based on the
    desired region size, they will be filtered based on this value instead
    :return: None
    """
    IDs, mat, first_spliced = make_dist_mat(distances, max_dist, lengths,
                                            scale_matrix = scale_matrix)
    np.savetxt(file_name, mat, delimiter = "\t")
    if ID_file_name:
        rw.write_list(IDs, ID_file_name)
    if first_spliced_file_name:
        rw.write_list(first_spliced, first_spliced_file_name)

def write_exon_starts(valid_junctions, file_name, exons, limit, add_chr = False, from_end = False, centre = False):
    """
    Given a list of intron IDs, write the coordinates of the first _limit_ bp
    of the downstream exon into a BED file.
    :param valid_junctions: list of intron IDs
    :param file_name: output file name
    :param exons: exons dictionary
    :param limit: how many nucleotides from the exon start to consider
    :param add_chr: if True, "chr" is prefixed to chromosome names
    :param from_end: if True, the last _limit_ bp are written to file instead
    :param centre: if True, the region will be centred on the exon boundary of interest rather
    than being in the beginning/at the end of the exon
    :return: None
    """
    plus = "+"
    if from_end:
        plus = "-"
    with open(file_name, "w") as out:
        for intron_name in valid_junctions:
            intron = intron_name.split(".")
            trans = intron[0]
            exon_number = int(intron[1]) + 1
            curr_exon = exons[trans][exon_number]
            chrom = curr_exon[0]
            if add_chr:
                chrom = "chr{0}".format(chrom)
            strand = curr_exon[6]
            if strand == plus:
                # convert to base 0
                start = curr_exon[3] - 1
                end = start + limit
                if centre:
                    start = int(start - limit/2)
                    end = int(end - limit/2)
            else:
                start = curr_exon[4] - limit
                end = curr_exon[4]
                if centre:
                    start = int(start + limit/2)
                    end = int(end + limit/2)
            out.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(chrom, start, end, intron_name, ".", strand))

def write_intron_starts(valid_junctions, file_name, exons, limit, add_chr = False):
    """
    Given a list of intron IDs, write the coordinates of the first _limit_ bp
    of the intron into a BED file.
    :param valid_junctions: list of intron IDs
    :param file_name: output file name
    :param exons: exons dictionary
    :param limit: how many nucleotides from the intron start to consider
    :param add_chr: if True, "chr" is prefixed to chromosome names
    :return: None
    """
    with open(file_name, "w") as out:
        for intron_name in valid_junctions:
            intron = intron_name.split(".")
            trans = intron[0]
            intron_number = int(intron[1])
            ups_exon = exons[trans][intron_number]
            chrom = ups_exon[0]
            if add_chr:
                chrom = "chr{0}".format(chrom)
            strand = ups_exon[6]
            if strand == "+":
                # convert to base 0
                start = ups_exon[4]
                end = start + limit
            else:
                start = ups_exon[3] - limit - 1
                end = ups_exon[3] - 1
            out.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(chrom, start, end, intron_name, ".", strand))

def write_read_lengths(input_file, output_file):
    '''
    Take a BED file of NGS reads (with the read sequence in column 12)
    and write and output file with the lengths of the reads.
    :param input_file: a BED file of reads with the read sequence in column 12 (base 1)
    :param output_file: a text file with the lengts of the reads, with one length per row
    :return: None
    '''
    with open(input_file) as bed_file, open(output_file, "w") as o_file:
        for line in bed_file:
            line = line.split("\t")
            length = len(line[11].rstrip("\n"))
            o_file.write("{0}\n".format(length))

def write_si_pos(valid_junctions, file_name, exons, add_chr=False, curr_exon=False):
    """
    Given a list of intron IDs, write the coordinates of the last nucleotide in
    the previous exon into a BED file.
    :param valid_junctions: list of intron IDs
    :param file_name: output file name
    :param exons: exons dictionary
    :param add_chr: if True, "chr" is prefixed to chromosome names
    :param curr_exon: if True, the last nucleotide of the CURRENT exon is written to file,
    not of the upstream one.
    :return: None
    """
    with open(file_name, "w") as out:
        for intron_name in valid_junctions:
            intron = intron_name.split(".")
            trans = intron[0]
            intron_number = int(intron[1])
            ups_exon = exons[trans][intron_number]
            if curr_exon:
                ups_exon = exons[trans][intron_number + 1]
            chrom = ups_exon[0]
            if add_chr:
                chrom = "chr{0}".format(chrom)
            strand = ups_exon[6]
            if strand == "+":
                # convert to base 0
                start = ups_exon[4] - 1
                end = ups_exon[4]
            else:
                start = ups_exon[3] - 1
                end = ups_exon[3]
            out.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(chrom, start, end, intron_name, ".", strand))


def main():
    description = "Record splicing distance."
    args = hk.parse_arguments(description, ["input_file", "gtf", "output_folder", "trans_active_file", "window_size", "intron_window_size", "outsuffix", "leave_terminal"], ints = [4, 5], flags = [7])
    input_file, gtf, output_folder, trans_active_file, window_size, intron_window_size, outsuffix, leave_terminal = args.input_file, args.gtf, args.output_folder, args.trans_active_file, args.window_size, args.intron_window_size, args.outsuffix, args.leave_terminal

    if outsuffix == "None":
        outsuffix = ""

    bare_input_path = input_file.split("/")[-1]
    bed = "{0}.bed".format(input_file[:-4])
    # hk.convert2bed(input_file, bed)

    # get descriptive stats of the reads
    length_file = "{0}/{1}_read_lengths.txt".format(output_folder, bare_input_path[:-4])
    write_read_lengths(bed, length_file)

    # read in CDS coordinates
    exons = rw.read_gtf(gtf, "CDS", gene=False)
    # only leave transcriptionally active genes (one isoform per gene)
    trans_active_genes = rw.read_many_fields(trans_active_file, "\t")[1:]
    # pull out the column with transcript IDs
    trans_active_genes = [i[3] for i in trans_active_genes]
    exons = {i: exons[i] for i in exons if i in trans_active_genes}
    terminal_suff = "_with_terminal"
    if not leave_terminal:
        # remove last exons
        exons = {i: exons[i][:-1] for i in exons}
        terminal_suff = ""
    # prepare exon-exon junctions
    exon_junctions_file = "{0}_exon_junctions{1}{2}.bed".format(gtf[:-4], outsuffix, terminal_suff)
    all_junctions = co.extract_3ss(exons, exon_junctions_file)

    out_bed = "{0}/{1}_first_{2}_bp{3}{4}.bed".format(output_folder, bare_input_path[:-4], window_size, outsuffix, terminal_suff)
    write_exon_starts(all_junctions, out_bed, exons, window_size, add_chr=True)
    out_bed_end = "{0}/{1}_last_{2}_bp{3}{4}.bed".format(output_folder, bare_input_path[:-4], window_size, outsuffix, terminal_suff)
    write_exon_starts(all_junctions, out_bed_end, exons, window_size, add_chr=True, from_end=True)
    intron_bed = "{0}/{1}_first_{2}_intronic_bp{3}{4}.bed".format(output_folder, bare_input_path[:-4], intron_window_size, outsuffix, terminal_suff)
    write_intron_starts(all_junctions, intron_bed, exons, intron_window_size, add_chr=True)
    out_bed = "{0}/{1}_first_centred_{2}_bp{3}{4}.bed".format(output_folder, bare_input_path[:-4], window_size, outsuffix, terminal_suff)
    write_exon_starts(all_junctions, out_bed, exons, window_size, add_chr=True, centre=True)
    out_bed_end = "{0}/{1}_last_centred_{2}_bp{3}{4}.bed".format(output_folder, bare_input_path[:-4], window_size, outsuffix, terminal_suff)
    write_exon_starts(all_junctions, out_bed_end, exons, window_size, add_chr=True, from_end=True, centre=True)
    out_bed_si = "{0}/{1}_si_pos{2}{3}.bed".format(output_folder, bare_input_path[:-4], outsuffix, terminal_suff)
    write_si_pos(all_junctions, out_bed_si, exons, add_chr=True)
    out_bed_si_current = "{0}/{1}_si_pos_current{2}{3}.bed".format(output_folder, bare_input_path[:-4], outsuffix, terminal_suff)
    write_si_pos(all_junctions, out_bed_si_current, exons, add_chr=True, curr_exon=True)
    # check which junctions are associated with a splicing intermediate read
    snr_bed = "{0}_snr.bed".format(bed[:-4])
    co.snr_bed(bed, snr_bed)
    si_counts_bed = "{0}/{1}_si_counts{2}{3}.bed".format(output_folder, bare_input_path[:-4], outsuffix, terminal_suff)
    co.intersect_bed(out_bed_si, snr_bed, force_strand=True, hit_count=True, no_dups=False, output_file=si_counts_bed)
    si_counts_current_bed = "{0}/{1}_si_counts_current{2}{3}.bed".format(output_folder, bare_input_path[:-4], outsuffix, terminal_suff)
    co.intersect_bed(out_bed_si_current, snr_bed, force_strand=True, hit_count=True, no_dups=False, output_file=si_counts_current_bed)

    # filter out reads that don't overlap exon-exon junctions
    exon_junction_bed = "{0}_exon_junctions{1}{2}.bed".format(input_file[:-4], outsuffix, terminal_suff)
    co.intersect_bed(bed, exon_junctions_file, write_both=True,
                     output_file=exon_junction_bed,
                  force_strand=True, no_dups=False)

    spliced_bed = "{0}_spliced{1}{2}.bed".format(input_file[:-4], outsuffix, terminal_suff)
    unspliced_bed = "{0}_unspliced{1}{2}.bed".format(input_file[:-4], outsuffix, terminal_suff)
    sr_distances = {}
    ur_distances = {}
    found_count = 0
    file_size = hk.line_count(exon_junction_bed)

    # will store all the intron names for which there are
    # either spliced or unspliced reads
    valid_junctions = []
    with open(exon_junction_bed) as file, open(spliced_bed, "w") as sfile, open(unspliced_bed, "w") as ufile:
        for pos, line in enumerate(file):

            if pos % 100000 == 0:
                print("{0}/{1}".format(pos, file_size))
                print("Found {0} spliced reads.".format(found_count))
                print("\n")

            line = line.split("\t")

            # reads that end at the last nucleotide of an exon
            intermediate_read = NGS.check_intermediate_read(line, exons)
            intron_name = line[20]

            if not intermediate_read:

                # check that it ends within the exon just downstream of
                # the 3' ss that is being analyzed

                in_dwns_exon = NGS.check_position_in_exon(line, exons)

                if in_dwns_exon:

                    # 'spliced', 'unspliced' or 'None' (=can't analyze)
                    read_type = NGS.analyze_cigar(line, overhang = 5)

                    if read_type:
                        if intron_name not in valid_junctions:
                            valid_junctions.append(intron_name)
                        splice_dist = NGS.get_splice_dist(line)
                        if read_type == "S":
                            sfile.write("\t".join([str(i) for i in line]))
                            found_count = found_count + 1
                            sr_distances = update_dist_dict(intron_name, sr_distances, splice_dist)
                        else:
                            ufile.write("\t".join([str(i) for i in line]))
                            ur_distances = update_dist_dict(intron_name, ur_distances, splice_dist)

    print("Proportion of spliced reads: {0}.".format(found_count/(pos + 1)))

    # for each valid junction, calculate the length of the exonic sequence
    # afterwards, so that you wouldn't consider intronic sequence in the distance
    # matrix
    lengths_dict = co.get_lengths(exons, valid_junctions)

    write_dist_mat(sr_distances, window_size,
                   "{0}/{1}_spliced_read_distances_{2}{3}{4}.txt".format(output_folder, bare_input_path[:-4], window_size, outsuffix, terminal_suff),
                   lengths_dict,
                   "{0}/{1}_spliced_read_{2}_intron_names{3}{4}.txt".format(output_folder, bare_input_path[:-4], window_size, outsuffix, terminal_suff),
                   "{0}/{1}_spliced_read_first_spliced{2}{3}.txt".format(output_folder, bare_input_path[:-4], outsuffix, terminal_suff))

    write_dist_mat(ur_distances, window_size,
                   "{0}/{1}_unspliced_read_distances_{2}{3}{4}.txt".format(output_folder, bare_input_path[:-4], window_size, outsuffix, terminal_suff),
                   lengths_dict,
                   "{0}/{1}_unspliced_read_{2}_intron_names{3}{4}.txt".format(output_folder, bare_input_path[:-4], window_size, outsuffix, terminal_suff),
                   "{0}/{1}_unspliced_read_first_unspliced{2}{3}.txt".format(output_folder, bare_input_path[:-4], outsuffix, terminal_suff))

if __name__ == "__main__":
    main()