'''Functions to help process NGS, especially NET-seq data.'''

import coord_ops as co
import csv
import housekeeping as hk
import re

def analyze_cigar(bed_line, overhang = 0):
    """
    Take a line from a BED file that results from a -wb intersect of
    a BED file of reads with a BED file of 3' splice sites and determine
    whether the read corresponds to a spliced read, an unspliced read,
    or 'unknown'.
    :param bed_line: line from BED file
    :param overhang: minimum overhang to either side of intron
    :return: "S", "U" or None for the different read types.
    """
    cigar = bed_line[7]
    strand = bed_line[5]
    read_start = int(bed_line[1])
    read_end = int(bed_line[2])
    # position of the first nucleotide of the downstream exon
    three_prime = int(bed_line[18])
    if strand == "+":
        # if the read starts/ends at the first position of the exon
        # and doesn't reach into the intron
        if three_prime == read_start:
            return None
    elif strand == "-":
        if three_prime + 1 == read_end:
            return None
    # parse cigar
    # leave out first/last element because will be empty string
    letters = re.split("\d+", cigar)[1:]
    numbers = re.split("[A-Z]", cigar)[:-1]
    # calculate
    # check if the cigar matches a spliced read
    if "N" in letters:
        # check if the intron length is as expected:
        intron_index = letters.index("N")
        intron_length = int(numbers[intron_index])
        if bed_line[-3] == "." or intron_length == int(bed_line[-3]):
            # check if the intron position is correct
            # don't count nucleotides inserted in the read
            dist_to_intron = None
            if strand == "+":
                dist_to_intron_upstream = sum([int(i) for pos, i in enumerate(numbers[: intron_index]) if letters[: intron_index][pos] != "I"])
                ref_dist_to_intron_upstream = (three_prime - intron_length) - read_start
                # downstream distance to check overhangs later
                ref_dist_to_intron_downstream = read_end - three_prime
            elif strand == "-":
                dist_to_intron_upstream = sum([int(i) for pos, i in enumerate(numbers[intron_index + 1:]) if letters[intron_index + 1:][pos] != "I"])
                ref_dist_to_intron_upstream = read_end - (three_prime + intron_length + 1)
                ref_dist_to_intron_downstream = (three_prime + 1) - read_start
            if dist_to_intron_upstream == ref_dist_to_intron_upstream:
                if ref_dist_to_intron_downstream >= overhang:
                    if ref_dist_to_intron_upstream >= overhang:
                        return("S")
    else:
        if strand == "+":
            ref_dist_to_intron_downstream = read_end - three_prime
            if ref_dist_to_intron_downstream >= overhang:
                return("U")
        else:
            ref_dist_to_intron_downstream = (three_prime + 1) - read_start
            if ref_dist_to_intron_downstream >= overhang:
                return("U")

def check_intermediate_read(line, exons):
    '''
    Check if a BED file line corresponds to an intermediate read.
    :param line: line from BED file (NET-seq reads intersected with 3' splice sites)
    :param exons: a CDS dictionary
    :return: 'True' if the read is intermediate, 'False' otherwise
    '''
    trans_name = line[-4].split(".")[0]
    current_trans = exons[trans_name]
    if line[-2] == "+":
        exon_ends = [i[4] for i in current_trans]
        if int(line[2]) in exon_ends:
            return True
    else:
        # -1 to go from base 1 to base 0
        exon_ends = [int(i[3]) - 1 for i in current_trans]
        if int(line[1]) in exon_ends:
            return True
    return False

def check_intron_end_peak(line, exons):
    '''
    Check if an intersect BED file line maps to the intron end peak.
    :param line: line from BED file (NET-seq reads intersected with 3' splice sites)
    :param exons: a CDS dictionary
    :return: 'True' if the 5' end of the read maps to the intron end peak, 'False' otherwise
    '''
    trans_name = line[-4].split(".")[0]
    current_trans = exons[trans_name]
    if line[-2] == "+":
        # the condition is to filter out cases where the peak would overlap with the preceding exon
        # note that current_trans[pos] refers to the PRECEDING exon
        # -1 to go from base 1 to base 0
        intron_end_peaks = [(i[3] - 25 - 1, i[3] - 1 - 1) for pos, i in enumerate(current_trans[1:]) if current_trans[pos][4] <= (i[3] - 25)]
        read_start = int(line[1])
        for peak in intron_end_peaks:
            if read_start >= peak[0] and read_start < peak[1]:
                return True
    else:
        intron_end_peaks = [(i[4] + 1, i[4] + 25) for pos, i in enumerate(current_trans[1:]) if current_trans[pos][3] > (i[4] + 25)]
        read_end = int(line[2])
        for peak in intron_end_peaks:
            if read_end > peak[0] and read_end <= peak[1]:
                return True
    return False

def check_position_in_exon(bed_line, exons, filter_start = None, filter_end = None):
    """
    Check that the read ends in the exon downstream from the 3'ss that
    it overlaps with.
    :param bed_line: line from BED file (NET-seq reads intersected with 3'ss)
    :param exons: CDS dictionary
    :param filter_start: if specified, give False if read ends before this exonic position
    :param filter_end: if specified, give False if read ends after this exonic position
    :return: 'True' if read end in relevant exon, 'False' otherwise
    """
    intron_ID = bed_line[-4].split(".")
    trans_name = intron_ID[0]
    intron_ID = int(intron_ID[1])
    # the exon downstream from intron i is the exon i + 1
    curr_exon = exons[trans_name][intron_ID + 1]
    strand = bed_line[5]
    if strand == "+":
        read_end = int(bed_line[2])
        # converting to base 0
        exon_start = (curr_exon[3] - 1)
        if filter_start:
            exon_start = exon_start + filter_start
        if read_end > exon_start:
            exon_end = curr_exon[4]
            if filter_end:
                exon_end = exon_start + filter_end
            if read_end <= exon_end:
                return True
    else:
        read_end = int(bed_line[1])
        exon_end = (curr_exon[3] - 1)
        if filter_end:
            exon_end = curr_exon[4] - filter_end
        if read_end >= exon_end:
            exon_start = curr_exon[4]
            if filter_start:
                exon_start = exon_start - filter_start
            if read_end < exon_start:
                return True
    return False

def density_per_transcript(exon_file, polII_bed, out_file, bed_input=False):
    """
    Given a GTF file of exon coordinates and a BED file of Pol II positions
    (or other types of reads),
    calculate the read density per transcript and store in file.
    :param exon_file: GTF file with exon coordinates
    :param polII_bed: BED file with reads
    :param out_file: tab-delimited file with transcript IDs in columns one
    and densities in column two
    :param bed_input: if True, BED (rather than GTF) input is assumed
    :return: dictionary of read densities
    """
    # get read counts per exon
    hits_per_exon = co.intersect_bed(exon_file, polII_bed, force_strand=True,
                     no_dups=False, hit_count = True)
    trans_regex = re.compile("(?<=transcript_id \")[A-Za-z0-9]+")
    out_dict = {}
    for exon in hits_per_exon:
        # parse out the exon name
        if bed_input:
            trans_ID = exon[3].split(".")[0]
        else:
            trans_ID = re.search(trans_regex, exon[8]).group(0)
        if trans_ID not in out_dict:
            out_dict[trans_ID] = {"length": 0, "read_count": 0}
        if bed_input:
            out_dict[trans_ID]["length"] = out_dict[trans_ID]["length"] + \
                                           (int(exon[2]) - int(exon[1]))
        else:
            out_dict[trans_ID]["length"] = out_dict[trans_ID]["length"] + \
                                           (int(exon[4]) - int(exon[3]) + 1)
        out_dict[trans_ID]["read_count"] = out_dict[trans_ID]["read_count"] + int(exon[-1])

    with open(out_file, "w") as file:
        file.write("transcript\treads_per_nt\n")
        for trans in out_dict:
            out_dict[trans] =  out_dict[trans]["read_count"]/out_dict[trans]["length"]
            file.write("{0}\t{1}\n".format(trans, str(out_dict[trans])))
    return(out_dict)

def get_splice_dist(bed_line):
    """
    Calculate the splicing distance: the distance between the 3' splice site
    and the 3' end of the read.
    :param bed_line: line from BED file
    :return: splicing distance in nucleotides
    """
    strand = bed_line[5]
    three_prime = int(bed_line[18])
    if strand == "+":
        read_end = int(bed_line[2])
        return(read_end - three_prime)
    else:
        read_end = int(bed_line[1]) - 1
        return(three_prime - read_end)

def get_spliced_reads(reads_file, exon_junctions_file, spliced_bed, unspliced_bed, exons = False, overhang = 0, intron_end_peak = False, remove_chr = False, filter_start = None, filter_end = None):
    """
    Detect spliced/unspliced reads and write them into separate BED files.
    :param reads_file: BED file with reads
    :param exon_junctions_file: BED file with the first positions of exons
    :param spliced_bed: output BED file for spliced reads
    :param unspliced_bed: output BED file for unspliced reads
    :param exons: exon coordinate dictionary
    :param overhang: the splicing status of the reads will only be called if they're at least
    __overhang__ nucleotides away from the splice site.
    :param intron_end_peak: if a file name has been specified,
    reads that map to intronic sequence within nucleotides
    -25 to -2 before the start of the exon will be stored in a separate file
    :param remove_chr: if True, remove "chr" from the start of chromosome names
    :param filter_start: if specified, only report reads that end after this exonic nucleotide
    :param filter_end: if specified, only report reads that end before this exonic nucleotide
    :return: None
    """
    exon_junction_bed = "{0}_exon_junctions_intersect.bed".format(reads_file[:-4])
    co.intersect_bed(reads_file, exon_junctions_file, write_both=True,
                     output_file=exon_junction_bed,
                  force_strand=True, no_dups=False)
    file_size = hk.line_count(exon_junction_bed)
    if intron_end_peak:
        intron_end_peak_file = open(intron_end_peak, "w")
        intron_end_writer = csv.writer(intron_end_peak_file, delimiter = "\t")
    found_count = 0
    with open(exon_junction_bed) as file, open(spliced_bed, "w") as sfile, open(unspliced_bed, "w") as ufile:
        reader = csv.reader(file, delimiter="\t")
        spliced_writer = csv.writer(sfile, delimiter="\t")
        unspliced_writer = csv.writer(ufile, delimiter="\t")
        for pos, line in enumerate(reader):

            if pos % 100000 == 0:
                print("{0}/{1}".format(pos, file_size))
                print("Found {0} spliced reads.".format(found_count))
                print("\n")

            # reads that end at the last nucleotide of an exon
            intermediate_read = False
            if exons:
                intermediate_read = check_intermediate_read(line, exons)
            intron_name = line[20]

            if not intermediate_read:

                # check that it ends within the exon just downstream of
                # the 3' ss that is being analyzed

                in_dwns_exon = True
                if exons:
                    in_dwns_exon = check_position_in_exon(line, exons, filter_start = filter_start, filter_end = filter_end)

                if in_dwns_exon:

                    if remove_chr:
                        line[0] = line[0].lstrip("chr")
                    carry_on = True
                    if intron_end_peak:
                        intron_end_peak_boolean = check_intron_end_peak(line, exons)
                        if intron_end_peak_boolean:
                            intron_end_writer.writerow(line[:17])
                            carry_on = False
                    if carry_on:

                        # 'spliced', 'unspliced' or 'None' (=can't analyze)
                        read_type = analyze_cigar(line, overhang = overhang)

                        if read_type:
                            # :17 because we only want the read informatuion, not the splice junction
                            # information
                            if read_type == "S":
                                spliced_writer.writerow([str(i) for i in line[:17]])
                                found_count = found_count + 1
                            else:
                                unspliced_writer.writerow([str(i) for i in line[:17]])
    if intron_end_peak:
        intron_end_peak_file.close()