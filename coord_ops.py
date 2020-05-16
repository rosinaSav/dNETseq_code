'''
Functions for various operations on genomic coordinates.
'''

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
import copy
import csv
import housekeeping as hk
import numpy as np
import random
import read_and_write as rw

def bed4_to_bed6(inbed, outbed):
    """
    Given a BED4 file, convert it to BED6, using the coordinates to make the name
    and with a "." in the score field.
    :param inbed: BED4
    :param outbed: BED6
    :return: None
    """
    with open(inbed) as infile, open(outbed, "w") as outfile:
        reader = csv.reader(infile, delimiter = "\t")
        writer = csv.writer(outfile, delimiter = "\t")
        for line in reader:
            new_line = [line[0], line[1], line[2], "{0}_{1}_{2}_{3}".format(line[0], line[1], line[2], line[3]), ".", line[3]]
            writer.writerow(new_line)

def extend_intervals(input_file, output_file, left_shift, right_shift, remove_chr=False, add_chr=False, names_file=None, three_prime = False):
    """
    Given a BED file, make a new BED file with intervals that start _left_shift_ nt upstream
    of the interval starts in the original file and end _right_shift_ nt to the right.
    Note that for intervals on the negative strand, right and left will be
    reversed.
    :param input_file: input BED file
    :param output_file: output BED file
    :param left_shift: distance between old and new interval start
    :param right_shift: distance between old interval start and new interval end
    :param remove_chr: if True, remove "chr" from the chromosome name in the output
    :param add_chr: if True, prefix "chr" to the chromosome name in the output
    :param names_file: if specified, then reads will only be processed
    if the ID is in the specified file
    :param three_prime: if True, extend around interval ends instead
    :return: None
    """
    remove_counter = 0
    names = []
    if names_file:
        names = rw.read_as_string(names_file).split("\n")
    plus = "+"
    if three_prime:
        plus = "-"
        temp_left = left_shift
        left_shift = right_shift
        right_shift = temp_left
    with open(input_file) as bed, open(output_file, "w") as out_bed:
        reader = csv.reader(bed, delimiter="\t")
        writer = csv.writer(out_bed, delimiter="\t")
        for line in reader:
            if (not names_file) or (line[3] in names):
                if names_file:
                    names.remove(line[3])
                template = line.copy()
            # make a BED interval starting _left_shift_ nt before the 5' end and ending _right_shift_ nt after it.
                if len(line) >= 6:
                    template = line[:6]
                    if line[5] == plus:
                        template[1] = int(line[1]) - left_shift
                        template[2] = int(line[1]) + right_shift
                    else:
                        template[1] = int(line[2]) - right_shift
                        template[2] = int(line[2]) + left_shift
                # write the interval into a BED file, ignoring cases where the
                # read is so close to the start of the chromosome that you end up with a
                # negative coordinate
                    if template[1] >= 0:
                        if remove_chr:
                            template[0] = template[0].lstrip("chr")
                        if add_chr:
                            template[0] = "chr{0}".format(template[0])
                        writer.writerow(template)
                    else:
                        remove_counter = remove_counter + 1
    print("Removed because would have exceeded the chromosome: {0}.".format(remove_counter))

def extract_3ss(exons, bed):
    '''
    Given a dictionary of exons, extract the coordinates of the 3'ss
    and write to a BED file, along with the distance to the 5'ss.
    '''
    all_junctions = []
    with open(bed, "w") as file:
        for gene in sorted(exons.keys()):
            curr_exons = exons[gene]
            #ignore single-exon genes
            if len(curr_exons) > 1:
                # the list index in which 5' and 3' coordinates of
                # exons can be found
                start = 3
                end = 4
                strand = curr_exons[0][6]
                if strand not in ["+", "-"]:
                    raise Exception("Unexpected strand information: {0}".format(strand))
                if strand == "-":
                    start = 4
                    end = 3
                # because we want a pretty BED file in the end
                chrom = "chr" + curr_exons[0][0]

                for pos, exon in enumerate(curr_exons[:-1]):
                    # the adding and subtracting one positions everything to
                    # the first and last nucleotide of the intron for the plus strand
                    # but messes everything up for the minus strand
                    # so that will be corrected when calculating the distance
                    five_prime = curr_exons[pos][end] + 1
                    intron_end = curr_exons[pos + 1][start] - 1
                    if strand == "+":
                        distance = intron_end - five_prime + 1
                    else:
                        distance = five_prime - intron_end - 3
                    # also convert to base 0 for BED
                    # note that what is written down is the first
                    # base of the downstream exon
                    bed_line = [chrom, intron_end,
                                  intron_end + 1,
                                  "{0}.{1}".format(gene, pos),
                                  distance, strand]
                    file.write("\t".join([str(i) for i in bed_line]))
                    file.write("\n")
                    all_junctions.append("{0}.{1}".format(gene, pos))
    return(all_junctions)

def extract_exon_junctions(exons, bed):
    '''
    Given a dictionary of exons, extract the coordinates of the junctions and write to .bed
    '''
    with open(bed, "w") as file:
        for trans in sorted(exons.keys()):
            curr_exons = exons[trans]
            if len(curr_exons) > 1:
                strand = curr_exons[0][6]
                chrom = "chr" + curr_exons[0][0]
                first = 5
                second = 3
                if strand == "-":
                    curr_exons.reverse()
                    first = 3
                    second = 5
                for pos, exon in enumerate(curr_exons[:-1]):
                    intron_start = [chrom, exon[4] - 1, exon[4],
                                "{0}.{1}.{2}".format(trans, pos, first),
                                    ".", strand]
                    file.write("\t".join([str(i) for i in intron_start]))
                    file.write("\n")
                    intron_end = [chrom, curr_exons[pos + 1][3] - 1,
                                  curr_exons[pos + 1][3],
                                  "{0}.{1}.{2}".format(trans, pos, second),
                                  ".", strand]
                    file.write("\t".join([str(i) for i in intron_end]))
                    file.write("\n")

def get_coverage(regions_file, reads_file, output_file_name):
    """
    Given a BED file with regions and a BED/BAM file of reads,
    count how many reads overlap each region and output in new BED file.
    :param regions_file: BED file
    :param reads_file: BED/BAM file
    :param output_file_name: name for output BED file
    :return: None
    """
    hk.run_process(["bedtools", "coverage", "-a", regions_file, "-b",
                    reads_file, "-counts", "-s"], file_for_output=output_file_name)

def get_exon_and_intron_coords(exons, with_ups_intron=False):
    """
    Given a set of GTF lines corresponding to the exons in a transcript,
    return a list of tuples with the 0-BASED start and end coordinates of the
    exons and introns in the trasncript.
    :param exons: list of GTF lines (list of lists)
    :param with_ups_intron: if True, group introns together with the downstream exon
    :return: list of tuples
    """
    if exons[0][6] == "-":
        exons.reverse()
    out_list = []
    if with_ups_intron:
        if exons[0][6] == "+":
            out_list.append((exons[0][3] - 1, exons[0][4]))
            for pos, exon in enumerate(exons[1:]):
                out_list.append((exons[pos][4], exon[4]))
        else:
            for pos, exon in enumerate(exons[:-1]):
                out_list.append((exon[3] - 1, exons[pos + 1][3] - 1))
            out_list.append((exons[-1][3] - 1, exons[-1][4]))
    else:
        for pos, exon in enumerate(exons[:-1]):
            out_list.append((exon[3] - 1, exon[4]))
            out_list.append((exon[4], exons[pos + 1][3] - 1))
        out_list.append((exons[-1][3] - 1, exons[-1][4]))
    return(out_list)

def get_exon_number(exons, valid_junctions):
    """
    Given a dictionary of exons and a list of exon-exon junctions to consider,
    calculate the exon numbers of the relevant transcripts.
    :param exons: dictionary with transcript IDs as keys and a list of exon lines from
    a GTF file as values
    :param valid_junctions: list of valid junction IDs
    :return: output dictionary
    """
    out_dict = {}
    for junction in valid_junctions:
        trans = junction.split(".")[0]
        out_dict[junction] = len(exons[trans])
    return(out_dict)

def get_exon_rank(exons, exon_starts):
    """
    Given a dictionary of all exons, grouped by transcripts,
    and a second dictionary where exon-exon junctions IDs are
    the keys and exon start coordinates are the values,
    calculate the exon rank of each exon.
    :param exons: dictionary with transcript IDs as keys and a list of exon lines from
    a GTF file as values
    :param exon_starts: a dictionary where exon-exon junctions IDs are the keys and exon
    start coordinates are the values
    :return: output dictionary for ranks counted from the transcript starts and another
    dictionary for ranks counted from the transcript ends
    """
    out_dict_start = {}
    out_dict_end = {}
    for junction in exon_starts:
        trans = junction.split(".")[0]
        exon_start = exon_starts[junction]
        if exon_start[5] == "+":
            start = int(exon_start[1]) + 1
            start_index = 3
        else:
            start = int(exon_start[2])
            start_index = 4
        rank_start = [pos for pos, i in enumerate(exons[trans]) if i[start_index] == start]
        if len(rank_start) == 1:
            rank_start = rank_start[0]
        else:
            raise Exception("Found more/fewer than exactly one exon!")
        rank_end = len(exons[trans]) - rank_start - 1
        out_dict_start[junction] = rank_start
        out_dict_end[junction] = rank_end
    return out_dict_start, out_dict_end

def get_flanking_intron_sizes(exons):
    """
    Given a dictionary with transcript IDs as keys and GTF lines as values,
    calculate the sizes of the upstream and downstream flanking introns
    of each exon.
    :param exons: dictionary with transcript IDs as keys and GTF lines of exon coordinates as values
    :return: a dictionary of up- and downstream intron sizes
    """
    intron_sizes = {}
    for trans in exons:
        #single-exon genes have no introns
        if len(exons[trans]) < 2:
            intron_sizes[trans] = None
        else:
            strands = set([i[6] for i in exons[trans]])
            # this is to detect trans-splicing and similar weirdness
            if len(strands) > 1:
                intron_sizes[trans] = {}
                intron_sizes[trans]["upstream"] = "trans_spliced"
                intron_sizes[trans]["downstream"] = "trans_spliced"
            else:
                intron_sizes[trans] = {}
                if exons[trans][0][6] == "+":
                    current_intron_sizes = [exons[trans][i+1][3] - exons[trans][i][4] - 1 for i in range(len(exons[trans]) - 1)]
                elif exons[trans][0][6] == "-":
                    current_intron_sizes = [exons[trans][i][3] - exons[trans][i+1][4] - 1 for i in range(len(exons[trans]) - 1)]
                #sanity check
                if len([i for i in current_intron_sizes if i < 1]) != 0:
                    print(exons[trans])
                    print(current_intron_sizes)
                    raise Exception("Returned negative intron size!")
                #the upstream and downstream flanking introns are the same, except that you set either the first or the last one to None
                intron_sizes[trans]["upstream"] = [None]
                intron_sizes[trans]["upstream"].extend(current_intron_sizes)
                intron_sizes[trans]["downstream"] = current_intron_sizes
                intron_sizes[trans]["downstream"].append(None)
    return(intron_sizes)

def get_introns(exons, upstream = None, downstream = None):
    '''
    Infer intron coordinates from a set of exon coordinates.
    With the upstream and downstream options, you can only get the specified
    number of basepairs in the intron (either immediately downstream
    or upstream from an exon) rather than getting the whole exon.
    '''
    # if you've been handed a list with the coordinates for a single transcript
    # rather than a dictionary
    list_input = False
    if type(exons) == list:
        exons = {"dummy_trans": exons}
        list_input = True
    if upstream and downstream:
        print("You cannot restrict both the upstream and downstream distance to the exon!")
        raise Exception
    introns = {}
    for trans in exons:
        curr_exons = copy.deepcopy(exons[trans])
        if len(curr_exons) > 1:
            upstream_local = upstream
            downstream_local = downstream
            introns[trans] = []
            #if the transcript is on the antisense strand, downstream becomes upstream
            #and the other way around
            if curr_exons[0][6] == "-":
                upstream_local = downstream
                downstream_local = upstream
                curr_exons.reverse()
            #determine the intron start and end coordinates
            starts = [i[4] + 1 for i in curr_exons[:-1]]
            ends = [i[3] - 1 for i in curr_exons[1:]]
            #trim coordinates if necessary
            if upstream_local:
                starts = [i - upstream_local + 1 for i in ends]
            if downstream_local:
                ends = [i + downstream_local - 1 for i in starts]
            coords = zip(starts, ends)
            template = curr_exons[0]
            #format intron coordinates
            for coord in coords:
                current_intron = template.copy()
                current_intron[2] = "intron"
                current_intron[3] = coord[0]
                current_intron[4] = coord[1]
                introns[trans].append(current_intron)
            if current_intron[6] == "-":
                introns[trans].reverse()
    if list_input:
        try:
            introns = introns["dummy_trans"]
        # this will happen if it's a single-exon gene
        except KeyError:
            introns = None
    return(introns)

def get_lengths(exons, valid_junctions, intronic=False):
    """
    Make a dictionary with intron IDs as keys and the lengths of the downstream exons
    as values.
    :param exons: a dictionary with transcript IDs as keys and a list of lists of
    exon coordinates (in GTF format) as values
    :param valid_junctions: a list of valid intron IDs
    :param intronic: return intron lengths instead
    :return: dictionary of downstream exon lengths
    """
    # when returning exon lengths, you want to get the downstream exon lengths
    # so you will be adding 1 to the intron position to get the exon position
    # but you don't need to do that when returning intron lengths
    to_add = 1
    if intronic:
        exons = get_introns(exons)
        to_add = 0
    out_dict = {}
    # the intron IDs are formed as 'trans_ID.intron_pos', with the intron
    # positions in base 0
    for intron in valid_junctions:
        intron_split = intron.split(".")
        curr_exon = exons[intron_split[0]][int(intron_split[1]) + to_add]
        length = curr_exon[4] - curr_exon[3] + 1
        out_dict[intron] = length
    return out_dict

def get_sequence(coords, genome_seq, impose_strand = False, bed_input = False, strip_chr = False, add_chr = False):
    '''
    Get the sequence corresponding to a GTF interval from a genome sequence.
    If a list of features is given, the resulting sequences will be concatenated in the order
    in which they appear
    in the list (reversed order if impose_strand is True and the feature is on the - strand).
    OPTIONS
    impose_strand: if True, the reverse complement of the genome sequence will be returned for
    features on the - strand.
    False by default.
    '''
    if bed_input:
        strand_col = 5
    else:
        strand_col = 6
    #check whether only a single feature or a list of features was supplied. Convert to list if it's the former.
    if type(coords[0]) == str:
        coords = [coords]
    complete_sequence = []
    if impose_strand and coords[0][strand_col] == "-":
        coords = list(reversed(coords))
    for i in coords:
        chr_name = i[0]
        if strip_chr:
            chr_name = chr_name.lstrip("chr")
        if add_chr:
            chr_name = "chr{0}".format(chr_name)
        #pyfaidx works with 0-based coordinates
        try:
            if bed_input:
                sequence = genome_seq[chr_name][(int(i[1])):(int(i[2]))]
            else:
                sequence = genome_seq[chr_name][(i[3]-1):(i[4])]
        except KeyError:
            print("{0} not in genome fasta!".format(chr_name))
            return(None)
        complete_sequence.append(str(sequence).upper())
    sequence = "".join(complete_sequence)
    if impose_strand:
        if coords[0][strand_col] == "-":
            sequence = Seq(sequence, IUPAC.unambiguous_dna)
            sequence = sequence.reverse_complement()
            sequence = str(sequence)
    return(sequence)

def get_transcripts(gtf, out_file, add_chr = False):
    """
    Given a GTF file that has exon coordinates (among others),
    make an output BED file with transcript coordinates.
    :param gtf: input GTF file
    :param out_file: output BED file name
    :param add_chr: if True, prefix "chr" to chromosome names
    :return: None
    """
    exons = rw.read_gtf(gtf, "exon", gene = False)
    with open(out_file, "w") as file:
        out_writer = csv.writer(file, delimiter = "\t")
        for trans in sorted(list(exons.keys())):
            starts = [i[3] for i in exons[trans]]
            ends = [i[4] for i in exons[trans]]
            # just any exon to get those fields that will be the same for all of them
            template = exons[trans][0]
            if add_chr:
                chrom = "chr{0}".format(template[0])
            else:
                chrom = template[0]
            # convert to BED
            to_write = [chrom, min(starts) - 1, max(ends), trans, ".", template[6]]
            out_writer.writerow(to_write)


def get_upstream_intron_size(exons, exon_ranks, downstream = False):
    """
    Given a dictionary of exon coordinates and a list of exon start coordinates (so you'd
    know which exon to consider), calculate the sizes of upstream introns.
    :param exons: dictionary with transcript IDs as keys and lists of GTF lines of exon
    coordinates as values
    :param exon_ranks: dictionary with junction IDs as keys and downstream exon ranks from start
    as values
    :param downstream: if True, get downstream intron size instead
    :return: dictionary with junction IDs as keys and upstream intron sizes as values
    """
    upstream = "upstream"
    if downstream:
        upstream = "downstream"
    flanking_intron_sizes = get_flanking_intron_sizes(exons)
    out_dict = {}
    for junction in exon_ranks:
        trans = junction.split(".")[0]
        if flanking_intron_sizes[trans][upstream] == "trans_spliced":
            out_dict[junction] = 0
        else:
            out_dict[junction] = flanking_intron_sizes[trans][upstream][exon_ranks[junction]]
    return(out_dict)

def intersect_bed(bed_file1, bed_file2, use_bedops = False, overlap = False, overlap_rec = False, write_both = False, sort = False, output_file = None,
                             force_strand = False, force_opposite_strand = False, no_name_check = False, no_dups = True, chrom = None, intersect = False, hit_count = False, bed_path = None, intersect_bam=None,
                  write_zero = False, write_bed = False, exclude = False):
    '''Use bedtools/bedops to intersect coordinates from two bed files.
    Return those lines in bed file 1 that overlap with intervals in bed file 2.
    OPTIONS
    output_file: write output to this file
    use_bedops: use bedops rather than bedtools. Certain options are only valid with one of the two, see below.
    overlap: minimum oxverlap required as a fraction of the intervals in bed file 1 (EX: 0.8 means that the
    overlap has to be at least 80% of the intervals in bed file 1).
    overlap_rec: require that the overlap as a fraction of the intervals in file 2 be at least as high as
    the threshold indicated in -f.
    write_both: if True, return not only the interval from bed file 1 but, tagged onto the end, also the
    interval from bed file 2 that it overlaps (only
    valid when using bedtools).
    exclude: if True, report intervals that DON'T overlap
    sort: sort bed files before taking the intersection
    force_strand: check that the feature and the bed interval are on the same strand (only valid with bedtools)
    force_opposite_strand: if True, check that the feature and the interval are on OPPOSITE strands
    no_name_check: if set to False, checks whether the chromosome names are the same in the too bed files (only valid with bedtools)
    no_dups: if True, only returns each interval once. If set to false, intervals in bed file 1 that overlap several intervals in
    bed file 2 will be returned several times (as many times as there are overlaps with different elements in bed file 2)
    chrom: limit search to a specific chromosome (only valid with bedops, can help in terms of efficiency)
    intersect: rather than returning the entire interval, only return the part of the interval that overlaps an interval in bed file 2.
    hit_count: for each element in bed file 1, return the number of elements it overlaps in bed file 2 (only valid with bedtools)
    intersect_bam: intersect a bam file with a bed file. Requires bam file to be called first
    write_zero: like write_both but also write A intervals that don't overlap with any B intervals,
    write_bed: when intersecting a bam file, write output as bed.'''
    if force_strand and force_opposite_strand:
        raise Exception("force_strand and force_opposite_strand can't both be True")
    hk.make_dir("temp_data/")
    temp_file_name = "temp_data/temp_bed_file{0}.bed".format(random.random())
    #have it write the output to a temporary file
    if use_bedops:
        bedtools_output = run_bedops(bed_file1, bed_file2, force_strand, force_opposite_strand, write_both, chrom, overlap, sort, output_file = temp_file_name, intersect = intersect, hit_number = hit_count, no_dups = no_dups, intersect_bam = intersect_bam, overlap_rec = overlap_rec)
    else:
        bedtools_output = run_bedtools(bed_file1, bed_file2, force_strand, force_opposite_strand, write_both, chrom, overlap, sort, no_name_check, no_dups, output_file = temp_file_name, intersect = intersect, hit_number = hit_count, bed_path = bed_path, intersect_bam = intersect_bam, write_zero = write_zero, overlap_rec = overlap_rec, write_bed = write_bed, exclude = exclude)
    #move it to a permanent location only if you want to keep it
    if output_file:
        hk.run_process(["mv", temp_file_name, output_file])
    else:
        bedtools_output = rw.read_many_fields(temp_file_name, "\t")
    hk.remove_file(temp_file_name)
    return(bedtools_output)

def merge_bed(in_bed, out_bed, distance):
    """
    Strand-specifically merge a BED files.
    Sorts the input BED file first.
    :param in_bed: input file name
    :param out_bed: output file name
    :param distance: maximum distance between
    two elements that are to be merged
    :return: None
    """
    sorted = hk.run_process(["sortBed", "-i", in_bed])
    # the -c and the -o are so that the strand would end up in the
    # right column
    hk.run_process(["mergeBed", "-s", "-d", distance, "-c", "4,5,6", "-o", "distinct,distinct,distinct"], input_to_pipe = sorted, file_for_output=out_bed)

def parse_3ss(junction_file):
    '''
    Parse a BED file of 3' splice sites into a dictionary.
    :param junction_file: BED file of alternating rows of 3' splice sites,
    with distance to the 5' splice site in column 7.
    :return: dictionary with gene names and exon numbers as keys,
    and 3' ss positions with the distance as values.
    '''
    junction_dict = {}
    with open(junction_file) as file:
        for pos, line in enumerate(file):
            line = line.split("\t")
            if len(line) > 0:
                distance = int(line[4])
                junction_dict[line[3]] = [int(line[1]), distance]
    return(junction_dict)

def parse_exon_junctions(junction_file):
    '''
    Parse a BED file of exon junction coordinates into a dictionary.
    :param junction_file: BED file of alternating rows of 5' and 3' splice sites.
    :return: dictionary with 5' splice sites as keys and 3' splice sites as values.
    '''
    junction_dict = {}
    with open(junction_file) as file:
        curr_lines = []
        for pos, line in enumerate(file):
            line = line.split("\t")
            if (pos > 0) and (pos%2 != 0):
                curr_lines.append(line)
                start_line = [i for i in curr_lines if i[3][-1] == "5"][0]
                end_line = [i for i in curr_lines if i[3][-1] == "3"][0]
                name = start_line[3]
                if int(start_line[1]) > int(end_line[1]):
                    junction_dict["{0}_{1}_{2}".format(name, start_line[0], start_line[1])] = int(end_line[1])
                else:
                    junction_dict["{0}_{1}_{2}".format(name, start_line[0], start_line[1])] = int(end_line[1])
                curr_lines = []
            else:
                curr_lines.append(line)
    return(junction_dict)

def peak_pos_in_exon(exon_starts_file, peaks_file, from_end = False, reads_file = False, reads_mode = False):
    """
    Given a set of exons and a set of peaks, make a dictionary with the peaks overlapping each exon.
    :param exon_starts_file: BED file with the starting regions of those exons that have been chosen for study
    :param peaks_file: BED file of peaks
    :param from_end: if True, distances will be calculated from the ends of exons rather than the starts
    :param reads_mode: if True, assume that output is reads and not peaks. The difference
    is that peaks are flat: every position that overlaps with a peak is a 1. Whereas
    with reads, if a position overlaps with more than one read, then you should count
    it as more than one.
    :return: dictionary with junction IDs as keys and a list of all the
    positions that overlap a peak (relative to the start of the exon,
    counting in the direction of transcription). Also return a second dictionary with the centres of the peaks
    (median rounded down to nearest integer).
    """
    # intersect the exon starts and the peaks
    intersect_output = "{0}_{1}".format(exon_starts_file[:-4], peaks_file.split("/")[-1])
    intersect_bed(exon_starts_file, peaks_file, write_both=True, output_file=intersect_output, force_strand=True,
                  no_dups=False, write_zero=True)
    plus = "+"
    if from_end:
        plus = "-"
    out_dict = {}
    with open(intersect_output) as file:
        for line in file:
            line = line.split("\t")
            junction = line[3]
            hk.add_key(junction, [], out_dict)
            # if this exon overlaps with peaks
            if line[6] != ".":
                out_dict[junction].append([])
                peak_start = int(line[7])
                peak_end = int(line[8])
                if line[5] == plus:
                    exon_start = int(line[1])
                    # loop over the nucleotides in the peak
                    for nt in range(peak_start, peak_end):
                        out_dict[junction][-1].append(nt - exon_start)
                else:
                    exon_end = int(line[2])
                    for nt in range(peak_start, peak_end):
                        out_dict[junction][-1].append(exon_end - nt - 1)
    # calculate the centres of the peaks
    out_dict_centres = {i: [int(np.median(j)) for j in out_dict[i]] for i in out_dict}
    # break up separate peaks and just record all the positions once
    if reads_mode:
        out_dict = {i: sorted(list(hk.flatten(out_dict[i]))) for i in out_dict}
    else:
        out_dict = {i: sorted(list(set(hk.flatten(out_dict[i])))) for i in out_dict}
    hk.remove_file(intersect_output)
    return(out_dict, out_dict_centres)

def run_bedops(A_file, B_file, force_strand = False, force_opposite_strand = False, write_both = False, chrom = None, overlap = None, sort = False, output_file = None, intersect = False, hit_number = None, no_dups = False, overlap_rec = None, intersect_bam = None):
    '''
    See intersect_bed for details.
    '''
    if intersect:
        command = "--intersect"
    else:
        command = "--element-of"
    if sort:
        sort_bed(A_file, A_file)
        sort_bed(B_file, B_file)
    bedops_args = ["bedops", "--chrom", "foo", command, "1", A_file, B_file]
    if overlap:
        bedops_args[4] = overlap
    if chrom:
        bedops_args[2] = chrom
        if intersect:
            del bedops_args[4]
    else:
        del bedops_args[1:3]
        if intersect:
            del bedops_args[2]
    if force_strand:
        print("Bedops can't search by strand! Either use bedtools or separate input data by strand!")
        raise Exception
    if force_opposite_strand:
        print("Bedops can't search by strand! Either use bedtools or separate input data by strand!")
        raise Exception
    if write_both:
        print("Bedops can't write both features!")
        raise Exception
    if hit_number:
        print("Bedops hasn't been set up to count the number of overlapping elements. Use bedtools!")
        raise Exception
    if no_dups:
        print("Bedops doesn't print duplicates by default!")
    if overlap_rec:
        print("Bedops hasn't been set up to filter by overlap in second file!")
    if intersect_bam:
        print("Use bedtools to intersect bam and bed!")
        raise Exception
    bedops_output = hk.run_process(bedops_args, file_for_output = output_file)
    return(bedops_output)

def run_bedtools(A_file, B_file, force_strand = False, force_opposite_strand = False, write_both = False, chrom = None, overlap = None, sort = False, no_name_check = False, no_dups = True, hit_number = False, output_file = None, intersect = False, bed_path = None, overlap_rec = None, intersect_bam = None, write_zero = None, write_bed = False, exclude = False):
    '''
    See intersect_bed for details.
    '''
    if write_zero:
        write_option = "-wao"
    elif hit_number:
        write_option = "-c"
    elif write_both:
        write_option = "-wo"
    else:
        write_option = "-wa"
    if sort:
        sort_bed(A_file, A_file)
        sort_bed(B_file, B_file)
    bedtools_args = ["bedtools", "intersect", "-a", A_file,"-b", B_file, write_option]
    if intersect:
        del bedtools_args[-1]
    if overlap:
        bedtools_args.extend(["-f", str(overlap)])
    if overlap_rec:
        bedtools_args.append("-r")
    if force_strand:
        bedtools_args.append("-s")
    elif force_opposite_strand:
        bedtools_args.append("-S")
    if no_name_check:
        bedtools_args.append("-nonamecheck")
    if no_dups:
        bedtools_args.append("-u")
    if chrom:
        print("Bedtools cannot be restricted to a single chromosome. Use bedops!")
        raise Exception
    if hit_number and no_dups:
        print("When counting hits, each interval in the first bed file is only reported once by default. Set no_dups to False!")
        raise(Exception)
    if bed_path:
        bedtools_args[0] = "{0}{1}".format(bed_path, bedtools_args[0])
    if exclude:
        bedtools_args.append("-v")
    if intersect_bam:
        if A_file[-4:] != ".bam":
            print("Bam file must be called first")
            raise Exception
        if B_file[-4:] == ".bam":
            print("BAM file must be called first")
            raise Exception
        bedtools_args = ["intersectBed", write_option, "-abam", A_file, "-b", B_file]
        if write_bed:
            bedtools_args.append("-bed")
    try:
        bedtools_output = hk.run_process(bedtools_args, file_for_output = output_file)
    except FileNotFoundError:
        bedtools_args[0] = "intersectBed"
        bedtools_output = hk.run_process(bedtools_args, file_for_output = output_file)
    return(bedtools_output)

def snr_bed(in_bed, out_bed, five_prime_most=False, shift_plus=0, shift_minus=0):
    """
    Given a BED file, make a new BED file that only has the 3' most nucleotide
    of each read.
    :param in_bed: input BED file
    :param out_bed: output BED file
    :param five_prime_most: if True, the 5'-most nucleotide is written instead.
    :param shift_plus: by how much to shift the positions on the plus strand
    (minus strand if five_prime_most==True)
    :param shift_minus: by how much to shift the positions on the minus strand
    (plus strand if five_prime_most==True)
    :return: None
    """
    plus = "+"
    if five_prime_most:
        plus = "-"
    with open(in_bed) as infile, open(out_bed, "w") as outfile:
        reader = csv.reader(infile, delimiter="\t")
        writer = csv.writer(outfile, delimiter="\t")
        for line in reader:
            outline = line.copy()
            if line[5] == plus:
                outline[1] = str(int(line[2]) - 1 + shift_plus)
                outline[2] = str(int(line[2]) + shift_plus)
            else:
                outline[2] = str(int(line[1]) + 1 + shift_minus)
                outline[1] = str(int(line[1]) + shift_minus)
            outline = outline[:6]
            writer.writerow(outline)

def sort_bed(input_file_name, output_file_name):
    '''
    Sort a bed file.
    '''
    #This is done via a temp file because that way you can specify the same file as input and output file and thus
    #overwrite the unsorted file with the sorted one.
    temp_file_name = "temp_data/temp_sorted_bed{0}.bed".format(random.random())
    hk.run_process(["sort-bed", input_file_name], file_for_output = temp_file_name)
    hk.run_process(["mv", temp_file_name, output_file_name])
    hk.remove_file(temp_file_name)

def trim_sequence(sequence, phase):
    """
    Trim a sequence so that it would start and end with a full codon.
    :param sequence: sequence string
    :param phase: integer specifying the phase of the sequence. Note that the phase should be entered
    like in a GTF file and is then converted internally into an in-house system.
    Phase 0 means that the sequence starts with the first base of a codon.
    Phase 1: the sequence starts with the third base of a codon (there is one extra base).
    Phase 2: the sequence starts with the second base of a codon (there are two extra bases).
    :return: trimmed sequence string
    """
    conversion_dict = {0:0, 1:2, 2:1}
    phase = conversion_dict[phase]
    #trim the start based on the phase information that was given
    if phase == 0:
        pass
    elif phase == 1:
        sequence = sequence[2:]
    elif phase == 2:
        sequence = sequence[1:]
    else:
        raise Exception("Invalid phase information!")
    #trim the end based on the length of the remaining sequence.
    if len(sequence)%3 == 0:
        return(sequence)
    elif len(sequence)%3 == 1:
        return(sequence[:-1])
    else:
        return(sequence[:-2])

def write_intron_lariat_pos_from_exons(exons, outbed, add_chr = False):
    """
    Write the last nucleotide of each intron into a BED file.
    :param exons: exon dictionary
    :param outbed: output file name
    :param add_chr: if True, "chr" will be prefixed to chromosome names.
    :return: None
    """
    introns = get_introns(exons)
    with open(outbed, "w") as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        for trans in introns:
            for pos, intron in enumerate(introns[trans]):
                strand = intron[6]
                if strand == "+":
                    start = intron[4] - 1
                    end = intron[4]
                else:
                    start = intron[3] - 1
                    end = intron[3]
                chrom = intron[0]
                if add_chr:
                    if chrom[0] not in ["G", "K"]:
                        chrom = "chr{0}".format(chrom)
                writer.writerow([chrom, start, end, "{0}.{1}".format(trans, pos), ".", strand])

def write_intron_start_peak_from_exons(exons, outbed, add_chr = False, alt_start = None, alt_end = None):
    """
    Write nts 20-40 of each intron into a BED file.
    :param exons: exon dictionary
    :param outbed: output file name
    :param add_chr: if True, "chr" will be prefixed to chromosome names.
    :param alt_start: if specified, start region from this nucleotide instead.
    :param alt_end: if specified, end region at this nucleotide instead.
    :return: None
    """
    begin = 20
    finish = 40
    if alt_start:
        begin = alt_start
    if alt_end:
        finish = alt_end
    introns = get_introns(exons)
    with open(outbed, "w") as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        for trans in introns:
            for pos, intron in enumerate(introns[trans]):
                if (intron[4] - intron[3] + 1) > finish:
                    strand = intron[6]
                    if strand == "+":
                        start = intron[3] + (begin - 1)
                        end = intron[3] + (finish - 1)
                    else:
                        start = intron[4] - finish
                        end = intron[4] - begin
                    chrom = intron[0]
                    if add_chr:
                        if chrom[0] not in ["G", "K"]:
                            chrom = "chr{0}".format(chrom)
                    writer.writerow([chrom, start, end, "{0}.{1}".format(trans, pos), ".", strand])

def write_si_pos_from_exons(exons, outbed, add_chr = False):
    """
    Write the last nucleotide of each exon into a BED file.
    :param exons: exon dictionary
    :param outbed: output file name
    :param add_chr: if True, "chr" will be prefixed to chromosome names.
    :return: None
    """
    with open(outbed, "w") as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        for trans in exons:
            for pos, exon in enumerate(exons[trans]):
                strand = exon[6]
                if strand == "+":
                    start = exon[4] - 1
                    end = exon[4]
                else:
                    start = exon[3] - 1
                    end = exon[3]
                chrom = exon[0]
                if add_chr:
                    if chrom[0] not in ["G", "K"]:
                        chrom = "chr{0}".format(chrom)
                writer.writerow([chrom, start, end, "{0}.{1}".format(trans, pos), ".", strand])
