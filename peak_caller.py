'''Peak caller for detecting regions where NET-seq read density is significantly
higher than expected by chance for a given transcript. For further documentation,
run python3 peak_caller.py --help. Requires BEDtools (built with v2.29.0) and
numpy (built with v1.17.2).'''

import coord_ops as co
import csv
import housekeeping as hk
import numpy as np
import read_and_write as rw
import scipy.signal

def filter_peaks(in_peak_bed, read_bed, read_count_file, out_peak_bed, min_reads_per_peak, min_peak_length, stats_file, no_PCR_filter = False):
    """
    Filter a BED file of peaks by removing peaks that overlap fewer than
    _min_reads_per_peak_ reads or are shorter than _min_peak_length_.
    Also write some stats in _stats_file_.
    :param in_peak_bed: BED file with peak coordinates
    :param read_bed: BED file with read coordinates
    :param read_count_file: text file with one row per transcript, containing all the significant positions
    along with their read counts
    :param out_peak_bed: file for filtered BED
    :param min_reads_per_peak: minimum number of reads per peak
    :param min_peak_length: minimum length of peak
    :param stats_file: file for the output stats
    :param no_PCR_filter: if True, no filtering of potential PCR duplicates will be performed
    :return: None
    """
    PCR_threshold = 0.9
    if no_PCR_filter:
        PCR_threshold = 1
    # parse read counts per significant position
    read_counts = {}
    with open(read_count_file) as file:
        for line in file:
            line = line.rstrip("\n").split("\t")
            curr_counts = [i.split(":") for i in line[1:]]
            curr_counts = {int(i[0]): int(i[1]) for i in curr_counts}
            read_counts[line[0]] = curr_counts
    # intersect the peaks with the original reads to count the
    # number of reads per peak
    intersect_file = "{0}_{1}_intersect.bed".format(in_peak_bed[:-4], read_bed.split("/")[-1][:-4])
    co.intersect_bed(in_peak_bed, read_bed, force_strand=True,
                     write_both=True, no_dups=False, output_file=intersect_file,
                     hit_count=True)
    lengths = []
    counts = []
    with open(intersect_file) as in_file, open(out_peak_bed, "w") as out_file, open(stats_file, "w") as stats_file:
        stats_file.write("transcript\tlength\tread_count\n")
        for line in in_file:
            line = line.rstrip("\n").split("\t")
            read_count = int(line[6])
            if read_count >= min_reads_per_peak:
                end = int(line[2])
                start = int(line[1])
                length = end - start
                if length >= min_peak_length:
                    # check that the peak hasn't been merged between more than one transcript
                    trans = line[3]
                    if "," not in trans:
                        # check that the most enriched position doesn't account for more than 90%
                        # of the density
                        curr_read_counts = [read_counts[trans][i] if i in read_counts[trans] else 0 for i in range(start, end)]
                        maximum = np.max(curr_read_counts)
                        all_pos = np.sum(curr_read_counts)
                        if maximum/all_pos <= PCR_threshold:
                            # so the file could be visualized in a genome browser
                            line[4] = "100"
                            out_file.write("{0}\n".format("\t".join(line)))
                            stats_file.write("{0}\t{1}\t{2}\n".format(trans, length, read_count))
                            lengths.append(length)
                            counts.append(read_count)
        print("Found a total of {0} peaks.".format(len(lengths)))
        print("The median peak length is {0}.".format(np.median(lengths)))
        print("The median peak read count is {0}.".format(np.median(counts)))

def generate_random_counts(read_counts, iterations):
    """
    Given a numpy array of read counts at each position in a transcript,
    generate _iterations_ random arrays where the reads have been shuffled
    and concatenate them all into one big array.
    :param read_counts: numpy array where each element corresponds to one position
    in a transcript and contains the read count at that position
    :param iterations: number of times that you should shuffle the reads
    :return: numpy array
    """
    pos_number = read_counts.shape[0]
    out_array = np.zeros(pos_number * iterations)
    read_number = int(np.sum(read_counts))
    for iter in range(iterations):
        curr_read_positions = np.random.randint(0, high=pos_number, size=read_number, dtype="int")
        curr_counts = [np.count_nonzero(curr_read_positions == i) for i in range(pos_number)]
        out_array[(iter * pos_number): ((iter + 1) * pos_number)] = curr_counts
    return(out_array)

def generate_random_windowed_counts(read_counts, iterations, window_size, no_slide, exclude = False):
    """
    Given a numpy array of read counts at each position in a transcript,
    generate _iterations_ random arrays where the reads have been shuffled,
    with read counts averaged over a sliding window of size window_size,
    and concatenate them all into one big array.
    :param read_counts: numpy array where each element corresponds to one position
    in a transcript and contains the read count at that position
    :param iterations: number of times that you should shuffle the reads
    :param window_size: size of sliding window to use
    :param no_slide: if True, use sequential windows of window_size bp rather than a sliding window
    :param exclude: optionally accepts a tuple (relative start and end coordinate). The positions
    between these coordinates will be ignored in the randomization.
    :return: numpy array
    """
    if exclude:
        read_counts = read_counts[list(range(exclude[0])) + list(range(-(read_counts.shape[0] - exclude[1]), 0))]
    pos_number = read_counts.shape[0]
    out_array = np.zeros(pos_number * iterations)
    read_number = int(np.sum(read_counts))
    for iter in range(iterations):
        curr_read_positions = np.random.randint(0, high=pos_number, size=read_number, dtype="int")
        curr_counts = np.bincount(curr_read_positions, minlength=pos_number)
        curr_counts = window_read_count_array(curr_counts, window_size, no_slide)
        out_array[(iter * pos_number): ((iter + 1) * pos_number)] = curr_counts
    return(out_array)

def get_reads_per_pos(reads_file, transcript_bed):
    """
    Given a BED file of reads and a BED file of transcript coordinates,
    make a dictionary with transcript IDs as keys and number of reads per position,
    as well as the absolute coordinates of the nucleotides, as values.
    :param reads_file: BED file with read coordinates
    :param transcript_bed: BED file with transcript coordinates
    :return: dictionary with numbers of reads per position
    """
    # intersect the transcripts and the reads, so you'd have an output file where
    # the transcript coordinates are followed by the overlapping read
    intermediate_file = "{0}_{1}read_per_pos_intermediate.bed".format(reads_file[:-4], transcript_bed.split("/")[-1][:-4])
    co.intersect_bed(transcript_bed, reads_file, force_strand=True,
                   write_both=True, no_dups = False, write_zero = False, output_file=intermediate_file)
    reads_per_pos = {}
    total = hk.line_count(intermediate_file)
    print("Calculating the number of reads per position in each transcript...")
    with open(intermediate_file, newline="") as file:
        file_reader = csv.reader(file, delimiter="\t")
        for pos, line in enumerate(file_reader):
            if pos%100000 == 0:
                print("{0}/{1}".format(pos, total))
            # prefix the chromosome and the strand to the transcript name cause you'll
            # need it later
            trans_name = line[3]
            trans_name = "{0}.{1}.{2}".format(line[0], line[5], trans_name)
            reads_per_pos = hk.add_key(trans_name, {"reads": {}}, reads_per_pos)
            strand = line[5]
            if strand == "+":
                position = int(line[8]) - 1
            else:
                position = int(line[7])
            reads_per_pos[trans_name]["reads"] = hk.add_key(position, 0, reads_per_pos[trans_name]["reads"])
            reads_per_pos[trans_name]["reads"][position] = reads_per_pos[trans_name]["reads"][position] + 1
            reads_per_pos[trans_name] = hk.add_key("coords", (int(line[1]), int(line[2])), reads_per_pos[trans_name])
    hk.remove_file(intermediate_file)
    return reads_per_pos

def get_sig_pos_trans(read_counts, exons, window_size, iterations, significance_threshold, no_slide, exclude_focal, with_ups_intron):
    """
    Return positions in a transcript that are significantly enriched for reads (using a sliding
    window).
    :param read_counts: an array containing the read counts at each position in a transcript
    :param exons: list of lists of GTF lines of the exon coordinates associated with the transcript
    :param window_size: size of sliding window
    :param iterations: number of iterations to perform when calculating the empirical distribution
    :param significance_threshold: significance threshold
    :param no_slide: if True, sequential windows of window_size bp will be used instead of a sliding window
    :param exclude_focal: if True, the significance threshold will be set separately for each exon
    and intron, excluding that exon or intron from the significance calculations
    :param with_ups_intron: if True, when excluding focal, exclude not just the focal exon but also the upstream
    intron
    :return:
    """
    windowed_counts = window_read_count_array(read_counts, window_size=window_size, no_slide=no_slide)
    if not exclude_focal:
        random_windowed_counts = generate_random_windowed_counts(read_counts, iterations,
                                                                 window_size=window_size, no_slide=no_slide)
        threshold = np.percentile(random_windowed_counts, (1 - significance_threshold) * 100)
        if exons[0][6] == "+":
            coords = (exons[0][3], exons[-1][4])
        else:
            coords = (exons[-1][3], exons[0][4])
        sig_positions = sig_pos(coords, windowed_counts, threshold)
    else:
        # if it's a single-exon gene, don't look for peaks
        if len(exons) == 1:
            return([])
        intervals = co.get_exon_and_intron_coords(exons, with_ups_intron=with_ups_intron)
        sig_positions = []
        cumul = 0
        for interval in intervals:
            length = interval[1] - interval[0]
            interval_end = cumul + length
            random_windowed_counts = generate_random_windowed_counts(read_counts, iterations,
                                                                     window_size=window_size, no_slide=no_slide, exclude=(cumul, interval_end))
            threshold = np.percentile(random_windowed_counts, (1 - significance_threshold) * 100)
            sig_positions.extend(sig_pos(interval, windowed_counts[cumul: interval_end], threshold))
            cumul = interval_end
    return(sig_positions)

def make_read_count_array(trans_reads):
    """
    Given a reads per position dictionary, as output by reads_per_pos,
    make a numpy array with the number of reads for each position
    in the transcript.
    :param trans_reads: one element from output dictionary from reads_per_pos
    :return: numpy array
    """
    start = trans_reads["coords"][0]
    end = trans_reads["coords"][1]
    read_counts = np.zeros(end - start, dtype = "int")
    for pos, nt in enumerate(range(start, end)):
        if nt in trans_reads["reads"]:
            read_counts[pos] = trans_reads["reads"][nt]
    return read_counts

def make_reads_dict_from_array(read_counts, coords):
    """
    Given a read count array, as well as the transcript start and end coordinate,
    make a new dictionary of read counts. Used after randomization of read counts in a negative
    control.
    :param read_counts: array where each position has the number of reads at that position
    in the transcript
    :param coords: tuple of the chromosomal start and end cooridnate of the transcript (0-based)
    :return: dictionary of read counts per position
    """
    new_reads_dict = {}
    for pos, nt in enumerate(range(coords[0], coords[1])):
        if read_counts[pos] != 0:
            new_reads_dict[nt] = int(read_counts[pos])
    return(new_reads_dict)

def sig_pos(coords, windowed_counts, threshold):
    '''
    Determine which positions are significantly enriched in reads in a transcript based on whether
    they exceed a previously determined threshold.
    :param coords: 2-tuple with start and end coordinate of transcript
    :param windowed_counts: array of read counts at the positions, avergaed over a window
    :param threshold: threshold that a position would need to exceed to be significant
    :return: list of significant positions
    '''
    sig_positions = [nt for pos, nt in enumerate(range(coords[0], coords[1])) if windowed_counts[pos] > threshold]
    return sig_positions

def window_read_count_array(read_counts, window_size, no_slide):
    """
    Given an array where each element gives the read count at a give position,
    pass a sliding window over each exon and make an array where each position is
    the average of the window centered at that position.
    :param read_counts: array of read count at successive positions
    in a transcript
    :param window_size: Size of the sliding window. Must be an odd number!
    :param no_slide: if True, use sequential windows of window_size bp rather than a sliding window
    :return: numpy array
    """
    if sum(read_counts) == 0:
        return([0 for i in read_counts])
    shift = int(window_size / 2)
    temp_array = np.pad(read_counts, shift, mode="edge")
    if no_slide:
        temp_out = np.zeros(len(temp_array))
        for position in range(0, len(temp_out), window_size):
            curr_mean = np.mean(temp_array[position: position + window_size])
            temp_out[position: position + window_size] = curr_mean
        out_array = temp_out[shift: (len(temp_out) - shift)]
    else:
        out_array = scipy.signal.fftconvolve(temp_array, np.ones((window_size,))/window_size, mode="valid")
    return(out_array)

def write_raw_peaks(reads_per_pos, output_bed, counts_file, exons, significance_threshold = 0.05, iterations = 10, min_read_count = 5, window_size = 10, neg_control = False, no_slide = False, exclude_focal = False, with_ups_intron = False):
    """
    Given the numbers of reads at different positions in different transcripts,
    determine significantly enriched positions based on permuting the read
    positions.
    :param reads_per_pos: A dictionary of read positions in transcripts,
    as output by get_reads_per_pos.
    :param output_bed: File for output.
    :param counts_file: File that will contain the read counts for each significant position.
    :param exons: Dictionary with transcript IDs as keys and exon GTF coordinates as keys.
    :param significance_threshold: significance threshold for considering
     a position as enriched for reads.
    :param iterations: number of random iterations to perform for each
    transcript.
    :param min_read_count: minimum total read count for a transcript to be considered.
    :param window_size: size of the window over which to average read counts
    :param neg_control: if True, the reads within each transcript will be randomly shuffled before analysis
    :param no_slide: if True, use sequential windows of 5 bp rather than a sliding window
    :param exclude_focal: if True, the threshold for significance will be set separately for each exon/intron
    :param with_ups_intron: if True, then each exon will be considered as a unit with its upstream intorn when
    calling peaks (presumes exclude_focal is True).
    :return: None
    """
    new_reads_file = None
    if neg_control:
        new_reads_file = "{0}_new_reads.bed".format(output_bed[30:-4])
        new_polII_bed = open(new_reads_file, "w", newline="")
        polII_writer = csv.writer(new_polII_bed, delimiter="\t")
    print("\nLooking for significant positions...")
    with open(output_bed, "w") as file, open(counts_file, "w") as counts_f:
        peak_counter = 0
        total_trans = len(reads_per_pos)
        for pos, transcript in enumerate(reads_per_pos):
            counts_write = [transcript.split(".")[-1]]
            if pos % 1000 == 0:
                print("{0}/{1} transcripts processed.".format(pos, total_trans))
            current_trans = reads_per_pos[transcript]
            read_counts = make_read_count_array(current_trans)
            # if there are such few reads in the transcript,
            # that any peak would have too few reads to be considered,
            # then don't even analyze it
            if np.sum(read_counts) >= min_read_count:
                transcript = transcript.split(".")
                # for a negative control, randomly shuffle the reads
                if neg_control:
                    # make a new read count array where the read positions have been randomized
                    read_counts = generate_random_counts(read_counts, 1)
                    # also make a new reads dictionary (will be used in the last filtering step
                    # at the end of the function)
                    current_trans["reads"] = make_reads_dict_from_array(read_counts, current_trans["coords"])
                    # also write to a new reads bed file
                    template = [transcript[0], 0, 0, transcript[-1], ".", transcript[1]]
                    for nt in current_trans["reads"]:
                        to_write = template.copy()
                        to_write[1] = nt
                        to_write[2] = nt + 1
                        for occurrence in range(current_trans["reads"][nt]):
                            polII_writer.writerow(to_write)
                sig_positions = get_sig_pos_trans(read_counts, exons[transcript[2]], window_size, iterations, significance_threshold, no_slide, exclude_focal=exclude_focal, with_ups_intron=with_ups_intron)
                template = [transcript[0], None, None, transcript[2], "100", transcript[1]]
                new_sig_pos = []
                for sig_pos in sig_positions:
                    # i.e. if the read count isn't 0. This is to trim empty edges.
                    # if there are any 0 positions in the middle of the peak,
                    # that's fine because they'll still be part of the peak
                    # because of the merge.
                    if sig_pos in current_trans["reads"]:
                        new_sig_pos.append(sig_pos)
                        curr_line = template.copy()
                        curr_line[1] = str(sig_pos)
                        curr_line[2] = str(sig_pos + 1)
                        file.write("{0}\n".format("\t".join(curr_line)))
                        peak_counter = peak_counter + 1
                        read_count = current_trans["reads"][sig_pos]
                        counts_write.append("{0}:{1}".format(sig_pos, read_count))
            counts_write = "\t".join(counts_write)
            counts_f.write("{0}\n".format(counts_write))
    print("Found {0} significant positions!".format(peak_counter))
    if neg_control:
        new_polII_bed.close()
    return(new_reads_file)

def main():

    description = "Call peaks in a BED file of NET-seq reads."
    help_info = ["BED file (at least a BED6) with NET-seq reads. Should be single-nucleotide resolution (each BED region is the 3' end of a read.).",
                 "Ensembl GTF file for the relevant species. Ensure that chromosome names are formatted the same way in both the GTF and the BED file with reads!",
                 "BED file with the coordinates of the transcripts to analyze. Only the name field is read, hence the others can hold placeholders. The name field must contain transcript IDs from the GTF file.",
                 "Name of the output file (BED file with peak coordinates).",
                 "Alpha value for calling a position as having a significantly higher local read denisty than expected by chance. Default: 0.01.",
                 "Merge distance: adjacent peaks will be merged if they are closer than this many nucleotides. Default: 21.",
                 "Minimum reads per peak. Default: 10.",
                 "The number of times the read position randomization should be performed for each transcript. Higher values make the significance calculation (marginally) more robust, however, they also make the programme very slow. Default: 5.",
                 "Minimum length of a peak in nucleotides. Default: 5.",
                 "Size of the sliding window to use when calculating the local read density. It may be sensible to set this to the same value as the merge distance. Should be an odd integer. Default: 21",
                 "The analysis will be performed this many times, with the output files numbered. Useful for running many negative control simulations at once. Default: 1.",
                 "Read positions will be shuffled within each transcript before analysis. This should disrupt any signal and should give a flat peak density profile.",
                 "Instead of a sliding window, adjacent non-overlapping windows will be used when calculating the local read density.",
                 "When calling peaks in a given exon/intron, do not include that exon/intron in the read position randomization.",
                 "When --exclude_focal is set, count an exon and its upstream intron as a single unit (except for the first exon).",
                 "Don't filter out likely PCR duplicates (peaks where more than 90%% of the reads come from a single nucleotide position).)"]
    defaults = {4: 0.01, 5: 21, 6: 10, 7: 5, 8: 5, 9: 21, 10: 1}
    args = hk.parse_arguments(description, ["reads_file", "gtf", "trans_active_file", "output_file", "significance_threshold", "merge", "min_reads_per_peak", "iterations", "min_peak_length", "window_size", "runs", "neg_control", "no_slide", "exclude_focal", "with_ups_intron", "no_PCR_filter"], floats = [4], ints = [5, 6, 7, 8, 9, 10], flags = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15], detailed_help = help_info, defaults=defaults)
    reads_file, gtf, trans_active_file, output_file, significance_threshold, merge, min_reads_per_peak, iterations, min_peak_length, window_size, runs, neg_control, no_slide, exclude_focal, with_ups_intron, no_PCR_filter = args.reads_file, args.gtf, args.trans_active_file, args.output_file, args.significance_threshold, args.merge, args.min_reads_per_peak, args.iterations, args.min_peak_length, args.window_size, args.runs, args.neg_control, args.no_slide, args.exclude_focal, args.with_ups_intron, args.no_PCR_filter

    print("Merge distance: {0}".format(merge))
    print("Minimum number of reads per peak: {0}".format(min_reads_per_peak))
    print("Minimum peak length: {0}".format(min_peak_length))
    print("Window size: {0}".format(window_size))
    print("Significance level: {0}".format(significance_threshold))
    print("Randomization iterations to perform: {0}".format(iterations))
    print("Runs: {0}".format(runs))

    neg_str = ""
    if neg_control:
        neg_str = "_neg_control"

    slide_str = ""
    if no_slide:
        slide_str = "_no_slide"
    intron_str = ""
    if with_ups_intron:
        intron_str = "w_ups_intr"

    # 0. make a BED file with the coordinates of transcripts

    transcripts_file = "{0}_transcripts.bed".format(gtf[:-4])
    co.get_transcripts(gtf, transcripts_file, add_chr = True)
    exons = rw.read_gtf(gtf, "exon")

    # 1. intersect the two files, loop over the result and make a
    # dictionary of reads per pos for each transcript, which has reads

    reads_per_pos = get_reads_per_pos(reads_file, transcripts_file)
    # only leave transcriptionally active genes (one isoform per gene)
    trans_active_genes = rw.read_many_fields(trans_active_file, "\t")[1:]
    # pull out the column with transcript IDs
    trans_active_genes = [i[3] for i in trans_active_genes]
    reads_per_pos = {i: reads_per_pos[i] for i in reads_per_pos if i.split(".")[-1] in trans_active_genes}

    for sim in range(runs):

        print("**********{0}**********".format(sim))

        # 2. for each transcript, randomly reshuffle the reads and calculate the
        # nth percentile depending on what the significance threshold is
        # keep positions that are higher than that threshold and write to BED file

        raw_peak_bed = "{0}_{1}_raw_peaks{2}_{3}_{4}{5}{6}{7}_{8}_sim.bed".format(reads_file[:-4], gtf.split("/")[-1][:-4], iterations, min_reads_per_peak, window_size, neg_str, intron_str, slide_str, sim)
        read_count_file = "{0}_{1}_read_counts{2}_{3}{4}{5}_{6}_sim.txt".format(reads_file[:-4], gtf.split("/")[-1][:-4], iterations, window_size, neg_str, intron_str, sim)
        new_reads_file = write_raw_peaks(reads_per_pos, raw_peak_bed, read_count_file, exons, iterations = iterations, min_read_count=min_reads_per_peak, window_size=window_size, neg_control=neg_control, no_slide=no_slide, exclude_focal=exclude_focal, with_ups_intron=with_ups_intron)
        if neg_control:
            reads_file = new_reads_file

        # 3. merge peaks

        merged_peak_bed = "{0}_{1}_merged_peaks{2}_{3}_{4}{5}{6}{7}_{8}_sim.bed".format(reads_file[:-4], gtf.split("/")[-1][:-4], iterations, window_size, merge, neg_str, slide_str, intron_str, sim)
        co.merge_bed(raw_peak_bed, merged_peak_bed, merge)
        print("Before filtering, there are {0} peaks.".format(hk.line_count(merged_peak_bed)))

        # 4. filter out peaks that don't have enough reads or are too short.
        # Write final results to file and also write a stats file with the size,
        # read count and overlapping transcript of the peaks

        stats_file = "{0}_stats_{1}_sim.txt".format(output_file[:-4], sim)
        filter_peaks(merged_peak_bed, reads_file, read_count_file, "{0}_{1}_sim.bed".format(output_file[:-4], sim), min_reads_per_peak, min_peak_length, stats_file, no_PCR_filter=no_PCR_filter)

if __name__ == "__main__":
    main()