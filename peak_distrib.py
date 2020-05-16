'''
Build meta-profile of peak/read densities.

ARGUMENTS
peaks_file: BED file with peak coordinates
GTF: Ensembl GTF file for the species in question
exon_starts_file: BED file with the coordinates of the regions
for which the meta-profile is desired. The name field needs to be formatted as XXXXX.N where
XXXXX is the Ensembl transcript ID and N is the exon number (counting
from 0, coding exons only).
output_file: name of the output_file (matrix with regions in rows, nucleotides
in columns, and a 1 in cells that overlap with a peak)
reads_file: BED file with the corresponding NET-seq reads (this is because the script also
outputs a file with the coverage of each region, calculated using bedtools coverage)
from_end: if True, the meta-profile will be reversed, so that the first position
in the output matrix corresponds to the last position of each region
intronic: if True, the region is from the intron upstream of the exon rather than the
exon itself
limit: how long of a region to analyze, starting from the start of the BED regions
nts_before_start: if the region contains both an exonic or an intronic part,
this variable holds the length of the intronic region. This is important to allow
filtering by exon size.
noncoding: if True, it means that the analysis will consider all exons, not just coding ones
reads_mode: if True, it means that the input BED file holds reads rather than peaks. Important
because with reads, only the 3'most position is analysed.
'''

import coord_ops as co
import housekeeping as hk
import read_and_write as rw
from splice_distance import write_dist_mat

def main():

    description = "Record the distribution of peaks for different exons."
    args = hk.parse_arguments(description, ["peaks_file", "gtf", "exon_starts_file", "output_file", "reads_file", "from_end", "intronic", "limit", "nts_before_start", "noncoding", "reads_mode"], flags = [5, 6, 9, 10], ints = [7, 8])
    peaks_file, gtf, exon_starts_file, output_file, reads_file, from_end, intronic, limit, nts_before_start, noncoding, reads_mode = args.peaks_file, args.gtf, args.exon_starts_file, args.output_file, args.reads_file, args.from_end, args.intronic, args.limit, args.nts_before_start, args.noncoding, args.reads_mode

    if noncoding:
        exons = rw.read_gtf(gtf, "exon", gene=False)
    else:
        exons = rw.read_gtf(gtf, "CDS", gene=False)

    # the 3' ss that will be analyzed
    valid_junctions = rw.read_many_fields(exon_starts_file, "\t")
    # pull out the column with transcript IDs
    valid_junctions = [i[3] for i in valid_junctions]

    lengths_dict = co.get_lengths(exons, valid_junctions, intronic=intronic)
    if nts_before_start:
        lengths_dict = {i: lengths_dict[i] + nts_before_start for i in lengths_dict}

    coverage_file_name = "{0}_{1}_coverage.bed".format(exon_starts_file[:-4], reads_file.split("/")[-1][:-4])
    co.get_coverage(exon_starts_file, reads_file, coverage_file_name)

    peak_distances_all, peak_centres = co.peak_pos_in_exon(exon_starts_file, peaks_file, from_end = from_end, reads_mode = reads_mode)

    write_dist_mat(peak_distances_all, limit, output_file, lengths_dict, "{0}_intron_names.txt".format(output_file[:-4]), None)

    write_dist_mat(peak_centres, limit, "{0}_centres.txt".format(output_file[:-4]), lengths_dict,
                                                               "{0}_centres_intron_names.txt".format(output_file[:-4]), None)

if __name__ == "__main__":
    main()