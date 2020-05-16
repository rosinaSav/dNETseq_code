'''
Given a BED file of NET-seq (or other) reads, filter out reads that end at the
last nucleotide of an exon (putative splicing intermediates) or the last
nucleotide of an intron (putative intron lariats).
'''
import coord_ops as co
import csv
import housekeeping as hk
import read_and_write as rw

def check_3prime_match(infile, outfile):
    """
    Given a BEDtools intersect file between two BED files of genomic coordinates
    (where the two intersecting regions have been written out),
    filter out overlaps where the 3' ends of the two regions don't match
    exactly. Also convert chrN chromosome labels to just N.
    :param infile: input BEDtools intersect file
    :param outfile: the file with records where the 3' ends don't match
    exactly removed
    :return: None
    """
    with open(infile) as ifile, open(outfile, "w") as ofile:
        reader = csv.reader(ifile, delimiter = "\t")
        writer = csv.writer(ofile, delimiter = "\t")
        # the intersect file has all the reads that overlap the relevant position
        # so you need to pick out the ones where the 3' end matches exactly
        for line in reader:
            if line[5] == "+":
                if line[2] == line[19]:
                    line[0] = line[0].lstrip("chr")
                    writer.writerow(line)
            else:
                if line[1] == line[18]:
                    line[0] = line[0].lstrip("chr")
                    writer.writerow(line)

def main():
    description = "Given a BED file of reads, filter out reads whose " \
                  "3' end maps to the last nucleotide of an intron or" \
                  "the last nucleotide of an exon."
    args = hk.parse_arguments(description, ["reads_file", "gtf", "outfile"])
    reads_file, gtf, outfile = args.reads_file, args.gtf, args.outfile

    print("Getting intron lariat positions...")

    # read in exon coordinates
    exons = rw.read_gtf(gtf, element="exon", gene=False)
    # make a BED file with the last positions of introns
    intron_lariat_bed = "{0}_intron_lariat_pos_all_exons.bed".format(reads_file[:-4])
    co.write_intron_lariat_pos_from_exons(exons, intron_lariat_bed, add_chr = True)

    # intersect the reads with intron lariat positions
    intron_lariat_intersect_file_name = "{0}_intersect_with_intron_lariat_pos_all_exons.bed".format(reads_file[:-4])
    co.intersect_bed(reads_file, intron_lariat_bed, force_strand=True, write_both=True, no_dups=False, output_file=intron_lariat_intersect_file_name)
    hk.remove_file(intron_lariat_bed)
    intron_lariat_reads_file = "{0}_intron_lariat_reads_all_exons.bed".format(reads_file[:-4])
    # check that the reads end exactly at intron lariat positions
    check_3prime_match(intron_lariat_intersect_file_name, intron_lariat_reads_file)
    hk.remove_file(intron_lariat_intersect_file_name)

    # write BED with the last positions of exons
    splice_intermediate_bed = "{0}_splice_intermediate_pos_all_exons.bed".format(reads_file[:-4])
    co.write_si_pos_from_exons(exons, splice_intermediate_bed, add_chr = True)

    print("Getting splice intermediate positions.")

    # intersect the reads with splice intermediate positions
    splice_intermediate_intersect_file_name = "{0}_intersect_with_SI_pos_all_exons.bed".format(reads_file[:-4])
    co.intersect_bed(reads_file, splice_intermediate_bed, force_strand=True, write_both=True, no_dups=False, output_file=splice_intermediate_intersect_file_name)
    hk.remove_file(splice_intermediate_bed)
    SI_reads_file = "{0}_SI_reads_all_exons.bed".format(reads_file[:-4])
    # check that the reads end exactly at the end of the exon
    check_3prime_match(splice_intermediate_intersect_file_name, SI_reads_file)
    hk.remove_file(splice_intermediate_intersect_file_name)

    print("Concatenating the two files.")

    # concatenate the IL and SI read files so you could exclude both in one go
    combined_file = "{0}_SI_and_IL_reads_all_exons.bed".format(reads_file[:-4])
    hk.run_process(["cat", SI_reads_file, intron_lariat_reads_file], file_for_output=combined_file)

    hk.remove_file(SI_reads_file)
    hk.remove_file(intron_lariat_reads_file)

    # do an exclusive intersect, requiring 1.0 overlap for both A and B, to remove the
    # putative intron lariat reads from the main reads file
    co.intersect_bed(reads_file, combined_file, overlap=1, overlap_rec=1, force_strand=True, no_dups=False, exclude=True, output_file=outfile)

    hk.remove_file(combined_file)

if __name__ == "__main__":
    main()