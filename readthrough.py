import coord_ops as co
import csv
import housekeeping as hk
import read_and_write as rw
from splice_distance import write_dist_mat

def main():
    description = "Prepare a BED file with the TES coordinates of transcriptionally" \
                  "active genes and make a metagene of reads within this region."

    args = hk.parse_arguments(description, ["trans_act_file", "gtf", "start_coord", "end_coord", "outname", "reads_file"], ints = [2, 3])
    trans_act_file, gtf, start_coord, end_coord, outname, reads_file = args.trans_act_file, args.gtf, args.start_coord, args.end_coord, args.outname, args.reads_file

    trans_act_genes = []
    with open(trans_act_file) as f:
        reader = csv.reader(f, delimiter = "\t")
        for line in reader:
            trans_act_genes.append(line[3])

    exons = rw.read_gtf(gtf, "exon")
    CDSs = rw.read_gtf(gtf, "CDS")

    exons = {i: exons[i] for i in exons if i in trans_act_genes}
    # protein-coding only
    exons = {i: exons[i] for i in exons if i in CDSs}

    ds_500 = "{0}_ds_500.bed".format(outname[:-4])
    with open(outname, "w") as out, open(ds_500, "w") as out_ds:
        writer = csv.writer(out, delimiter="\t")
        writer_ds = csv.writer(out_ds, delimiter="\t")
        for trans in exons:
            strand = exons[trans][0][6]
            chrom = "chr{0}".format(exons[trans][0][0])
            if strand == "+":
                TES = exons[trans][-1][4]
                new_start = TES - start_coord
                new_end = TES + end_coord
                new_start_ds = TES
                new_end_ds = TES + 500
            else:
                TES = exons[trans][-1][3]
                new_start = TES - start_coord - 1
                new_end = TES + start_coord - 1
                new_start_ds = TES - 500 - 1
                new_end_ds = TES - 1
            writer.writerow([chrom, new_start, new_end, trans, "0", strand])
            chrom = chrom.lstrip("chr")
            writer_ds.writerow([chrom, new_start_ds, new_end_ds, trans, "0", strand])

    intersect = "{0}_ds500_intersect.bed".format(outname[:-4])
    transcripts_file = "{0}_transcripts.bed".format(gtf[:-4])
    co.intersect_bed(ds_500, transcripts_file, write_both = True, force_strand=False, no_dups = False, output_file=intersect)

    co.get_transcripts(gtf, transcripts_file, with_detail=True)
    mapping = co.transcript_mapping(transcripts_file)

    to_exclude = []
    with open(intersect) as int_file:
        reader = csv.reader(int_file, delimiter = "\t")
        for line in reader:
            strand = line[5]
            curr_gene = mapping[line[3]]
            other_gene = mapping[line[9]]
            if curr_gene != other_gene:
                to_exclude.append(line[3])

    filtered_out_name = "{0}_filt.txt".format(outname[:-4])
    with open(filtered_out_name, "w") as filt_f:
        for name in to_exclude:
            filt_f.write("{0}\n".format(name))

    final_out_name = "{0}_distrib.bed".format(outname[:-4])

    distances = co.peak_pos_in_exon(outname, reads_file, from_end = True, reads_mode = True)[0]
    write_dist_mat(distances, start_coord + end_coord, final_out_name, None, "{0}_names.txt".format(final_out_name[:-4]), None)


if __name__ == "__main__":
    main()