from collections import Counter
import coord_ops as co
import csv
import housekeeping as hk

def main():
    description = "Write out a BED file with the region surrounding the TSS for a set of genes."
    args = hk.parse_arguments(description, ["genes_file", "gtf", "outfile", "start_coord", "end_coord"], ints = [3,4])
    genes_file, gtf, outfile, start_coord, end_coord = args.genes_file, args.gtf, args.outfile, args.start_coord, args.end_coord

    need_to_seek = False
    if genes_file[-3:] != "bed":
        # it means that you got a list of gene symbols rather than a
        # BED file with coordinates
        need_to_seek = True
        transcript_file = "{0}_transcripts.gtf".format(gtf[:-4])
        co.get_transcripts(gtf, transcript_file, add_chr=False, with_detail=False, output_gtf=True)


    with open(genes_file) as gf, open(outfile, "w") as of:
        reader = csv.reader(gf, delimiter = "\t")
        writer = csv.writer(of, delimiter = "\t")
        for line in reader:
            if need_to_seek:
                gene = line[0]
                print(gene)
                possibilities = hk.run_process(["grep", "gene_symbol \"\"{0}\"\"".format(gene), transcript_file])
                possibilities = [i.split("\t") for i in possibilities.split("\n")[:-1]]
                chrom = "chr{0}".format(possibilities[0][0])
                strand = possibilities[0][6]
                starts = [i[3] for i in possibilities]
                ends = [i[4] for i in possibilities]
                if strand == "+":
                    counts = Counter(starts)
                    start = int(counts.most_common()[0][0])
                    end = ends[starts.index(str(start))]
                    length = int(end) - start
                    new_end_coord = min(length, end_coord)
                    new_start = str(start - start_coord - 1)
                    new_end = str(start + (new_end_coord - 1))
                elif strand == "-":
                    counts = Counter(ends)
                    end = int(counts.most_common()[0][0])
                    start = starts[ends.index(str(end))]
                    length = end - int(start)
                    new_end_coord = min(length, end_coord)
                    new_start = str(end - new_end_coord)
                    new_end = str(end + start_coord)
                new_line = [chrom, new_start, new_end, gene, ".", strand]
                writer.writerow(new_line)
            else:
                if line[0] != "chrom":
                    length = int(line[2]) - int(line[1])
                    curr_end_coord = min(length, end_coord)
                    if line[5] == "+":
                        new_start = str(int(line[1]) - start_coord)
                        new_end = str(int(line[1]) + curr_end_coord)
                    elif line[5] == "-":
                        new_start = str(int(line[2]) - curr_end_coord)
                        new_end = str(int(line[2]) + start_coord)
                    else:
                        raise Exception("Invalid strand!")
                    new_line = line.copy()
                    new_line[1] = new_start
                    new_line[2] = new_end
                    # to make it a BED6
                    new_line = new_line[:-2]
                    writer.writerow(new_line)



if __name__ == "__main__":
    main()