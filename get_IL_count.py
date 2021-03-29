import coord_ops as co
import housekeeping as hk

def main():
    description = "Make a file with intron lariat read counts per exon."
    args = hk.parse_arguments(description, ["intron_lariat_file", "regions_file"])
    intron_lariat_file, regions_file = args.intron_lariat_file, args.regions_file

    # the intron_lariat_file contains only those reads whose
    # 3' ends map to the last position of an intron
    snr_name = "{0}_snr.bed".format(intron_lariat_file[:-4])
    co.snr_bed(intron_lariat_file, snr_name)

    co.intersect_bed(regions_file, snr_name, force_strand=True, hit_count=True, no_dups=False, output_file="{0}_il_counts.bed".format(regions_file[:-4]))

if __name__ == "__main__":
    main()