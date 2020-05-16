'''
Make a PPM of the sequence just around the 5' ends of NET-seq reads.
'''

import coord_ops as co
import housekeeping as hk
import nucleotide_comp as nc
import numpy as np

def extract_true_and_control_string(fasta_name, true_indices, control_indices):
    """
    Given a FASTA file and two tuples of indices, make two NUMPY arrays with the sequence
    between the indices in each of the FASTA entries.
    :param fasta_name: input FASTA file
    :param true_indices: tuple indicating the start and end (0-based) of the first segment
    :param control_indices: tuple giving the start and end of the second segment
    :return: two numpy arrays with rows corresponding to FASTA entries and columns to positions
    in the segment.
    """
    expected_length_true = true_indices[1] - true_indices[0]
    expected_length_control = control_indices[1] - control_indices[0]
    # Get the number of lines in FASTA and divide by 2 to get the number of sequences
    fasta_length = hk.line_count(fasta_name)/2
    # Pre-allocate two arrays, one for the true sequence and one for the control
    occ_mat_true = np.empty((int(fasta_length), int(true_indices[1] - true_indices[0])), dtype="str")
    occ_mat_control = np.empty((int(fasta_length), int(control_indices[1] - control_indices[0])), dtype="str")
    pos_in_fasta = 0
    error_counter = 0
    with open(fasta_name) as fasta:
        for line in fasta:
            if line[0] != ">":
                true_string = line[true_indices[0]:true_indices[1]]
                control_string = line[control_indices[0]:control_indices[1]]
                if (len(true_string) != expected_length_true) or (len(control_string) != expected_length_control):
                    error_counter = error_counter + 1
                for pos in range(len(true_string)):
                    occ_mat_true[pos_in_fasta, pos] = true_string[pos]
                    occ_mat_control[pos_in_fasta, pos] = control_string[pos]
                pos_in_fasta = pos_in_fasta + 1
    print("Errors: {0}.".format(error_counter))
    return(occ_mat_true, occ_mat_control)

def PPM_wrapper(occ_mat, bases, output_file):
    """
    Given an occurrence matrix, make a PPM, print it to screen and save it to file.
    :param occ_mat: 2D numpy array with sequences in rows and positions in columns
    :param bases: bases to be considered
    :param output_file: output file name
    :return: PPM as 2D numpy array
    """
    print("Matrix shape:")
    print(np.shape(occ_mat))
    PPM = nc.make_PPM(occ_mat, bases)
    print("PPM:")
    print("\t".join(bases))
    print(PPM)
    print("\n")
    print("Average base composition:")
    print("\t".join(bases))
    print(np.sum(PPM, axis=0)/np.shape(PPM)[0])
    print("\n")
    np.savetxt(output_file, PPM, delimiter="\t")
    return(PPM)

def main():

    # Get arguments.
    description = "Check if nucleotide composition at the 5' ends of NET-seq reads is biased."
    args = hk.parse_arguments(description, ["input_file", "output_file", "genome_fasta", "gtf", "three_prime"], flags = [4])
    input_file, output_file, genome_fasta, gtf, three_prime = args.input_file, args.output_file, args.genome_fasta, args.gtf, args.three_prime

    # Convert to .bed, if not already .bed
    if input_file[-3:] != "bed":
        print("Converting input file to .bed...")
        input_file_new_name = "{0}bed".format(input_file[:-3])
        hk.convert2bed(input_file, input_file_new_name)
        input_file = input_file_new_name

    # Make an extended version of each read that extends 5 nt 5prime and 35 nt 3prime
    print("Extending the reads...")
    suffix = ""
    if three_prime:
        suffix = "_three_prime"
    temp_bed = "{0}_extended_for_bias{1}.bed".format(input_file[:-4], suffix)
    co.extend_intervals(input_file, temp_bed, 5, 35, remove_chr=True, add_chr=False, three_prime=three_prime)

    # Make a FASTA file from the BED file.
    print("Extracting sequences...")
    fasta_name = "{0}fasta".format(temp_bed[:-3])
    hk.run_process(["fastaFromBed", "-bed", temp_bed, "-fi", genome_fasta, "-fo", fasta_name, "-s"])
    print("Number of lines in FASTA:")
    print(hk.run_process(["wc", "-l", fasta_name]))

    # Store the sequences at -5:+5 and 30:40 in a 2D array
    print("Storing sequences in arrays...")
    occ_mat_true, occ_mat_control = extract_true_and_control_string(fasta_name, (0, 10), (30, 40))

    # Make a PPM for either column
    bases = ["A", "T", "C", "G"]
    print("Making PPMs...\n")
    print("TRUE:")
    PPM_wrapper(occ_mat_true, bases, "{0}.true".format(output_file))
    print("CONTROL:")
    PPM_wrapper(occ_mat_control, bases, "{0}.control".format(output_file))

if __name__ == "__main__":
    main()