'''Module that contains functions for reading and writing files.'''

import csv
import housekeeping as hk
import re
import sys

def read_fasta(input_file, multiline = False, return_dict = False):
    """
    Given a fasta file return a first list containing the sequence identifiers and a
    second list containing the sequences (in the same order). Also possible to return a
    dictionary with sequence IDs as keys and sequences as values
    :param input_file: input FASTA
    :param multiline: if True, the FASTA has multiline sequences
    :param return_dict: return a dictionary rather than two lists
    :return: either two lists or a dictionary
    """
    with open(input_file) as file:
        input_lines = file.readlines()
    #if the sequences are split over multiple lines
    if multiline:
        for line in range(len(input_lines) - 1):
            if input_lines[line][0] != ">" and input_lines[line + 1][0] != ">":
                input_lines[line] = input_lines[line].rstrip("\n")
        input_lines = "".join(input_lines)
        input_lines = input_lines.split("\n")
    #remove blank lines
    #it must be > 1 rather than > 0 because of the newline character
    input_lines = [i.rstrip("\n") for i in input_lines if len(i) > 1]
    names = [i.lstrip(">") for i in input_lines if i[0] == ">"]
    sequences = [i for i in input_lines if i[0] != ">"]
    if len(sequences) != len(names):
        print("Problem extracting data from fasta file!")
        print(len(sequences))
        print(len(names))
        raise Exception
    if len(sequences) == 0:
        print("No sequences were extracted!")
        raise Exception
    if return_dict:
        dictionary = {names[pos]: sequences[pos] for pos in range(len(names))}
        return(dictionary)
    else:
        return(names, sequences)

def read_as_string(file_name):
    '''
    Open a file, read in the contents as a string, close the file.
    :param file_name: file name
    :return: output string
    '''
    with open(file_name) as file:
        return("".join(file))

def read_gtf(file_name, element, gene = False, filter_parameter = None, filter_value = None):
    '''
    Read in all rows that contain coordinates of type _element_ from gtf file.
    Will make a dictionary with either gene or transcript IDs as keys
    (depending on what _gene_ is set to).
    Will also convert start and end coordinates to integers.
    If filter_parameter and filter-value have been specified, only lines with the specified value for
    that parameter will be returned.
    '''
    filter_pattern = ""
    if filter_parameter:
        if not filter_value:
            raise Exception("If filter_parameter has been specified, then filter_value must be too!")
        filter_pattern = "{0} \"{1}\"".format(filter_parameter, filter_value)
    output = {}
    if gene:
        pattern = re.compile("(?<=gene_id \")[\d\w]*")
    else:
        pattern = re.compile("(?<=transcript_id \")[\d\w]*")
    # check if you're on linux or Mac to know whether you need the -P flag
    platform = sys.platform
    if platform == "linux" or platform == "linux2":
        relevant_lines = hk.run_process(["grep", "-P", r"\t{0}\t".format(element), file_name]).rstrip("\n").split("\n")
    elif platform == "darwin":
        relevant_lines = hk.run_process(["grep", r"\t{0}\t".format(element), file_name]).rstrip("\n").split("\n")
    for line in relevant_lines:
        line = line.split("\t")
        if len(line) > 1:
            # if you need to filter by parameter
            if filter_pattern in line[8]:
                # convert start and end coordinates to integers
                line[3] = int(line[3])
                line[4] = int(line[4])
                # get identifier (transcript or gene)
                idn = re.search(pattern, line[8]).group()
                if idn not in output:
                    output[idn] = []
                output[idn].append(line)
    return(output)

def read_many_fields(input_file, delimiter, skip_header = False):
    '''
    Read a csv/tsv/... into a list of lists with each sublist corresponding to one line.
    '''
    file_to_read = open(input_file)
    field_reader = csv.reader(file_to_read, delimiter = delimiter)
    lines = []
    for i in field_reader:
        lines.append(i)
    file_to_read.close()
    if skip_header:
        lines = lines[1:]
    return(lines)

def write_list(list_to_write, file_name):
    '''
    Given a list, write each element to a separate line of an output file.
    :param list_to_write: list
    :param file_name: name of output file
    :return: None
    '''
    list_to_write = [str(i) for i in list_to_write]
    with open(file_name, "w") as file:
        for element in list_to_write[:-1]:
            file.write("{0}\n".format(element))
        file.write(list_to_write[-1])

def write_fasta(seqs, names, outfile_name):
    """
    Write a set of sequences to FASTA format.
    :param seqs: list of strings
    :param names: list of strings
    :param outfile_name: string
    :return: None
    """
    with open(outfile_name, "w") as ofile:
        for pos, seq in enumerate(seqs):
            ofile.write(">{0}\n".format(names[pos]))
            ofile.write("{0}\n".format(seq))

def write_gtf(indict, outfile, add_chr=False):
    """
    Given a dictionary where the values are lists of lists of lines from a GTF file,
    write them to a new GTF file.
    :param indict: dictionary of lists of lists of GTF lines
    :param outfile: output file name (string)
    :param add_chr: if True, append "chr" to chromosome names
    :return: None
    """
    with open(outfile, "w") as file:
        for element in indict:
            for line in indict[element]:
                if add_chr:
                    # only modify chromosomes
                    if line[0][0] not in ["K", "G"]:
                        line[0] = "chr{0}".format(line[0])
                file.write("\t".join([str(i) for i in line]))
                file.write("\n")

def write_many_fields(to_write, file_name, delimiter):
    """
    Write a list of lists to a delimited file.
    :param to_write: input list
    :param file_name: file to write the list of lists to
    :param delimiter: field separator
    :return: None
    """
    with open(file_name, "w") as file:
        field_writer = csv.writer(file, delimiter = delimiter)
        for line in to_write:
            field_writer.writerow(line)
