'''
Author: Rosina Savisaar.
Module that contains generic utility functions that make life a bit easier.
'''

import argparse
import itertools as it
import multiprocessing as mp
import os
import subprocess

def add_key(key, value, dict):
    """
    Check if a key is already present in a dictionary and, if not,
    add it with the specified value.
    :param key: key to add
    :param value: value for key
    :param dict: dictionary
    :return: dictionary with key added
    """
    if key not in dict:
        dict[key] = value
    return dict

def convert2bed(input_file_name, output_file_name):
    '''
    Converts an input file (sam, bam, gtf, gff...) to a bed file using bedops.
    '''
    extension = get_extension(input_file_name, 3, ["sam", "bam", "gtf", "gff"])
    bed_data = run_process(["convert2bed", "--input={0}".format(extension)],
                           file_for_input = input_file_name,
                           file_for_output = output_file_name)
    print("Converted data from {0} to bed.".format(extension))

def flatten(structured_list):
    """
    Flatten a structured list.
    :param structured_list: a list of lists
    :return: a flat list
    """
    flat_list = list(it.chain(*structured_list))
    return(flat_list)

def get_extension(file_name, extension_length, valid_list = None):
    '''
    Determine the extension at the end of a file name.
    '''
    extension = file_name[-extension_length:]
    if valid_list:
        if extension not in valid_list:
            raise Exception("File format must be included in {0}!".format(valid_list))
    return(extension)

def line_count(file):
    '''
    Count the number of lines in a file.
    '''
    #not using wc -l because I want the number of lines, not the number of newlines.
    output = run_process(["grep", "-c", "^", file])
    return(int(output))

def list_to_dict(input_list, key_index, value_index, floatify = False, intify = False):
    """
    Convert list to dictionary
    :param input_list: list to convert
    :param key_index: index for dictionary keys
    :param value_index: index for dictionary values
    :param floatify: convert value to float
    :param intify: convert value to integer
    :return: dictionary
    """
    out_dict = {}
    for elem in input_list:
        if floatify:
            elem[value_index] = float(elem[value_index])
        if intify:
            elem[value_index] = int(elem[value_index])
        out_dict[elem[key_index]] = elem[value_index]
    return(out_dict)

def make_dir(path):
    '''
    Create new directory if it doesn't already exist
    '''
    if not os.path.exists(path):
        os.mkdir(path)

def parse_arguments(description, arguments, floats = None, flags = None, defaults = None, ints = None, detailed_help = None):
    '''
    Use argparse to parse a set of input arguments from the command line.
    '''
    if not defaults:
        defaults = {}
    if not floats:
        floats = []
    if not flags:
        flags = []
    if not ints:
        ints = []
    parser = argparse.ArgumentParser(description = description)
    for pos, argument in enumerate(arguments):
        if detailed_help:
            help_info = detailed_help[pos]
        else:
            help_info = argument
        if pos in floats:
            curr_type = float
        elif pos in ints:
            curr_type = int
        else:
            curr_type = str
        if pos in flags:
            if pos in defaults:
                parser.add_argument("--{0}".format(argument), action="store", type = curr_type, help=help_info, nargs='?', default=defaults[pos])
            else:
                parser.add_argument("--{0}".format(argument), action="store_true", help=help_info)
        else:
            parser.add_argument(argument, type = curr_type, help = help_info)
    args = parser.parse_args()
    return(args)

def print_elements(input_list):
    '''
    Take a list and print out the elements separated by carriage returns.
    '''
    for i in input_list:
        print(i)
    print("\n")

def remove_file(file_name):
    '''
    Remove a file, if it exists.
    '''
    try:
        os.remove(file_name)
    except FileNotFoundError:
        pass

def run_in_parallel(input_list, args, func, kwargs_dict = None, workers = None, onebyone = False):
    """
    Take an input list, divide into chunks and then apply a function to each of the chunks in parallel.
    :param input_list: a list of the stuff you want to parallelize over (for example a list of gene names)
    :param args: a list of arguments to the function
    :param func: the function
    :param kwargs_dict: a dictionary of any keyword arguments the function might take
    :param workers: number of parallel processes to launch
    :param onebyone: if True, allocate one element from input_list to each process
    :return:
    """
    if not workers:
        #divide by two to get the number of physical cores
        #subtract one to leave one core free
        workers = int(os.cpu_count()/2 - 1)
    elif workers == "all":
        workers = os.cpu_count()
    #in the list of arguments, I put in "foo" for the argument that corresponds to whatever is in the input_list
    # because I couldn't be bothered to do something less stupid
    arg_to_parallelize = args.index("foo")
    if not onebyone:
        #divide input_list into as many chunks as you're going to have processes
        chunk_list = [input_list[i::workers] for i in range(workers)]
    else:
        #each element in the input list will constitute a chunk of its own.
        chunk_list = input_list
    pool = mp.Pool(workers)
    results = []
    #go over the chunks you made and laucnh a process for each
    for i in chunk_list:
        current_args = args.copy()
        current_args[arg_to_parallelize] = i
        if kwargs_dict:
            process = pool.apply_async(func, tuple(current_args), kwargs_dict)
        else:
            process = pool.apply_async(func, tuple(current_args))
        results.append(process)
    pool.close()
    pool.join()
    return(results)

def run_process(arguments, return_string = True, input_to_pipe = None, return_error = False, file_for_input = None, file_for_output = None, univ_nl = True, shell = False, verbose = False):
    '''
    Run a command on the command line.
    '''
    if file_for_input:
        input_file = open(file_for_input)
        stdin_src = input_file
    else:
        stdin_src = subprocess.PIPE
    if file_for_output:
        output_file = open(file_for_output, "w")
        stdout_dest = output_file
    else:
        stdout_dest = subprocess.PIPE
    arguments = [str(i) for i in arguments]
    if shell:
        arguments = " ".join(arguments)
    if verbose:
        print("Running command: {0}".format("".join(arguments)))
    process = subprocess.Popen(arguments, shell = shell, stdout = stdout_dest, stderr = subprocess.PIPE,
                               stdin = stdin_src, universal_newlines = univ_nl)
    if input_to_pipe:
        stdout, stderr = process.communicate(input_to_pipe)
    else:
        stdout, stderr = process.communicate()
    if file_for_input:
        input_file.close()
    if file_for_output:
        output_file.close()
    return_code = process.poll()
    if return_code != 0:
        print("Process failed!")
        print(" ".join(arguments))
        print(stderr)
        return("error")
    #if the process returns bytes but you want to get a string back.
    if return_string and type(stdout) == bytes:
        stdout = stdout.decode("utf-8")
    if return_error:
        return(stderr)
    else:
        return(stdout)

def update_counter(counter, step, message = None):
    '''
    Print out and update counter.
    '''
    if counter % step == 0:
        if message:
            print(message.format(counter))
        else:
            print(counter)
    counter = counter + 1
    return(counter)