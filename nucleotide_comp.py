'''Analysis of nucleotide/motif composition.'''

import coord_ops as co
import csv
import housekeeping as hk
import itertools as it
import numpy as np
import scipy.stats

def calc_nt_freqs(motifs, order, return_array = False, no_phase = False, count = False, alt_bases = None):
    '''
    Calculate the order-nucleotide frequencies of a set of motifs.
    no_phase if you only want to consider the first reading frame (no overlaps between kmers)
    returns numpy array if return_array, else returns dictionary
    '''
    all_kmers = generate_all_kmers(order, alt_bases=alt_bases)
    all_kmers_in_motifs = []
    if no_phase:
        phase_range = 1
    else:
        phase_range = order
    #loop over the different possible reading frames
    for phase in range(phase_range):
        for motif in motifs:
            current_strings = [motif[i:i+order] for i in range(phase, len(motif), order)]
            #so it wouldn't split the k-mer if it's halfway through the kmer when it reaches the end of the motif
            current_strings = [i for i in current_strings if (len(i) == order) and ("N" not in i)]
            all_kmers_in_motifs.extend(current_strings)
    kmer_number = len(all_kmers_in_motifs)
    if return_array:
        if count:
            result_array = np.array([all_kmers_in_motifs.count(kmer) for kmer in all_kmers])
        else:
            result_array = np.array([all_kmers_in_motifs.count(kmer)/kmer_number for kmer in all_kmers])
        if abs(np.sum(result_array) - 1) > 0.0001 :
            print(result_array)
            print("Frequencies don't sum up to 1!")
            raise Exception
        return(result_array)
    freqs_dict = {}
    if kmer_number:
        for kmer in all_kmers:
            if count:
                freqs_dict[kmer] = all_kmers_in_motifs.count(kmer)
            else:
                freqs_dict[kmer] = all_kmers_in_motifs.count(kmer)/kmer_number
    else:
        return(None)
    if not count:
        if abs(np.sum(list(freqs_dict.values())) - 1) > 0.0001:
                print(freqs_dict)
                print(motifs)
                print("Frequencies don't sum up to 1!")
                raise Exception
    return(freqs_dict)

def generate_all_kmers(k, alt_bases = False):
    """
    Generate all k-mers of a given k.
    :param k: motif length
    :param alt_bases: bases to use. If not specified, the four canonical bases will be used.
    :return: List of strings.
    """
    bases = ["A", "T", "C", "G"]
    if alt_bases:
        bases = alt_bases
    kmers = [''.join(i) for i in it.product(bases, repeat=k)]
    return(kmers)

def get_4fold_deg(sequence):
    """
    Determine where the fourfold degenerate positions are in a sequence.
    :param sequence: sequence string (starts and ends in a full codon)
    :return: list of fourfold degenerate position indices.
    """
    fourfold_pos = []
    for i in range(0,len(sequence),3):
        if sequence[i:i+3] in fourfold_deg_list:
            fourfold_pos.append(i+2)
    return(fourfold_pos)

def get_fourfold_deg_list():
    """
    Get a list of all the fourfold degenerate codons.
    :return: list of codons
    """
    codons = ["TGA","TAG","TAA","TTT", "TTC", "TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "ATT", "ATC", "ATA", "ATG", "GTT", "GTC", "GTA", "GTG", "TCT", "TCC", "TCA", "TCG", "CCT", "CCC", "CCA", "CCG", "ACT", "ACC", "ACA", "ACG", "GCT", "GCC", "GCA", "GCG", "TAT", "TAC", "CAT", "CAC", "CAA", "CAG", "AAT", "AAC", "AAA", "AAG", "GAT", "GAC", "GAA", "GAG", "TGT", "TGC", "TGG", "CGT", "CGC", "CGA", "CGG", "AGT", "AGC", "AGA", "AGG", "GGT", "GGC", "GGA", "GGG"]
    aas =  ["STOP","STOP","STOP","Phe", "Phe", "Leu", "Leu", "Leu", "Leu", "Leu", "Leu", "Ile", "Ile", "Ile", "Met", "Val", "Val", "Val", "Val", "Ser", "Ser", "Ser", "Ser", "Pro", "Pro", "Pro", "Pro", "Thr", "Thr", "Thr", "Thr", "Ala", "Ala", "Ala", "Ala", "Tyr", "Tyr", "His", "His", "Gln", "Gln", "Asn", "Asn", "Lys", "Lys", "Asp", "Asp", "Glu", "Glu", "Cys", "Cys", "Trp", "Arg", "Arg", "Arg", "Arg", "Ser", "Ser", "Arg", "Arg", "Gly", "Gly", "Gly", "Gly"]
    bases = ["A","T","C","G"]
    #create a list for storing all codons that have a fourfold degenerate site in their third position
    #loop over all the codons and store the first two bases as a string.
    #then loop over the four bases and append each of them in turn to the string containing the first two bases in the current codon.
    #check whether the amino-acid encoded for by the new codon is the same as that encoded for by the new codon you just made.
    #if this is the case for all four bases, add the current codon to the list you created in the beginning of this section.
    fourfold_deg_codons = []
    for i in range(len(codons)):
        first_two = "".join(codons[i][0:2])
        codon_change = False
        for j in bases:
            codon_comp = "".join([first_two,j])
            if aas[i] != aas[codons.index(codon_comp)]:
                codon_change = True
        if codon_change == False:
            fourfold_deg_codons.append(codons[i])
    return(fourfold_deg_codons)

def get_GC(sequence, alternative_bases = None):
    """
    Get the GC content of a sequence. Can also be used for other sets of bases than GC.
    :param sequence: string of nucleotides
    :param alternative_bases: list of characters to be used instead of G and C
    :return: GC content of sequence
    """
    sequence = sequence.upper()
    sequence = [i for i in sequence if i in ["A", "T", "C", "G"]]
    if alternative_bases:
        bases = alternative_bases
    else:
        bases = ["G","C"]
    counter = 0
    for base in bases:
        counter = counter + sequence.count(base)
    GC = counter/len(sequence)
    return(GC)

def get_GC4(sequence, phase, alternative_bases = None):
    """
    Get the GC content at the fourfold degenerate bases of a sequence.
    :param sequence: sequence string.
    :param phase: Phase of the sequence. Note that the phase should be entered
    like in a GTF file and is then converted internally into an in-house system.
    Phase 0 means that the sequence starts with the first base of a codon.
    Phase 1: the sequence starts with the third base of a codon (there is one extra base).
    Phase 2: the sequence starts with the second base of a codon (there are two extra bases).
    :param alternative_bases: List of characters to be used instead of G and C.
    :return: GC4 (a float)
    """
    #make sure the sequence starts and ends with full codons.
    sequence = co.trim_sequence(sequence, phase)
    #figure out where the 4-fold degenerate sites are.
    fourfold_deg_pos = get_4fold_deg(sequence)
    if not fourfold_deg_pos:
        return None
    if alternative_bases:
        bases = alternative_bases
    else:
        bases = ["G","C"]
    hits = 0
    for i in range(2,len(sequence),3):
        if sequence[i] in bases:
            if i in fourfold_deg_pos:
                hits = hits + 1
    GC4 = hits/len(fourfold_deg_pos)
    return(GC4)

def get_ss_strength(exons, genome_file, upstream = True, five = True, exonic = 3, intronic = 6):
    """
    Given a set of exons, get an estimate of splice site strength.
    :param exons: Dictionary of CDS lines.
    :param genome_file: File with genome sequence.
    :param upstream: evaluate the (5' or 3') splice site of the upstream intron (rather than downstream)
    :param five: evaluate the 5' splice site (rather than 3')
    :param exonic: how many nucleotides to include from the exon
    :param intronic: how many nucleotides to include from the intron
    :return: a dictionary with the splice site strength for each exon
    """
    # will contain the splice site strengths
    out_dict = {}
    # will contain the names of the exons so that later on, we'd know which
    # splice site strength value goes with which exon
    names = []

    # write splice site coordinates to GTF
    hk.make_dir("temp_data")
    temp_file_name = "temp_data/ss_sequences.gtf"
    with open(temp_file_name, "w") as temp_file:
        writer = csv.writer(temp_file, delimiter = "\t")
        for transcript in exons:
            curr_exons = exons[transcript]
            for pos, exon in enumerate(curr_exons):
                # don't analyze first exons
                if (pos != 0):
                    # cause you can't do the downstream intron of the last exon
                    if (upstream or (pos != len(curr_exons) - 1)):
                        if five:
                            if upstream:
                                template = curr_exons[pos - 1].copy()
                            else:
                                template = exon.copy()
                            if template[6] == "+":
                                template[3] = template[4] - exonic + 1
                                template[4] = template[4] + intronic
                            elif template[6] == "-":
                                template[4] = template[3] + exonic - 1
                                template[3] = template[3] - intronic
                        else:
                            if upstream:
                                template = exon.copy()
                            else:
                                template = curr_exons[pos + 1].copy()
                            if template[6] == "+":
                                template[4] = template[3] + exonic - 1
                                template[3] = template[3] - intronic
                            elif template[6] == "-":
                                template[3] = template[4] - exonic + 1
                                template[4] = template[4] + intronic
                        # this is for scaffolds etc.
                        if template[3] >= 0:
                            # so you'd know the order of the values in the MaxEntScan output
                            names.append("{0}.{1}".format(transcript, pos - 1))
                            writer.writerow(template)

    # make a FASTA with splice site sequences
    temp_fasta_file_name = "{0}.fasta".format(temp_file_name[:-4])
    hk.run_process(["bedtools", "getfasta", "-fi", genome_file, "-bed", temp_file_name, "-fo", temp_fasta_file_name, "-s"])
    # filter FASTA for Ns
    fasta_lines = []
    with open(temp_fasta_file_name) as fasta:
        for line in fasta:
            if line[0] == ">":
                curr_name = line
            else:
                if "N" not in line:
                    fasta_lines.append(curr_name)
                    fasta_lines.append(line)
    with open(temp_fasta_file_name, "w") as fasta:
        for line in fasta_lines:
            fasta.write(line)

    # run MaxEntScan on the FASTA
    # lazy hardcoded path, replace as appropriate...
    mes_direct = "/Users/rsavisaar/Software/MaxEntScan/fordownload"
    if five:
        cmd = "/Users/rsavisaar/Software/MaxEntScan/fordownload/score5.pl"
    else:
        cmd = "/Users/rsavisaar/Software/MaxEntScan/fordownload/score3.pl"
    temp_mes_file_name = "{0}_mes.txt".format(temp_file_name[:-4])
    hk.run_process(["perl", cmd, temp_fasta_file_name], file_for_output=temp_mes_file_name, verbose=True)
    hk.remove_file(temp_fasta_file_name)
    hk.remove_file(temp_file_name)

    # read in splice site scores and store in output directory
    with open(temp_mes_file_name, newline = "") as mes_file:
        reader = csv.reader(mes_file, delimiter = "\t")
        for pos, line in enumerate(reader):
            out_dict[names[pos]] = float(line[1])
    hk.remove_file(temp_mes_file_name)
    return(out_dict)

def kmer_enrichment_test(occurrence_number, probability, trials):
    """
    Performs a binomial test on whether or not a motif is occurring as expected by chance.
    :param occurrence_number: number of times the motif was observed
    :param probability: expected probability of occurrences
    :param trials: number of trials
    :return: a tuple consisting of the expetced number of occurrences, the ratio of observed to
    expected and a two-tailed binomial p-value
    """
    expected = trials * probability
    ratio = occurrence_number/expected
    p = scipy.stats.binom_test(occurrence_number, n=trials, p=probability, alternative="two-sided")
    return(expected, ratio, p)

def make_PPM(occ_mat, bases):
    """
    Given a numpy matrix of aligned sequences of equal length, make a position probability matrix
    (positions in rows, bases in columns).
    :param occ_mat: 2D numpy array
    :param bases: list of bases to consider
    :return: PPM (2D numpy array)
    """
    PPM = np.transpose(np.apply_along_axis(make_PPM_core, 0, occ_mat, bases))
    return(PPM)

def make_PPM_core(in_row, bases):
    """
    Given a sequence (represented as a 1D numpy array), count the numbers of each of the bases.
    :param in_row: 1D numpy array
    :return: base counts as 1D numpy array
    """
    out = np.array([np.count_nonzero(in_row == base) for base in bases])
    # divide by the sum rather than by the length of the row number of
    # occ_mat because you don't want to count anything that is not in bases
    # (i.e. Ns)
    out = np.divide(out, np.sum(out))
    return(out)

def markov(seq_frags, bases, k, coding = False):
    """
    Given a set of sequence fragments of equal length, build a 2nd-order Markov model of nucleotide transition
    probabilities. If the sequence is coding, also condition the model on position in the codon.
    :param seq_frags: list of strings
    :param bases: bases to consider
    :param k: order of Markov model
    :param coding: if True, the model will be conditioned also on codon position
    :return: two dictionaries, one with dinucleotide frequencies and one with the
    transition probabilities.
    """
    if coding:
        # only get the frequency at codon positions 1 and 2 (if k == 2, otherwise adjust)
        # NB! This doesn't necessarily mean that they are actually the first and second in a codon:
        # in this code, codon position 0 just means the same reading frame as the first position
        # in the fragment.
        # the k-mer frequencies dictionary is for final output, the k+1-mer one
        # is for calculating the transition probabilities
        start_seqs = hk.flatten([[i[j:j+k] for j in range(0, len(i), 3)] for i in seq_frags])
        freqs = calc_nt_freqs(start_seqs, order=k, alt_bases=bases)
        trans_probs = {i: {j: {k: 0 for k in range(3)} for j in bases} for i in freqs.keys()}
        # map where the k+1-mer will start if we're considering each end position
        frag_start_pos = {}
        for end_pos in range(3):
            frag_start_pos[end_pos] = (end_pos + 3 - k) % 3
        for pos in range(3):
            substrings = hk.flatten([[i[j:j+k+1] for j in range(frag_start_pos[pos], len(i), 3)] for i in seq_frags])
            substrings = [i for i in substrings if len(i) == (k + 1)]
            freqs_for_Markov = calc_nt_freqs(substrings, order=(k + 1), alt_bases=bases)
            curr_trans_probs = nt_freqs_to_transitions(freqs_for_Markov, bases)
            # insert values into correct places in the big dictionary
            for dint in curr_trans_probs:
                for base in curr_trans_probs[dint]:
                    trans_probs[dint][base][pos] = curr_trans_probs[dint][base]
    else:
        # just calculate the frequencies of the dinucleotides anywhere in the sequence fragments
        freqs = calc_nt_freqs(seq_frags, order=k, alt_bases=bases)
        freqs_for_Markov = calc_nt_freqs(seq_frags, order=(k + 1), alt_bases=bases)
        trans_probs = nt_freqs_to_transitions(freqs_for_Markov, bases)
    return(freqs, trans_probs)

def markov_prob(motif, start_freqs, trans_probs, k, coding=False):
    """
    Calculate the probability of a motif, given a Markov model.
    :param motif: string
    :param start_freqs: frequencies of motifs of order k - 1
    if k is the order of the Markov model. If the motif is from coding sequence,
    then only give the frequencies at the first two codon positions.
    :param trans_probs: Markov model as output by Markov()
    :param k: order of Markov model
    :param coding: if True, then the Markov model must be conditional
    on codon position
    :return: probability of the motif (float)
    """
    probs = np.empty(len(motif) - 1)
    probs[0] = start_freqs[motif[: k]]
    for codon_pos in range(len(motif)):
        # so you wouldn't get incomplete codons
        if codon_pos + 3 <= len(motif):
            curr_start = motif[codon_pos: codon_pos + k]
            if coding:
                probs[codon_pos + 1] = trans_probs[curr_start][motif[codon_pos + k]][(codon_pos + k)%3]
            else:
                probs[codon_pos + 1] = trans_probs[curr_start][motif[codon_pos + k]]
    prob = np.product(probs)
    return(prob)

def nt_freqs_to_transitions(frequencies, bases):
    """
    Convert a dictionary of k+1-mer frequencies into transition probabilities for a
    Markov model of order k.
    :param frequencies: dictionary where the keys are kmers and the
    values are relative frequencies of occurrence
    :param bases: nucleotide bases in alphabet
    :return:
    """
    # figure out k based on the frequency dictionary
    k = len(list(frequencies.keys())[0])
    starts = [i[: k - 1] for i in list(frequencies.keys())]
    transition_dict = {}
    # loop over all of the k-2 length fragments and make a temporary dictionary
    # for each
    for start in starts:
        transition_dict[start] = {}
        temp_dict = {}
        for base in bases:
            current_freq = frequencies["".join([start, base])]
            temp_dict[base] = current_freq
        total_fraction = sum(list(temp_dict.values()))
        if total_fraction == 0:
            transition_dict[start] = {base: 1/(len(bases)) for base in bases}
        else:
            transition_dict[start] = {base: temp_dict[base]/total_fraction for base in bases}
    return(transition_dict)

def PPM_diff(PPM1, PPM2):
    """
    Calculate the sum of squared differences between two PPMs.
    :param PPM1: 2D numpy array (first PPM)
    :param PPM2: 2D numpy array (second PPM)
    :return: float (sum of squared differences)
    """
    return(np.sum(np.square(np.subtract(PPM1, PPM2))))

fourfold_deg_list = get_fourfold_deg_list()
