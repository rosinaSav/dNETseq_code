This repository contains some of the Python code used to perform the analyses reported
in the dNET-seq paper (reference TBA). (The remainder of the Python code can be found in
https://github.com/kennyrebelo). Importantly, this repository is intended merely as
documentation of the methods. I haven't checked that the code runs on
systems other than my own and the documentation of dependencies is not systematic
(although an effort was made to make peak_caller.py a little more user-friendly than
the others). The code was developed using Python3.7 on a Mac (10.14.6).

However, if anybody does want to use any of this code and runs into problems,
please do drop me an e-mail (rosinasavisaar@gmail.com) and we'll figure it out!

###############################################

---Main scripts---

check_5prime_bias.py: Make a PPM of the sequence just around the 5' ends
of NET-seq reads.

filter_out_intron_lariats_and_splice_intermediates.py: Given a BED file of NET-seq
(or other) reads, filter out reads that end at the
last nucleotide of an exon (putative splicing intermediates)
or the last nucleotide of an intron (putative intron lariats).

make_stats_for_splice_sites.py: Make a file with various gene architecture
statistics for a set of introns.

mnase_bias.py: Create a control set that has the same nucleotide composition
at 5' ends as a true NET-seq dataset.

peak_caller.py: Peak caller for detecting regions where NET-seq read density is significantly
higher than expected by chance for a given transcript.

peak_distrib.py: Build meta-profile of peak/read densities.

splice_distance.py: Identify spliced.unspliced reads and build a meta-profile for either.

---Modules---

These are modules that contain sets of thematically related functions, called
upon in the main script. Some of the functions contained in these modules
are not actually used in the main scripts (because of analyses that didn't
make it into the final version of the paper.)

coord_ops.py: Various operations on genomic coordinates.

housekeeping.py: Generic utility functions that make life a bit easier.

NGS.py: Functions to help process NGS, especially NET-seq data.

nucleotide_comp.py: Analysis of nucleotide/motif composition.

read_and_write.py: Reading and writing files.

---Tests---

The repository also contains test scripts intended to be run using the
Python unittest module. The directory _tests_ contains various test input
and output files that are necessary to run the tests.
