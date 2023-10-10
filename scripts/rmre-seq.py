#!/usr/bin/env python
import sys
import os
import re
import argparse

# 2023 Ben Perry
# MIT License
# Copyright (c) 2023 Ben Perry
# Version: 1.0
# Email: ben.perry@agresearch.co.nz

def get_options():
    
    description = 'Utilities for processing, merging, and summarising rmre-seq data in the form of aligned .bed files'
    
    # Setting up the main parser object
    main_parser = argparse.ArgumentParser(description = description)
    
    utilities_parser = main_parser.add_subparsers(title = "utilities", 
                                                  dest = 'task')

    parser_methyl = utilities_parser.add_parser("bed-to-CH4",
                                        description="Convert alignment .bed file into .CH4.bed file")
    parser_methyl.add_argument('-i', '--input',
                        dest = 'bed_input',
                        type = str,
                        default = None,
                        required = True,
                        help = 'input .bed file or rmre-seq aligned reads.')
    parser_methyl.add_argument('-o', '--outfile',
                               dest = 'outfile',
                               type = str,
                               default = None,
                               required = True,
                               help = 'output file path for CH4 formatted .CH4.bed file.')
    

    parser_merge = utilities_parser.add_parser("merge-CH4",
                                        description="Merge multiple .CH4.bed files with/without a reference 'CCGG' gff set.")
    parser_merge.add_argument('-r', '--reference',
                        dest = 'reference',
                        type = str,
                        default = None,
                        required = False,
                        help = '''
                                reference "CCGG" motif gff file. Assumes generated with emboss fuzznuc and takes the format:
                                NC_056054.1     fuzznuc nucleotide_motif        85      88      4       +       .       ID=NC_056054.1.1;pattern=CCGG
                                ''')
    parser_merge.add_argument('-o', '--outfile',
                        dest = 'outfile',
                        type = str,
                        required = True,
                        help = 'output file path for the merged .CH4.bed.txt matrix.')
    
    
    args = main_parser.parse_args()

    return args



# def splice_sample_names(args):
#     """
#     for each file, parse the filename to get sampleID and locusID,
#     then parse the file contents and splice these in and list, print
#     updated summary table to stdout.
#     """
#     for path in args["input_files"]:
#         # parse the filename to get sampleID and locusID parts
#         # e.g. from HERB01-A07-J11-65-14_S27_EPSPS.summary.txt we want to be
#         # able to extract HERB01-A07-J11-65-14  EPSPS try each regular
#         # expression (not case-sensitive) and use the first match
#         # (usually there will be only one regexp and usually the default)
#         sample_name = None
#         locus = None

#         for regexp in args["regexps"]:
#             match = re.search(regexp, os.path.basename(path), re.IGNORECASE)
#             if match is not None:
#                 if len(match.groups()) == 2:
#                     (sample_name, locus) = match.groups()

#         if sample_name is None or locus is None:
#             raise Exception(
#                 "Sorry could not parse sampleID and locusID from %s using any of : %s" % (
#                     path, str(args["regexps"])))

#         # now open file, parse contents and splice in sample name and locus
#         # file contents are like :
#         # more /dataset/gseq_processing/active/bin/gtseq_prism/unit_test/HERB01-A07-J11-65-14_S27_ACCase1999.summary.txt
#         # 0.0000  ATTCCCATGAGCGGTCTGTTCCTCGTGCTGGGCAAGTCTGGTTTCCAGATTCTGCTACCAAGACAGCGCAGGCAATGTTGGACTTCAA;size=2072      *       *       *       *       *       *       *       *       0       0   00       0       0       *       N
#         # 0.0000  ATTCCCATGAGCGGTCTGTTTCTCGTGCTGGGCAAGTCTGGTTTCCAGATTCTGCTACCAAGACAGCGCAGGCAATGTTGGACTTCAA;size=42        *       *       *       *       *       *       *       *       0       0   00       0       0       *       N
#         with open(path, "r") as savs:
#             record_count = 0
#             field_count = 0
#             for record in savs:
#                 fields = re.split("\t", record.strip())
#                 if record_count == 0:
#                     field_count = len(fields)
#                 record_count += 1
#                 # sanity check file - records should all be the same length
#                 if len(fields) != field_count:
#                     raise Exception(
#                         "oops : in %s, first rec had %d fields, but record %d has %d fields" % (
#                             path, field_count, len(fields)))
#                 allele_and_size = re.split(";", fields[1])
#                 if len(allele_and_size) != 2:
#                     raise Exception(
#                         "oops : malformed sampleID field (%s) in record %d" %
#                         fields[1], record_count)
#                 (allele, size) = allele_and_size
#                 # output record , e.g. like CCCGATTGAGAAGGATGCCAAGGAGGAAGTAAAGCTCTTCTTGGGCAACGCTGGAACTGCAATGCGGCCATTGACGGCAGCTGTAGTAGCTGCTGGTGGA  size=1125  HERB01-A07-J11-65-14  EPSPS
#                 print(allele, size, sample_name, locus, sep="\t")


def main():
    get_options()

    # if opts["task"] == "splice_sample_names":
    #     splice_sample_names(opts)
    # else:
    #     raise Exception(
    #         "Oops... --task %(splice_sample_names)s is not supported yet !" % opts)
    return 0


if __name__ == "__main__":
    main()
