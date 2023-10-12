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
    
    description = "Utilities for processing, merging, and summarising rmre-seq data in the form of aligned .bed files"
    
    # Setting up the main parser object
    main_parser = argparse.ArgumentParser(description = description)
    
    # initiat the subparser
    utilities_parser = main_parser.add_subparsers(title = "utilities", 
                                                  dest = "task")

    # bed-to-CH4 arg parser
    parser_CH4 = utilities_parser.add_parser("bed-to-CH4",
                                        description = "Convert alignment .bed file into .CH4.bed file",
                                        help = "Convert alignment .bed file into .CH4.bed file")
    parser_CH4.add_argument("-i", "--input",
                        dest = "bed_input",
                        type = str,
                        default = None,
                        required = True,
                        help = "input .bed file or rmre-seq aligned reads.")
    parser_CH4.add_argument("-o", "--outfile",
                               dest = "CH4_out",
                               type = str,
                               default = None,
                               required = True,
                               help = "output file path for CH4 formatted .CH4.bed file.")
    parser_CH4.add_argument("-d", "--header", 
                            dest = "header_option", 
                            action='store_true',
                            help = "include header information in output file.")
    parser_CH4.add_argument("-q", "--mapq_score", 
                            dest = "mapq_score", 
                            default = 0,
                            type = int,
                            help = "score (mapq) threshold for reporting read greater than (default = 0). Note: premapped read counts are used for calculating proportions.")
    
    # merge-CH4 arg parser
    parser_merge = utilities_parser.add_parser("merge-CH4",
                                        description="Merge multiple .CH4.bed files with/without a reference 'CCGG' gff set.",
                                        help = "Merge multiple .CH4.bed files with/without a reference 'CCGG' gff set.")
    parser_merge.add_argument("-R", "--reference",
                        dest = "reference",
                        type = str,
                        default = None,
                        required = False,
                        help = '''
                                reference "CCGG" motif gff file. Assumes generated with emboss fuzznuc and takes the format:\n
                                NC_056054.1     fuzznuc nucleotide_motif        85      88      4       +       .       ID=NC_056054.1.1;pattern=CCGG
                                ''')
    parser_merge.add_argument("-i", "--CH4bed",
                        dest = "input",
                        required = True,
                        nargs = "*",
                        help = "input .CH4.bed files; one or more required.")
    parser_merge.add_argument("-o", "--outfile",
                        dest = "outfile",
                        type = str,
                        required = True,
                        help = "output file path for the merged .CH4.bed.txt matrix.")
    parser_merge.add_argument("-t", "--report",
                        dest = "report",
                        type = str,
                        required = True,
                        choices=['proportions', 'counts'],
                        help = "values to report in .CH4.bed.txt matrix.")

    
    # check for no args and print help
    if len(sys.argv)==1:
        main_parser.print_help()
        # parser.print_usage() # for just the usage line
        sys.exit(1)
        
    # parse it
    args = main_parser.parse_args()
    print()
    
    # ship it
    return args


def bedtoCH4(bed_file_path, threshold):
    """ ("alignment.bed") -> CH4bed[]
    :param .bed file: String; path to an input .bed file.
    :return: .CH4.bed dataframe 
    
    Parses a .bed file from rmre-seq alignment. Then returns a dataframe of
    all unique positions and their counts.
    """
    import pandas as pd
    # read in gzipped bed file
    print()
     
    bed_file = pd.read_csv(bed_file_path,
                           compression = "gzip",
                           sep = "\t",
                           names = ["chrom", "startChrom", "endChrom", "name", "score", "strand"],
                           dtype = {"chrom":str, "startChrom":int, "endChrom":int, "name":str, "score":int, "strand":str},
                           comment = "#")
    
    # log some metrics
    total_reads = len(bed_file)
    print("total reads in .bed file: " + str(total_reads))
    
    score_0 = len(bed_file[bed_file["score"] == 0])
    print("reads filtered due to mapping score == 0 : " + str(score_0))
    
    bed_file = bed_file[bed_file["score"] > threshold ]
    print("reads filtered with map score less than " + str(threshold) + ": " + str(total_reads - len(bed_file)))

    print("reads passing filter: " + str(len(bed_file)))
    
    
    # split by strand
    bed_file_plus = bed_file[bed_file["strand"] == "+"]
    bed_file_minus = bed_file[bed_file["strand"] == "-"]
    
    # drop full file
    del(bed_file)

    # uniq counts of startChrom for '+'
    CH4_file_plus = bed_file_plus.groupby(["chrom", "startChrom"]).size().reset_index(name='count')
    CH4_file_plus["strand"] = "+"
    CH4_file_plus["endChrom"] = CH4_file_plus["startChrom"] + 1
    CH4_file_plus["name"] = CH4_file_plus["chrom"].transform(str) + "_" + CH4_file_plus["startChrom"].transform(str)
    
    # uniq counts of endChrom for '-'
    CH4_file_minus = bed_file_minus.groupby(["chrom", "endChrom"]).size().reset_index(name='count')
    CH4_file_minus["strand"] = "-"
    CH4_file_minus["startChrom"] = CH4_file_minus["endChrom"] - 1
    CH4_file_minus["name"] = CH4_file_minus["chrom"].transform(str) + "_" + CH4_file_minus["startChrom"].transform(str)
    
    # append counts back together
    CH4_file = pd.concat([CH4_file_plus, CH4_file_minus], ignore_index=True, axis=0)
    
    # add "total" and proportion columns
    CH4_file["total"] = CH4_file["count"].sum()
    CH4_file["proportion"] = CH4_file["count"] / CH4_file["total"]
    
    # reorder and sort CH4_file
    CH4_file = CH4_file[["chrom", "startChrom", "endChrom", "name", "count", "strand", "total", "proportion"]]
    CH4_file = CH4_file.sort_values(by = ["chrom", "startChrom", "endChrom"], ascending = True)
    
    # ship it
    return CH4_file

def mergeCH4():
    #TODO: make this.
    return

def main():
    # parse CLI
    args = get_options()

    # check 'task' and run accordingly.
    if args.task == "bed-to-CH4":
        # parse bed and return .CH4.bed formatted datframe
        print("running bed-to-CH4 on: " + str(args.bed_input))
        methylation = bedtoCH4(bed_file_path = args.bed_input, threshold = args.mapq_score)
        
        # print .CH4.bed file
        print("writing .CH4.bed file to: " + str(args.CH4_out))
        methylation.to_csv(args.CH4_out, sep = "\t", header = args.header_option, compression = "gzip", index = False, float_format="%.10f")
        
        print("""
MIT License
copyright (c) 2023 Benjamin J Perry
version: 1.0.0
email: ben.perry@agresearch.co.nz
citation: TBD
              """)
        
        sys.exit(0)
    elif args["task"] == "merge-CH4":
        raise Exception("merge-CH4 is not implemented at this time.\n")
        sys.exit(1)
    else:
        raise Exception("task undefined. I'm not sure how you got to this point...\n")
        sys.exit(1)


if __name__ == "__main__":
    main()
