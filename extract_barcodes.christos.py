"""

Updated Nov 2021 for Christos' UTR MPRA

"""
import argparse
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import numpy as np
import re
from itertools import islice
from collections import defaultdict
from time import gmtime, strftime

##########################
###### VARIABLES #########
##########################


# this updated script assumes single end reads only
UPSTREAM_CONSTANT_SEQ="GATATTTTATTGCGGCCAGC"
VERY_DOWNSTREAM_CONSTANT_SEQ="GCGATCGCCTAGAATTACTG"

LENGTH_OF_DOWNSTREAM_TO_CHECK = 10
BARCODE_LENGTH = 8

##########################
###### FUNCTIONS #########
##########################

def parseIndex(index, header):
    """
    read in an index dataframe using pandas

    arguments
    ---------
    index: str, filename or path of index_df
    header: str, filename or path of the header (optional)

    returns
    -------
    index_df: pandas dataframe
    """

    # case where there is already a header (don't pass a separate header file through)
    if header == None:
        index_df = pd.read_table(index, sep="\t")

    # case where the header is in a separate file
    else:
        index_df = pd.read_table(index, sep="\t", header=None)
        header = pd.read_table(header, sep="\t", header=None)
        index_df.columns = list(header[0])

    return index_df

def createReverseComplement(seq_str):
    """
    one-liner to vectorize reverse complementing a string, to pass to pandas apply

    arguments
    ---------
    seq_str: string to be reverse complemented

    returns
    -------
    rev_comp: string that has been reverse complemented
    """

    rev_comp = str(Seq(seq_str).reverse_complement())
    return rev_comp

def expandIndexDetails(index_df):
    """
    function to add details to index that are necessary for parsing the reads

    arguments
    ---------
    index_df: pandas dataframe with index_df

    returns
    -------
    sequences_to_check_for: pandas dataframe with barcode, barcode rev comp, element, element rev comp, etc
    """

    # create dataframe with sequences to check for. 
    sequences_to_check_for = pd.DataFrame(columns = ["barcode", "element"])
    sequences_to_check_for["barcode"] = index_df["barcode"].str.upper()
    sequences_to_check_for["barcode_rev_comp"] = sequences_to_check_for["barcode"].apply(createReverseComplement)
    try:
        sequences_to_check_for["element"] = index_df["element"].str.upper()
        sequences_to_check_for["element_rev_comp"] = sequences_to_check_for["element"].apply(createReverseComplement)
    except:
        sequences_to_check_for["element"] = ""
        sequences_to_check_for["element_rev_comp"] = ""

    # subset the element to the first X bases you want to check
    sequences_to_check_for["element_rev_comp_sub"] = sequences_to_check_for["element_rev_comp"].str[0:LENGTH_OF_DOWNSTREAM_TO_CHECK]
    sequences_to_check_for["element_sub"] = sequences_to_check_for["element"].str[0:LENGTH_OF_DOWNSTREAM_TO_CHECK]

    # add column for very downstream constant region (which we actually don't check but keep it there for now)
    sequences_to_check_for["constant_downstream"] = VERY_DOWNSTREAM_CONSTANT_SEQ

    return sequences_to_check_for

def createIndexDict(sequences_to_check_for):
    """
    function to create the dictionary needed to search the reads for barcodes and element_rev_comp_sub

    arguments
    ---------
    sequences_to_check_for: pandas dataframe made via function above

    returns
    -------
    sequences_to_check_for_dict: dictionary keyed by barcode, with elements to check for as values.
        0th value is always what to check for in r1 (rev comp'd read) and 1st is what to check in r2 (not rev comp'd)

    """

    seqs_to_check_for_dict = dict(zip(sequences_to_check_for["barcode_rev_comp"], zip(sequences_to_check_for["element_rev_comp_sub"], sequences_to_check_for["constant_downstream"])))
    return seqs_to_check_for_dict

def extractBarcodes(r1_fastq, index_df):

    # first make the dictionary of sequences to check for
    sequences_to_check_for = expandIndexDetails(index_df)
    n_barcodes = len(sequences_to_check_for)
    print("checking for %s barcodes..." % (n_barcodes))
    seqs_to_check_for_dict = createIndexDict(sequences_to_check_for)

    # store the barcodes to look for as a set
    rev_comp_barcode_set = set(seqs_to_check_for_dict.keys())
    barcode_set = set(sequences_to_check_for["barcode"])

    # initialize dicts to store barcode counts
    # single_lenient: only match the constant region + barcode
    # single_strict: match above + additional 10 bp after constant region
    single_lenient = dict.fromkeys(barcode_set, 0)
    single_strict = dict.fromkeys(barcode_set, 0)

    # initialize a variable to count reads
    tot_reads = 0

    # note the time 
    start_time = strftime("%Y-%m-%d %H:%M:%S")
    print("parsing fastq file, starting at %s" % start_time)

    # open gzip'd fastq and iterate through it every 4 lines
    # make sure to deal with single end and paired end differently
    if r1_fastq.endswith(".gz"):
        f = gzip.open(r1_fastq, "r")
    else:
        f = open(r1_fastq, "r")
    seq_iterator = islice(f, 1, None, 4)

    for item in seq_iterator:

        # the most efficient way i could conceive of doing this is to iterate directly through the items in the iterator
        # if it's a single-end iterator, each item will only have 1 read
        r1 = item.strip("\n")

        # print "item: %s" % (item)
        # print ""
        # print "r1: %s" % (r1)
        # print "r2: %s" % (r2)
        # different strategies depending on the sequencing type/cloning step combination

        if UPSTREAM_CONSTANT_SEQ in r1:
            constant_start = r1.find(UPSTREAM_CONSTANT_SEQ)
            barcode_start = constant_start + len(UPSTREAM_CONSTANT_SEQ)
            barcode_end = barcode_start + BARCODE_LENGTH
            found_barcode = r1[barcode_start:barcode_end]

            if found_barcode in rev_comp_barcode_set:
                # at this step, convert it back from rev comp so we can have DFs of actual barcodes
                # and add 1 to lenient dictionary
                actual_barcode = str(Seq(found_barcode).reverse_complement())
                single_lenient[actual_barcode] += 1

                # continue with strict checks
                expected_r1 = seqs_to_check_for_dict[found_barcode][0]
                r1_elem_start = barcode_end
                r1_elem_end = r1_elem_start + LENGTH_OF_DOWNSTREAM_TO_CHECK
                found_r1 = r1[r1_elem_start:r1_elem_end]

                if expected_r1 == found_r1:
                    single_strict[actual_barcode] += 1

        tot_reads += 1


    # record end time
    end_time = strftime("%Y-%m-%d %H:%M:%S")
    print("done parsing fastq file (%s total reads), ended at %s" % (tot_reads, end_time))

    # turn everything back into df
    single_lenient_df = pd.DataFrame.from_dict(single_lenient, orient="index").reset_index()
    single_strict_df = pd.DataFrame.from_dict(single_strict, orient="index").reset_index()

    print("================ SINGLE-END READS, LENIENT ====================")
    try:
        single_lenient_df.columns = ["barcode", "count"]
        single_lenient_df_sum = single_lenient_df.sum(axis=0, numeric_only=True).iloc[0]
        single_lenient_barcodes = len(single_lenient_df[single_lenient_df["count"] > 0])
        print("**** Found %s unique barcodes with perfect matches to constant downstream element %s ***" % (single_lenient_barcodes, UPSTREAM_CONSTANT_SEQ))
        print("**** Those unique barcodes comprised %s total reads ***" % (single_lenient_df_sum))
        print("**** which corresponds to %s percent of total reads ****" % (float(single_lenient_df_sum)/tot_reads*100))
        print("")

    except ValueError:
        print("**** Found ZERO unique barcodes with perfect matches to constant downstream element %s ***" % (UPSTREAM_CONSTANT_SEQ))
        print("")

    print("================ SINGLE-END READS, STRICT ====================")
    try:
        single_strict_df.columns = ["barcode", "count"]
        single_strict_df_sum = single_strict_df.sum(axis=0, numeric_only=True).iloc[0]
        single_strict_barcodes = len(single_strict_df[single_strict_df["count"] > 0])
        print("**** Found %s unique barcodes with perfect matches to constant downstream element %s and %s extra bases ***" % (single_strict_barcodes, UPSTREAM_CONSTANT_SEQ, LENGTH_OF_DOWNSTREAM_TO_CHECK))
        print("**** Those unique barcodes comprised %s total reads ***" % (single_strict_df_sum))
        print("**** which corresponds to %s percent of total reads ****" % (float(single_strict_df_sum)/tot_reads*100))
        print("")

    except ValueError:
        print("**** Found ZERO unique barcodes with perfect matches to constant downstream element %s and %s extra bases ***" % (UPSTREAM_CONSTANT_SEQ, LENGTH_OF_DOWNSTREAM_TO_CHECK))
        print("")

    return single_lenient_df, single_strict_df

##########################
######    MAIN   #########
##########################

def main(strInput=None):

    # Define arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--r1_fastq", type=str, required=True, help="FASTQ of interest (r1)")
    parser.add_argument("-i", "--index", type=str, required=True, help="index file")
    parser.add_argument("-e", "--header", type=str, required=False, default=None, help="header file for index")
    parser.add_argument("-o", "--output_file", type=str, required=True, help="output file to store extracted seqs")

    # Collect arguments
    args = parser.parse_args(strInput.split()) if strInput else parser.parse_args()
    r1_fastq = args.r1_fastq
    index = args.index
    header = args.header
    output_lenient_file = args.output_file + "LENIENT_BARCODES.txt"
    output_strict_file = args.output_file + "STRICT_BARCODES.txt"

    index_df = parseIndex(index, header)

    final_lenient_barcodes, final_strict_barcodes = extractBarcodes(r1_fastq, index_df)

    # Write dataframe
    print("wrote single-end lenient barcodes to file %s" % (output_lenient_file))
    final_lenient_barcodes.to_csv(output_lenient_file, sep="\t", index=False)

    print("wrote single-end strict barcodes to file %s" % (output_strict_file))
    final_strict_barcodes.to_csv(output_strict_file, sep="\t", index=False)



if __name__ == "__main__":
    main()

