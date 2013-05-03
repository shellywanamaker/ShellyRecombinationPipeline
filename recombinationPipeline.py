#!/usr/bin/env python

from filterCandidates import *
from getCandidates import *
from bowtieFastqs import *


if __name__=="__main__":

    # ---- Play Nice with Unix
    signal.signal(signal.SIGPIPE,signal.SIG_DFL)

    # ---- Prepping for Run
    csv_files = [x for x in os.listdir(os.getcwd()) if ".csv" in x]

    if len(csv_files) == 0:
        print("\nCouldn't Find a .csv to work against!\n")
        sys.exit(1)

    elif len(csv_files) > 1:
        print("Found too many csv files in directory. Only expected one!")
        sys.exit(1)

    else:
        csv_file = csv_files[0]
        print("Using %s as the Gene Pairs csv" % (csv_file))

    # ---- Script start
    # ---- Instantiate Classes
    bowtieFastqs     = bowtieFastqs()
    loxData          = loxData(path)
    filterCandidates = filterData(csv_file)

    # ---- bowtieFastqs
    bowtieFastqs.bowtie()

    # ---- getCandidates Module
    loxData.slim_and_clean_sam_files(no_filter=False,harsh_filter=False)
    loxData.align2gff()
    loxData.getCandidateReads()

    # ---- filter Candidates
    filterCandidates.compareCandidateReads2Predicted()
