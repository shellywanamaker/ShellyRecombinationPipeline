#!/usr/bin/env python

from filterCandidates import *
from getCandidates import *

if __name__=="__main__":

    # ---- Play Nice with Unix
    signal.signal(signal.SIGPIPE,signal.SIG_DFL)

    # ---- Prepping for Run
    csv_files = [x for x in os.listdir(os.getcwd()) if ".csv" in x]

    if len(csv_files) == 0:
        print("\nCouldn't Find a .csv to Work against!\n")
        sys.exit(1)

    elif len(csv_files) > 1:
        print("Found too many csv files in directory. Only expected one!")
        sys.exit(1)

    else:
        csv_file = csv_files[0]
        print("Using %s as the Gene Pairs csv" % (csv_file))

    # ---- Script start
    path    = startupChecks()
    loxData = loxData(path)
    f       = filterData(csv_file)

    loxData.bowtie()
    loxData.slim_and_clean_sam_files()
    loxData.align2gff()
    loxData.getCandidateReads()

    f.compareCandidateReads2Predicted()
