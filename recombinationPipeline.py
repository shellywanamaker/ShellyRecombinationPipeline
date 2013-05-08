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

    # ---- Options for Bowtie
    genome_basename = "tair10"
    options         = "--local -p 3"
    indexes_folder  = "/home/shelly/bin/bowtie2/INDEXES/"


    # ---- Instantiate Classes
    bowtieFastqs     = bowtieFastqs()
    loxData          = loxData()
    filterCandidates = filterData(csv_file)

    # ---- bowtieFastqs.py
    bowtieFastqs.bowtie(options=options,indexes_folder=indexes_folder,genome_basename=genome_basename)

    # ---- getCandidates.py
    loxData.slim_and_clean_sam_files(no_filter=False,harsh_filter=True)
    loxData.align2gff()
    loxData.getCandidateReads()

    # ---- filterCandidates.py
    # NOTE: If modify accession numbers is False the script will still check the
    #       the closest accession number matches. It just won't modify genes to the
    #       closest one.
    filterCandidates.compareCandidateReads2Predicted(modify_accession_numbers=False)
