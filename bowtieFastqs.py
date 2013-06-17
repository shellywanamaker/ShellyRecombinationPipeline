#!/usr/bin/env python

import os
import sys
import signal
from subprocess import call


class bowtieFastqs(object):


    def __init__(self):
        
        try:
            path = sys.argv[1]
        except IndexError:
            print("No Folder given. Assuming Fastq's are in current directory.")
            path = os.getcwd()

        # Make a list of all the Absolute paths to the fastq's
        # Assume Run is pair end. Break run in to pair end tuples
        fastqs = [path + "/" + x for x in os.listdir(path) if ".fastq" in x]
        R1     = [x for x in fastqs if "_R1_" in x]
        R2     = [x for x in fastqs if "_R2_" in x]

        R1.sort()
        R2.sort()

        if len(fastqs) < 1:
            print("Could not find .fastq files in either the provided directory or in the current directory")
            sys.exit(1)

        elif len(R1) != len(R2):
            print("The number of R1's and R2's don't match")
            sys.exit(1)

        # End
        self.R1 = R1[:]
        self.R2 = R2[:]
        self.path = path


    def bowtie(self,options="--local -p 3",indexes_folder="/mnt/sculpin/data5/Shelly/bin/bowtie2/INDEXES/tair10.cDNA",genome_basename="tair10.cDNA"):
        """
        As long as Align output is Sam the script will still work.
        """
        if not os.path.isdir(indexes_folder):
            print("Could not find your INDEXES Folder: %s" % indexes_folder)

            while True:

                indexes_folder = raw_input("What is the path (abs or relative) to the Bowtie2 INDEXES: ")
                indexes_folder = os.path.abspath(os.path.expanduser(indexes_folder))

                if os.path.isdir(indexes_folder) and\
                        len([x for x in os.listdir(indexes_folder) if genome_basename in x]) > 0:

                    print("Looks like that will work!")
                    break

                elif os.path.isdir(indexes_folder):
                    print("I couldn't find a genome with a basename %s in %s" %(genome_basename,indexes_folder))
                    print("Try another folder")

                else:
                    print("Looks like that folder doesn't exist!")


        # Bowtie to Yeast and Tair10
        for genome in [genome_basename]:
            # More specific for options for each genome
            if genome == "yeast":
                options += " "

            # Bowtie R1
            indexes = os.path.join(indexes_folder,genome)

            print("Bowtie-ing R1 reads to %s" % genome)
            commandR1 = " ".join(["bowtie2",options,indexes,",".join(self.R1),"1> bowtie.R1.%s.sam 2> bowtie.R1.%s.stats" % (genome,genome)])
            call(commandR1,shell=True)

            # Bowtie R2
            print("Bowtie-ing R2 reads %s" % genome)
            commandR2 = " ".join(["bowtie2",options,indexes,",".join(self.R2),"1> bowtie.R2.%s.sam 2> bowtie.R2.%s.stats" % (genome,genome)])
            call(commandR2,shell=True)

        # # Loading Bowtied Yeast ReadIDs into memory
        # yeast_bowtie_output = [x for x in os.listdir(os.getcwd()) if "yeast" in x and "sam" in x]
        # readIDs_to_remove   = set()

        # for f in yeast_bowtie_output:
        #     print("\tLoading %f into Memory" % f)
        #     with open(f,"r") as input_file:
        #         for line in input_file:
        #             row = line.strip().split()

        #             readID    = row[0]
        #             alignment = row[2]

        #             if alignment != "*":
        #                 readIDs_to_remove.add(readID)

        # # Using these ReadID's parse the Tair10 sam files and remove readIDs
        # # that also bowtied to Yeast
        # print("Removing Yeast ReadIDs from Tair10 sam files")

        # tair_bowtie_output = [x for x in os.listdir(os.getcwd()) if ".sam" in x and "tair"]

        # for tair in tair_bowtie_output:
        #     tair = os.path.join("../known_positives/alignments/",tair)

        #     if "R1" in tair:
        #         output_file = open("bowtie.R1.no.yeast.sam","w")
        #     elif "R2" in tair:
        #         output_file = open("remove.R2.no.yeast.sam","w")

        #     with open(tair,"r") as t_file:
        #         for line in t_file:
        #             row = line.strip().split()

        #             readID    = row[0]
        #             alignment = row[2]

        #             if readID not in yeast_readIDs and alignment != "*":
        #                 output_file.write(line)

        #     output_file.close()


def check_PATH_for_program(f):
        """
        Check the unix $PATH environment for specified program
        """

        path = os.environ["PATH"].split(":")

        for p in path:

            if os.path.isfile(os.path.join(p,f)):
                return True

        return False



if __name__=="__main__":

    genome_basename = "tair10.cDNA"
    options         = "--local -p 3"
    indexes_folder  = "/mnt/sculpin/data5/Shelly/bin/bowtie2/INDEXES/tair10.cDNA" 

    # ---- Check for Bowtie2
    if not check_PATH_for_program("bowtie2"):
        print("\nBowtie2 is not in your $PATH. Can you modify your .bash_profile or .bashrc, please?\n")
        sys.exit(1)

    # ---- Run Script
    b = bowtieFastqs()
    b.bowtie(options=options,indexes_folder=indexes_folder,genome_basename=genome_basename)
