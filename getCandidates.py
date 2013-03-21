#!/usr/bin/env python

import os
import sys
import signal
from subprocess import call
from subprocess import check_output


class loxData(object):

    def __init__(self,path):

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

    def bowtie(self,options="--local -p 8",indexes="/home/shelly/bin/bowtie2/INDEXES/tair10"):
        """
        As long as Align output is Sam the script will still work.
        """

        # Bowtie R1
        print("Bowtie-ing R1 reads")
        commandR1 = " ".join(["bowtie2",options,indexes,",".join(self.R1),"1> bowtie.R1.out.sam 2> bowtie.R1.stats"])
        call(commandR1,shell=True)

        # Bowtie R2
        print("Bowtie-ing R2 reads")
        commandR2 = " ".join(["bowtie2",options,indexes,",".join(self.R2),"1> bowtie.R2.out.sam 2> bowtie.R2.stats"])
        call(commandR2,shell=True)

        # End
        self.R1sam = "bowtie.R1.out.sam"
        self.R2sam = "bowtie.R2.out.sam"

    def slim_and_clean_sam_files(self):

        print("Slimming Files")
        self.slim_sam_files(self.R1sam,"R1.slim")
        self.slim_sam_files(self.R2sam,"R2.slim")

        print("Removing Clones from Slim Files")
        for r in ["R1.slim","R2.slim"]:
            print("\tRemoving Clones from %s" % r)

            command = "cat %s | sort -u -k2,2n -k3,3n -k4,4n > %s.clean" % (r,r)
            call(command,shell=True)

        # End
        self.R1slim = "R1.slim.clean"
        self.R2slim = "R2.slim.clean"

    def slim_sam_files(self,input_file,output_file):

        with open(input_file,"r") as input_obj:

            with open(output_file,"w") as output:
                
                for i,line in enumerate(input_obj):
                    
                    row = line.strip().split()
                    
                    if i < 9:
                        continue
                    elif row[1] == "4":
                        continue
                    elif row[2] in ["mitochondria","chloroplast","ChrC","ChrM"]:
                        continue

                    readID         = row[0]
                    direction      = row[1] 
                    chromosome     = row[2].replace("Chr","").replace("chr","")
                    start_position = row[3]
                    end_position   = self.snipstring2End(row[5],start_position)
                    sequence       = row[9]

                    if direction == "0":
                        output.write(" ".join([readID,chromosome,start_position,end_position,"+",sequence]) + "\n")

                    elif direction == "16":
                        output.write(" ".join([readID,chromosome,start_position,end_position,"-",sequence]) + "\n")

    def snipstring2End(self,snip_string,start_position):
        """
        Given a SAM formated snip string and a start position this method will
        calculate the length the end position of the read
        """

        num = 0
        temp = []

        for char in snip_string:

            temp.append(char)

            if char == "S":
                temp = []
                continue

            elif char == "M":
                num += int("".join(temp[:-1]))
                temp = []

            elif char == "D":
                num += int("".join(temp[:-1]))
                temp = []

            elif char == "I":
                num -= int("".join(temp[:-1]))
                temp = []

        return str(int(start_position) + int(num))

    def align2gff(self):
        """
        This may change at a certain point. Be sure that this method can be ripped out.
        """
        print("Prepping Directory for Alignment/Gene Matching")
        make_gff_output_dir = "mkdir anno"
        call(make_gff_output_dir,shell=True)

        move_alignment_files = "mkdir alignments;mv *.sam alignments/;mv *.stats alignments/"
        call(move_alignment_files,shell=True)

        # Counts

        print("Running Alignment -> Gene matching")
        run_gff4joe = "/home/feeney/bin/gff4joe/gff4joe.pl %s" % (os.getcwd())
        call(run_gff4joe,shell=True)

        # Bring Back out the files to CWD. Remove anno directory
        print("Removing Temp GFF Folders")
        restructure_folder = "mv %s/anno/*.clean ./; rm -r anno" % (os.getcwd())
        call(restructure_folder,shell=True)

        # The end of this method only needs to provide R1 and R2 .slim.clean in the CWD

    def getCandidateReads(self):
        print("Getting Candidate Reads")

        print("\tSorting R1 by ReadID")
        sort_R1 = "sort -k1,1 R1.slim.clean > x; mv x R1.slim.clean"
        call(sort_R1,shell=True)

        print("\tSorting R2 by ReadID")
        sort_R2 = "sort -k1,1 R2.slim.clean > x; mv x R2.slim.clean"
        call(sort_R2,shell=True)

        print("\tJoining R1 and R2")
        join_R1_and_R2 = "join -j1 R1.slim.clean R2.slim.clean > R1R2"
        call(join_R1_and_R2,shell = True)

        print("\tFiltering R1R2: Only keepin reads where GeneA != GeneB")
        self.filterR1R2("R1R2")

        print("Cleaning up Temp files")

    def filterR1R2(self,R1R2):

        with open(R1R2,"r") as all_reads:

            with open("Candidate_Reads.txt","w") as candidate_reads:

                for i,line in enumerate(all_reads):

                    row = line.strip().split()
                    # print line

                    A_line_info = row[:6]
                    B_line_info = row[8:13]

                    geneAPos = row[6]
                    geneANeg = row[7]
                    geneBPos = row[13]
                    geneBNeg = row[14]

                    genes_to_choose = {"geneA":(geneAPos,geneANeg),"geneB":(geneBPos,geneBNeg)}

                    genes = self.chooseGeneAandGeneB(genes_to_choose)

                    geneAinfo = genes["geneA"]
                    geneBinfo = genes["geneB"]

                    geneA = geneAinfo.split(":")[0].split(".")[0]
                    geneB = geneBinfo.split(":")[0].split(".")[0]

                    if geneA != geneB:
                        candidate_reads.write(" ".join([" ".join(A_line_info),geneAinfo," ".join(B_line_info),geneBinfo + "\n"]))

    def chooseGeneAandGeneB(self,geneDict):
        for gene in geneDict:

            pos = geneDict[gene][0]
            neg = geneDict[gene][1]

            if pos != "Intergenic":

                if pos.split(":")[1] == "Three":
                    pos_val = (5,pos)

                elif pos.split(":")[1] == "CDS":
                    pos_val = (4,pos)

                elif pos.split(":")[1] == "Promoter":
                    pos_val = (3,pos)

                elif pos.split(":")[1] == "Post":
                    pos_val = (2,pos)
                
                elif pos.split(":")[1] == "pseudo_start":
                    pos_val = (2,pos)

                elif pos.split(":")[1] == "Five":
                    pos_val = (2,pos)

                elif pos.split(":")[1] == "pseudogenic_transcript":
                    pos_val = (2,pos)

                elif pos.split(":")[1] == "pseudo_end":
                    pos_val = (2,pos)

            else:
                pos_val = (1,pos)

            if neg != "Intergenic":

                if neg.split(":")[1] == "Three":
                    neg_val = (5,neg)

                elif neg.split(":")[1] == "CDS":
                    neg_val = (4,neg)

                elif neg.split(":")[1] == "Promoter":
                    neg_val = (3,neg)

                elif neg.split(":")[1] == "Post":
                    neg_val = (2,neg)

                elif neg.split(":")[1] == "pseudo_start":
                    neg_val = (2,neg)
                
                elif neg.split(":")[1] == "Five":
                    neg_val = (2,neg)
                
                elif neg.split(":")[1] == "pseudogenic_transcript":
                    neg_val = (2,neg)

                elif neg.split(":")[1] == "pseudo_end":
                    neg_val = (2,neg)


            else:
                neg_val = (1,neg)

            
            gene_tuple = max(pos_val,neg_val)
            geneDict[gene] = gene_tuple[1]

        return geneDict


def startupChecks():

    # Check for Bowtie2
    if not check_PATH_for_program("bowtie2"):
        print("\nBowtie2 is not in your $PATH. Can you modify your .bash_profile or .bashrc, please?\n")
        sys.exit(1)

    # Check for a .CSV file for a Candidate (if any) reads
    csvs = [x for x in os.listdir(os.getcwd()) if ".csv" in x]

    try:
        return sys.argv[1]
    except IndexError:
        return os.getcwd()

def check_PATH_for_program(f):
    """
    Check the unix $PATH environment for specified program
    """

    path = os.environ["PATH"].split(":")

    for p in path:

        if os.path.isfile(p + "/" + f):
            return True

    else:
        return False



if __name__== "__main__":

    # Play Nice with Unix
    signal.signal(signal.SIGPIPE,signal.SIG_DFL)

    # Script start
    path    = startupChecks()

    loxData = loxData(path)

    loxData.bowtie()
    loxData.slim_and_clean_sam_files()
    loxData.align2gff()
    loxData.getCandidateReads()
