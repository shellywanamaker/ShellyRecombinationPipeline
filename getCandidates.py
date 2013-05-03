#!/usr/bin/env python

import os
import sys
import signal
from subprocess import call
from subprocess import check_output


class loxData(object):

    # ---- Initiate Object
    def __init__(self):
        R1sam = [x for x in os.listdir(path) if "R1" in x and "sam" in x]
        R2sam = [x for x in os.listdir(path) if "R2" in x and "sam" in x]

        if len(R1sam) > 1 or len(R2sam) >1:
            print("I found too many Bowtie output files in here!")
            print("I only expected to find 1 R1.sam and 1 R2.sam")
            sys.exit(1)

        self.R1sam = os.path.realpath(R1sam[0])
        self.R2sam = os.path.realpath(R2sam[0])

        print self.R1sam,self.R2sam

    # ----  Remove non-aligned; ChrM and ChrC; low quality; and clones
    def slim_and_clean_sam_files(self,no_filter=False,harsh_filter=False):

        print("Slimming Files")
        self.slim_sam_files(self.R1sam,"R1.slim",no_filter=no_filter,harsh_filter=harsh_filter)
        self.slim_sam_files(self.R2sam,"R2.slim",no_filter=no_filter,harsh_filter=harsh_filter)

        print("Removing Clones from Slim Files")
        for r in ["R1.slim","R2.slim"]:
            print("\tRemoving Clones from %s" % r)

            command = "cat %s | sort -u -k2,2n -k3,3n -k4,4n > %s.clean" % (r,r)
            call(command,shell=True)

        # End
        self.R1slim = "R1.slim.clean"
        self.R2slim = "R2.slim.clean"

    # ---- ---- Helper Function for Slim and Clean Removal of low quality
    # ---- ---- and ChrM/ChrC performed Here
    def slim_sam_files(self,input_file,output_file,no_filter,hash_filter):
        """
        Removed:
        - all extra sam information that is not needed
        - remove reads that don't start at 1 and are less than 40bps in length
        """
        print("\tSlimming %s" % input_file)
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
                    quality        = int(row[4])
                    cigar_string   = row[5]
                    end_position   = self.snipstring2End(cigar_string,start_position)
                    sequence       = row[9]

                    # Choose a Filter
                    if no_filter:
                        if direction == "0":
                            to_write = " ".join([readID,chromosome,start_position,end_position,"+",sequence]) + "\n"

                        elif direction == "16":
                            to_write = " ".join([readID,chromosome,start_position,end_position,"-",sequence]) + "\n"
                    
                    # Basic Filter where matches must be of length 40 and start within 3bp of the start of the read
                    if not no_filter and not harsh_filter and self.findLength(cigar_string) >= 40 and self.getStartPos(cigar_string) <= 3:
                        if direction == "0":
                            to_write = " ".join([readID,chromosome,start_position,end_position,"+",sequence]) + "\n"

                        elif direction == "16":
                            to_write = " ".join([readID,chromosome,start_position,end_position,"-",sequence]) + "\n"
                    
                    # Same as basic except Reads must have no more than One Insert and Deletion
                    if not no_filter and hash_filter and self.findLength(cigar_string) >= 40 \
                        and self.getStartPos(cigar_string) <= 3 \
                        and len(deletions) <=1 and len(inserts) <= 1:

                        if direction == "0":
                            to_write = " ".join([readID,chromosome,start_position,end_position,"+",sequence]) + "\n"

                        elif direction == "16":
                            to_write = " ".join([readID,chromosome,start_position,end_position,"-",sequence]) + "\n"

                    else:
                        continue
                    

                    output.write(to_write)


    # ---- ---- Helper Functions for Slim Sam Files
    def snipstring2End(self,cigar_string,start_position):
        """
        Given a SAM formated snip string and a start position this method will
        calculate the length the end position of the read
        """

        return str(int(start_position) + self.findLength(cigar_string) - 1)


    @staticmethod
    def findLength(cigar_string):
        good_cigar = set()
        good_cigar.add("M")
        good_cigar.add("I")

        num = []
        length = 0
        for char in cigar_string:
            try:
                n = int(char)
                num.append(str(n))
            except ValueError:
                to_length = "".join(num)
                if char in good_cigar:
                    length += int(to_length)
                num = []

        return int(length)


    @staticmethod
    def getStartPos(cigar_string):
        """
        given a sam cigar string return the position of the First aligned Base pair
        in the alignment
        """
        temp =[]

        for char in cigar_string:
            try:
                n = int(char)
                temp.append(char)

            except ValueError:
                if char == "M":
                    return 1
                else:
                    return int("".join(temp)) + 1


    # ---- Get what genes the Tair10 alignments are
    def align2gff(self):
        """
        This may change at a certain point. Be sure that this method can be ripped out.
        """
        # New implementation
        # Read Chroms into memory

        # Find the dir where the script is installed
        abs_path_to_script = os.path.realpath(sys.argv[0])
        script_dir = os.path.split(abs_path_to_script)[0]
        chrom_annotations_dir = os.path.join(script_dir,"TAIR10_Chrom_Annotations")

        chromosome_files = [os.path.join(chrom_annotations_dir,str(x)) for x in range(1,6)]
        chromosomes = {str(x):{} for x in range(1,6)}
        
        for chromosome in chromosome_files:
            print("Loading Chromosome %s into Memory" % chromosome)
            chrom = os.path.basename(chromosome)

            with open(chromosome) as in_file:
                for line in in_file:
                    row = line.strip().split("\t")
                    position = row[0]
                    positive = row[1]
                    negative = row[2]

                    chromosomes[chrom][position] = (positive,negative)

        # Go through and match
        for tup in [("R1.slim.clean","R1.gff.out"),("R2.slim.clean","R2.gff.out")]:
            print("Aligning %s to the GFF" % tup[0])
            with open(tup[0],"r") as in_file:
                with open(tup[1],"w") as out_file:

                    for line in in_file:
                        row = line.strip().split()

                        aligned_chromosome = row[1]
                        aligned_start = row[2]

                        genes = chromosomes[aligned_chromosome][aligned_start]

                        to_write = [line.strip(),genes[0],genes[1]]
                        out_file.write(" ".join(to_write) + "\n")

        # The end of this method only needs to provide R1 and R2 .slim.clean in the CWD

    # ---- Get Candidate Reads by joining R1 and R2 together. Filter and keep genes
    # ---- where R1 != R2
    def getCandidateReads(self):
        print("Getting Candidate Reads")

        print("\tSorting R1 by ReadID")
        sort_R1 = "sort -k1,1 R1.gff.out > x; mv x R1.gff.out"
        call(sort_R1,shell=True)

        print("\tSorting R2 by ReadID")
        sort_R2 = "sort -k1,1 R2.gff.out > x; mv x R2.gff.out"
        call(sort_R2,shell=True)

        print("\tJoining R1 and R2")
        join_R1_and_R2 = "join -j1 R1.gff.out R2.gff.out > R1R2"
        call(join_R1_and_R2,shell = True)

        print("\tFiltering R1R2: Only keepin reads where GeneA != GeneB")
        self.filterR1R2("R1R2")

        print("Cleaning up Temp files")

    # ---- ---- Helper Function for getCandidateReads
    def filterR1R2(self,R1R2,debug=False):

        with open(R1R2,"r") as all_reads:

            with open("Candidate_Reads.txt","w") as candidate_reads:

                for i,line in enumerate(all_reads):

                    row = line.strip().split()

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

                    if debug:
                        print line.strip(),geneAinfo,geneBinfo

                    geneA = geneAinfo.split(":")[0].split(".")[0]
                    geneB = geneBinfo.split(":")[0].split(".")[0]

                    if geneA != geneB:
                        candidate_reads.write(" ".join([" ".join(A_line_info),geneAinfo," ".join(B_line_info),geneBinfo + "\n"]))

    # ---- ---- Helper Function for filterR1R2. Decides how to rank the two strands.
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


if __name__== "__main__":

    # ---- Play Nice with Unix
    signal.signal(signal.SIGPIPE,signal.SIG_DFL)

    # ---- Script start
    loxData = loxData()

    # ---- Run Methods
    loxData.slim_and_clean_sam_files(no_filter=False,harsh_filter=False)
    loxData.align2gff()
    loxData.getCandidateReads()
