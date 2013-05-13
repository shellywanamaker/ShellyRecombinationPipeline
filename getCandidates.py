#!/usr/bin/env python

import os
import sys
import signal
import re
import subprocess
from subprocess import call
from subprocess import check_output


class loxData(object):


    # ---- Initiate Object
    def __init__(self):
        path = os.getcwd()
        R1sam = [x for x in os.listdir(path) if "R1" in x and "sam" in x and not "multi" in x]
        R2sam = [x for x in os.listdir(path) if "R2" in x and "sam" in x and not "multi" in x]

        if len(R1sam) > 1 or len(R2sam) >1:
            print("I found too many Bowtie output files in here!")
            print("I only expected to find 1 R1.sam and 1 R2.sam")
            sys.exit(1)

        elif len(R1sam) == 0 or len(R2sam) == 0:
            print("Didn't Find any sam files!")
            print("Make sure your R1 & R2 files end in a .sam")
            sys.exit(1)

        self.R1sam = os.path.realpath(R1sam[0])
        self.R2sam = os.path.realpath(R2sam[0])

        # Check for Chrom Annotations
        # Find the dir where the script is installed
        # Use that to infer where the annotations are kept.
        annotations_folder_name    = "TAIR10_Chrom_Annotations"
        abs_path_to_script         = os.path.realpath(sys.argv[0])
        script_dir                 = os.path.split(abs_path_to_script)[0]

        self.chrom_annotations_dir = os.path.join(script_dir,annotations_folder_name)

        # Is there a Chrom Annotations zip in the dir?
        if not os.path.isdir(self.chrom_annotations_dir) and os.path.isfile(self.chrom_annotations_dir + ".zip"):
            # Unzip annotations folder
            print("Installing TAIR10 Chromosome Annotations in the Script Directory")
            print("This will happen only the first time the script is run on a server")
            print("\tUnzipping...")
            unzip = "unzip %s -d %s" % (self.chrom_annotations_dir + ".zip",script_dir)
            call(unzip,shell=True)
            print
            print("Now Resuming Script.")

        elif not os.path.isdir(self.chrom_annotations_dir) and not os.path.isfile(self.chrom_annotations_dir + ".zip"):
            print("Could not find either the TAIR10 Chrom annotations folder or chrom Annotations Zip")
            print("Please re-clone from GITHUB")
            sys.exit(1)


    # ----  Remove non-aligned; ChrM and ChrC; low quality; and clones
    def slim_and_clean_sam_files(self,no_filter=False,harsh_filter=False):

        if harsh_filter:
            print("Using Harsh Filters")

        R1Basename = os.path.basename(self.R1sam)
        R2Basename = os.path.basename(self.R2sam)

        print("Removing Multi-Aligned and Non-Aligned Reads per Bowtie2 Spec")
        for r in [(self.R1sam,R1Basename),(self.R2sam,R2Basename)]:
            print("\t Removing from %s" % r[0])
            # Remove any ambiguous reads, Keep only aligned and skip headers, sort on column 1
            command = """cat %s | sed '/XS:/d' | awk '{if($3 != "*" && NF > 5)print $0}' | sort -k1,1 > %s.no.multi""" % (r[0],r[1])
            subprocess.call(command,shell=True)

        self.R1sam = R1Basename + ".no.multi"
        self.R2sam = R2Basename + ".no.multi"

        # This step is where the extra SAM information is tossed and the
        # Snipstring filters are.
        print("Slimming Files")
        self.slim_sam_files(self.R1sam,"R1.slim.filtered",no_filter=no_filter,harsh_filter=harsh_filter)
        self.slim_sam_files(self.R2sam,"R2.slim.filtered",no_filter=no_filter,harsh_filter=harsh_filter)

        print("Sorting Slim Files")
        sort_R1 = "sort -k1,1 R1.slim.filtered > R1.slim.filtered.sorted"
        sort_R2 = "sort -k1,1 R2.slim.filtered > R2.slim.filtered.sorted"
        
        call(sort_R1,shell=True)
        call(sort_R2,shell=True)

        print("Joining Slim Files")
        # This is not as complicated as it seems
        # All it is saying is that we want to see every read we pass it
        # If One side didn't align then we are going to see "NA" in the place
        # the columns that should have been seen.
        command = """join -a1 -a2 -j 1 -o 0 1.2 1.3 1.4 1.5 1.6 1.7 2.2 2.3 2.4 2.5 2.6 2.7 -e "NA" %s %s > R1R2.no.genes""" % ("R1.slim.filtered.sorted","R2.slim.filtered.sorted")
        call(command,shell=True)

        print("Removing Clones")
        sort_and_remove_clones = "cat R1R2.no.genes |sort -u -k3,3n -k5,5 -k9,9n -k11,11 > x; mv x R1R2.no.genes.no.clones"
        call(sort_and_remove_clones,shell=True)

        # End
        self.R1R2 = "R1R2.no.genes.no.clones"


    # ---- ---- Helper Function for Slim and Clean Removal of low quality
    # ---- ---- and ChrM/ChrC performed Here
    def slim_sam_files(self,input_file,output_file,no_filter,harsh_filter):
        """
        Removed:
        - all extra sam information that is not needed
        - remove reads that don't start at 1 and are less than 40bps in length
        """
        bad_alignments = set(["mitochondria","chloroplast","ChrC","ChrM"])

        print("\tSlimming and Filtering %s" % input_file)
        with open(input_file,"r") as input_obj:

            with open(output_file,"w") as output:
                
                for i,line in enumerate(input_obj):
                    
                    row = line.strip().split()
                    
                    try:
                        if row[2] in bad_alignments:
                            continue
                    except IndexError:
                        line.strip()
                        continue

                    readID         = row[0]
                    direction      = row[1] 
                    chromosome     = row[2].replace("Chr","").replace("chr","")
                    start_position = row[3]
                    quality        = int(row[4])
                    cigar_string   = row[5]
                    end_position   = self.snipstring2End(cigar_string,start_position)
                    sequence       = row[9]

                    # Find number of Inserts and Deletions
                    deletions = [x for x in cigar_string if x == "D"]
                    inserts   = [x for x in cigar_string if x == "I"]

                    pos = " ".join([readID,chromosome,start_position,end_position,cigar_string,"+",sequence]) + "\n"
                    neg = " ".join([readID,chromosome,start_position,end_position,cigar_string,"-",sequence]) + "\n"

                    # Choose a Filter
                    if no_filter:
                        if direction == "0":
                            to_write = pos

                        elif direction == "16":
                            to_write = neg
                    
                    # Basic Filter where matches must be of length 40 and start within 3bp of the start of the read
                    if not no_filter and not harsh_filter and self.findLength(cigar_string) >= 40 and self.getStartPos(cigar_string) <= 3:
                        if direction == "0":
                            to_write = pos

                        elif direction == "16":
                            to_write = neg
                    
                    # Same as basic except Reads must have no more than One Insert and Deletion
                    if not no_filter and harsh_filter and self.findLength(cigar_string) >= 40 \
                        and self.getStartPos(cigar_string) <= 3 \
                        and len(deletions) <=1 and len(inserts) <= 1\
                        and self.at_least_one_x_length_match(cigar_string,40):

                        if direction == "0":
                            to_write = pos

                        elif direction == "16":
                            to_write = neg

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


    @staticmethod
    def at_least_one_x_length_match(cigar_string,min_length):
        """
        """
        good_cigar = set()
        good_cigar.add("M")
        num = []
        lengths = []


        for char in cigar_string:
            try:
                n = int(char)
                num.append(str(n))
            except ValueError:
                to_length = "".join(num)
                if char in good_cigar:
                    lengths.append(int(to_length))
                num = []

        max_alignment_length = max(lengths)

        if max_alignment_length >= min_length:
            return True
        else:
            return False


    # ---- Get what genes the Tair10 alignments are
    def align2gff(self,debug=False):
        """
        Annotated Chromosomes are read into memory.
        Then the R1R2 files is read one line at a time.
        Using the chromosomal position of the alignment the annotated chromsomes 
        return what genes are on the positive and negative strands.

        NOTE: If the annotations are not found it assumed they are zipped and will be unzipped
        """
        print("Prepping for GFF alignment Step.")

        if debug:
            self.R1R2 = "R1R2.no.genes.no.clones"

        chromosome_files = [os.path.join(self.chrom_annotations_dir,str(x)) for x in range(1,6)]
        chromosomes = {str(x):{} for x in range(1,6)}
        
        for chromosome in chromosome_files:
            print("\tLoading Chromosome %s into Memory" % chromosome)
            chrom = os.path.basename(chromosome)

            with open(chromosome) as in_file:
                for line in in_file:
                    row = line.strip().split("\t")
                    position = row[0]
                    positive = row[1]
                    negative = row[2]

                    chromosomes[chrom][position] = (positive,negative)

        # Go through and match
        print("Aligning %s to the GFF" % self.R1R2)
        with open(self.R1R2,"r") as in_file:
            with open("R1R2.gff.out","w") as out_file:

                for line in in_file:
                    row = line.strip().split()

                    readID           = row[0]
                    
                    geneA_chromosome = row[1]
                    geneA_start      = row[2]
                    geneA_info       = row[1:7]

                    geneB_chromosome = row[7]
                    geneB_start      = row[8]
                    geneB_info       = row[7:]

                    if geneA_chromosome != "NA":
                        geneAPos,geneANeg = chromosomes[geneA_chromosome][geneA_start]
                    else:
                        geneAPos,geneANeg = ("NotAligned","NotAligned")

                    if geneB_chromosome != "NA":
                        geneBPos,geneBNeg = chromosomes[geneB_chromosome][geneB_start]
                    else:
                        geneBPos,geneBNeg = ("NotAligned","NotAligned")

                    to_write = [readID] + geneA_info + [geneAPos,geneANeg] + geneB_info + [geneBPos,geneBNeg]
                    out_file.write(" ".join(to_write) + "\n")


    # ---- Get Candidate Reads by joining R1 and R2 together. Filter and keep genes
    # ---- where R1 != R2
    def getCandidateReads(self):
        print("Getting Candidate Reads")

        print("\tFiltering R1R2: Only keeping reads where GeneA != GeneB")
        self.filterR1R2("R1R2.gff.out")


    # ---- ---- Helper Function for getCandidateReads
    def filterR1R2(self,R1R2,debug=False):
        seen_genes = {}

        with open(R1R2,"r") as all_reads:

            with open("Candidate_Reads.txt","w") as candidate_reads:

                for i,line in enumerate(all_reads):

                    row = line.strip().split()

                    A_line_info = row[:7]
                    B_line_info = row[9:15]

                    geneAPos = row[7]
                    geneANeg = row[8]
                    geneBPos = row[15]
                    geneBNeg = row[16]

                    if geneAPos == "NotAligned" or geneBPos == "NotAligned":
                        continue

                    genes_to_choose = {"geneA":(geneAPos,geneANeg),"geneB":(geneBPos,geneBNeg)}

                    genes = self.chooseGeneAandGeneB(genes_to_choose)

                    geneAinfo = genes["geneA"]
                    geneBinfo = genes["geneB"]

                    if debug:
                        print line.strip(),geneAinfo,geneBinfo

                    geneA = geneAinfo.split(":")[0].split(".")[0]
                    geneB = geneBinfo.split(":")[0].split(".")[0]

                    # Count the amount of times we see genes
                    seen_genes[geneA] = seen_genes.get(geneA,0) + 1
                    seen_genes[geneB] = seen_genes.get(geneB,0) + 1

                    if geneA != geneB:
                        candidate_reads.write(" ".join([" ".join(A_line_info),geneAinfo," ".join(B_line_info),geneBinfo + "\n"]))

        # Print out Seen_Genes
        seen_genes_keys = seen_genes.keys()[:]
        seen_genes_keys.sort()

        print("Writing seen genes out to GENE_COUNTS.csv")
        with open("GENE_COUNTS.csv","w") as out_file:
            for gene in seen_genes_keys:
                out_file.write("%s,%s\n" % (gene,seen_genes[gene]))


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

    
    # ---- Removes temporary files in current working directory
    def cleanUp(self):
        """
        This will clean up non major intermediates and leave only the files after
        the largest steps.
        """

        print("Cleaning up Temp files")
        clean = "rm R1.slim* R2.slim* R1R2.no.genes* *.no.multi"
        call(clean,shell=True)


if __name__== "__main__":

    # ---- Play Nice with Unix
    signal.signal(signal.SIGPIPE,signal.SIG_DFL)

    # ---- Script start
    loxData = loxData()

    # ---- Run Methods
    # loxData.slim_and_clean_sam_files(no_filter=False,harsh_filter=True)
    # loxData.align2gff(debug=True)
    loxData.getCandidateReads()
    # loxData.cleanUp()
