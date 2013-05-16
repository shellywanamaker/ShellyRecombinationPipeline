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
    def __init__(self,genome="tair10"):
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

        self.R1sam  = os.path.realpath(R1sam[0])
        self.R2sam  = os.path.realpath(R2sam[0])
        self.genome = genome

        # Check for Chrom Annotations
        # Find the dir where the script is installed
        # Use that to infer where the annotations are kept.
        annotations_folder_name    = "Chromosome_Annotations"
        abs_path_to_script         = os.path.realpath(sys.argv[0])
        script_dir                 = os.path.split(abs_path_to_script)[0]
        
        # Current Implementation assumes that the Chromosome Annotations will
        # always exist but the tair10
        self.chrom_annotations_dir  = os.path.join(script_dir,annotations_folder_name)
        self.genome_annotations_dir = os.path.join(self.chrom_annotations_dir,genome)

        if not os.path.isdir(self.genome_annotations_dir):
            print("Can't Find the chromosome annotations for %s in %s" % (self.genome,self.chrom_annotations_dir))

            gff = [x for x in os.listdir(self.chrom_annotations_dir) if "gff" in x.lower() and genome in x.lower()]
            gff = gff[0]
   
            if not gff:
                print("Couldn't find a corresponding GFF for %s in %s" % (self.genome,self.chrom_annotations_dir))
                print("Please place one in %s and run again" % (self.chrom_annotations_dir))
                sys.exit(1)
      
            print("Using %s to create Chromosome Annotations" % (gff)) 
            print("This will only occur once for every new genome")

            command = "%s %s" % (os.path.join(script_dir,"parseGFF.py"),os.path.join(self.chrom_annotations_dir,gff))
            call(command,shell=True)

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

        chromosome_files = [os.path.join(self.genome_annotations_dir,gen) for gen in os.listdir(self.genome_annotations_dir)]
        chromosomes = {x:{} for x in os.listdir(self.genome_annotations_dir)}
        
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

                    try:
                        if geneA_chromosome != "NA":
                            geneAPos,geneANeg = chromosomes[geneA_chromosome][geneA_start]
                        else:
                            geneAPos,geneANeg = ("NotAligned","NotAligned")
                    
                    except KeyError:
                        geneAPos,geneANeg = ("Intergenic","Intergenic")

                    try:
                        if geneB_chromosome != "NA":
                            geneBPos,geneBNeg = chromosomes[geneB_chromosome][geneB_start]
                        else:
                            geneBPos,geneBNeg = ("NotAligned","NotAligned")
                    
                    except KeyError:
                        geneBPos,geneBNeg = ("Intergenic","Intergenic")

                    to_write = [readID] + geneA_info + [geneAPos,geneANeg] + geneB_info + [geneBPos,geneBNeg]
                    out_file.write(" ".join(to_write) + "\n")


    # ---- Get Candidate Reads by joining R1 and R2 together. Filter and keep genes
    # ---- where R1 != R2
    def getCandidateReads(self,debug=False):
        print("Getting Candidate Reads")

        print("\tFiltering R1R2: Only keeping reads where GeneA != GeneB")
        self.filterR1R2("R1R2.gff.out",debug=debug)


    # ---- ---- Helper Function for getCandidateReads
    def filterR1R2(self,R1R2,debug):
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

                    # If overlapping genes Choose the one with the highest ranking
                    # This sorting technique will be used for choosing overlapping strands
                    # and between GeneA and GeneB
                    geneAPos = sorted(geneAPos.split(";"),key = lambda x: self.ranking(x))[0]
                    geneANeg = sorted(geneANeg.split(";"),key = lambda x: self.ranking(x))[0]

                    geneBPos = sorted(geneBPos.split(";"),key = lambda x: self.ranking(x))[0]
                    geneBNeg = sorted(geneBNeg.split(";"),key = lambda x: self.ranking(x))[0]
                    
                    # Now Choose between strands
                    geneAinfo = sorted([geneAPos,geneANeg],key = lambda x: self.ranking(x))[0]
                    geneBinfo = sorted([geneBPos,geneBNeg],key = lambda x: self.ranking(x))[0]


                    if debug:
                        print line.strip(),geneAinfo,geneBinfo

                    geneA = geneAinfo.split(":")[0].split(".")[0]
                    geneB = geneBinfo.split(":")[0].split(".")[0]

                    # Count the amount of times we see genes
                    seen_genes[geneA] = seen_genes.get(geneA,0) + 1
                    seen_genes[geneB] = seen_genes.get(geneB,0) + 1

                    print geneA,geneB

                    if geneA != geneB:
                        candidate_reads.write(" ".join([" ".join(A_line_info),geneAinfo," ".join(B_line_info),geneBinfo + "\n"]))

        # Print out Seen_Genes
        seen_genes_keys = seen_genes.keys()[:]
        seen_genes_keys.sort()

        print("Writing seen genes out to GENE_COUNTS.csv")
        with open("GENE_COUNTS","w") as out_file:
            for gene in seen_genes_keys:
                out_file.write("%s,%s\n" % (gene,seen_genes[gene]))

    @staticmethod
    def ranking(full_gene_info):
        """
        """
        if full_gene_info == "Intergenic":
            return 20

        try:   
            section = full_gene_info.split(":")[1]
        
        except IndexError:
            section = full_gene_info

        if section == "CDS":
            return 1
        
        elif section == "three_prime_UTR":
            return 2
        
        elif section == "five_prime_UTR":
            return 3
        
        else:
            return 10


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
    #loxData.align2gff(debug=True)
    loxData.getCandidateReads()
    # loxData.cleanUp()
