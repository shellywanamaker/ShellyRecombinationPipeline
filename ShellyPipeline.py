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

    def bowtie(self,options="--local -p 2",indexes="/home/feeney/bin/bowtie2/INDEXES/tair10"):
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

        # Remove intermediate Slim Files
        remove_slims = "rm *.slim"
        call(remove_slims,shell=True)

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
        print("Cleaning up Temporary Files and Folders")
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
        # counts
        clean_temp = "rm *.clean"
        call(clean_temp,shell=True)

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

    def compareCandidateReads2Predicted(self,predictedPairsCSV=None,predictedGenesCSV=None):
        """
        """

        predictedGenes = set()
        predictedPairs = set()

        predictedGenesInReads = set()

        alignedGenes   = set()
        alignedPairs   = smartDict()

        predictedGenePairedWithUnpredictedGene                     = smartDict()
        predictedGenesInPredictedPair                              = smartDict()
        predictedGenesPairedWithPredcitedGeneButPairIsNotPredicted = smartDict()

        # Adding predicted Genes to a set of both Pairs and Genes themselves
        if predictedPairsCSV:
            with open(predictedPairsCSV) as predictedPairCSV_obj:

                for line in predictedPairCSV_obj:
                    row = line.strip().split(",")

                    geneA = row[0]
                    geneB = row[1]

                    if geneA == "" or geneB == "":
                        continue

                    predictedGenes.add(geneA)
                    predictedGenes.add(geneB)

                    predictedPairs.add((geneA,geneB))

        if predictedGenesCSV:
            print "poop"

        # Going Through Candidate Reads
        with open("Candidate_Reads.txt","r") as candidate_reads:

            for line in candidate_reads:
                row = line.strip().split()

                readID = row[0]
                geneA  = row[6].split(":")[0].split(".")[0]
                geneB  = row[12].split(":")[0].split(".")[0]

                alignedGenes.add(geneA)
                alignedGenes.add(geneB)

                # Checking the Pairs!
                if (geneA,geneB) not in alignedPairs and (geneB,geneA) not in alignedPairs:
                    alignedPairs.add(geneA,geneB)

                elif (geneA,geneB) in alignedPairs and (geneB,geneA) not in alignedPairs:
                    alignedPairs.add(geneA,geneB)

                elif (geneA,geneB) not in alignedPairs and (geneB,geneA) in alignedPairs:
                    alignedPairs.add(geneB,geneA)

                else:
                    print("Whoops! That shouldnt happen")
                    alignedPairs(geneA,geneB)

                # Checking the Individual Genes
                if   geneA in predictedGenes and geneB in predictedGenes and (geneA,geneB) in predictedPairs:
                    # both genes predicted and pair A,B is predicted
                    predictedGenesInPredictedPair.add(geneA,geneB)
                
                elif geneA in predictedGenes and geneB in predictedGenes and (geneB,geneA) in predictedPairs:
                    # both genes predicted and pair B,A is predicted
                    predictedGenesInPredictedPair.add(geneB,geneA)

                elif geneA in predictedGenes and geneB in predictedGenes and ((geneA,geneB) not in predictedPairs or (geneB,geneA) not in predictedPairs):
                    # both genes predicted pair is not predicted
                    # Write a check to see if they are close in acession number to another gene.
                    if geneA in predictedGenesPairedWithPredcitedGeneButPairIsNotPredicted and geneB in predictedGenesPairedWithPredcitedGeneButPairIsNotPredicted[geneA]:
                        predictedGenesPairedWithPredcitedGeneButPairIsNotPredicted.add(geneA,geneB)

                    elif geneB in predictedGenesPairedWithPredcitedGeneButPairIsNotPredicted and geneA in predictedGenesPairedWithPredcitedGeneButPairIsNotPredicted[geneB]:
                        predictedGenesPairedWithPredcitedGeneButPairIsNotPredicted.add(geneB,geneA)

                    else:
                        predictedGenesPairedWithPredcitedGeneButPairIsNotPredicted.add(geneA,geneB)

                elif geneA in predictedGenes and geneB not in predictedGenes:
                    # A is predicted but B is not
                    predictedGenePairedWithUnpredictedGene.add(geneA,geneB)

                elif geneA not in predictedGenes and geneB in predictedGenes:
                    # A is not predicted but B is.
                    predictedGenePairedWithUnpredictedGene.add(geneB,geneA)

                elif geneA not in predictedGenes and geneB not in predictedGenes:
                    # A is not and B is not
                    # This is already taken care of atht eh top! they are both added to aligned genes automatically
                    continue
                
                else:
                    print("Oh fuck I missed a use case")

        # ------------------------- Prepare Output ---------------------------- #
           
        predictedGenePairedWithUnpredictedGeneList = []

        for g1 in predictedGenePairedWithUnpredictedGene:

            for g2 in predictedGenePairedWithUnpredictedGene[g1]:
                tup = (g1,g2,predictedGenePairedWithUnpredictedGene[g1][g2])

                predictedGenePairedWithUnpredictedGeneList.append(tup)

        predictedGenePairedWithUnpredictedGeneList.sort()

        # Sorting predicted pair list
        predictedGenesInPredictedPairList = []

        for g1 in predictedGenesInPredictedPair:

            for g2 in predictedGenesInPredictedPair[g1]:

                tup = (g1,g2,predictedGenesInPredictedPair[g1][g2])

                predictedGenesInPredictedPairList.append(tup)

        predictedGenesInPredictedPairList.sort()

        # Predicted Genes but the pair is not predicted.
        predictedGenesPairedWithPredcitedGeneButPairIsNotPredictedList = []
        test1 = set()
        test2 = set()

        for g1 in predictedGenesPairedWithPredcitedGeneButPairIsNotPredicted:

            for g2 in predictedGenesPairedWithPredcitedGeneButPairIsNotPredicted[g1]:

                test1.add(g1)
                test2.add(g2)

                tup = (g1,g2,predictedGenesPairedWithPredcitedGeneButPairIsNotPredicted[g1][g2])

                predictedGenesPairedWithPredcitedGeneButPairIsNotPredictedList.append(tup)

        predictedGenesPairedWithPredcitedGeneButPairIsNotPredictedList.sort()

        # Since Sam files have all reads. Count raw sam count gives number of input reads
        print("Counting Reads...")
        count_input_reads_command = "wc -l ./alignments/bowtie.R1.out.sam | awk '{print $1}'"
        self.count_of_reads = check_output(count_input_reads_command,shell=True).strip()

        bowtie_R1_alignments_command = """cat ./alignments/bowtie.R1.out.sam | awk '{count++; if(count > 9 && $3 != "*") print $0}' | wc -l | awk '{print $1}'"""
        self.bowtie_R1_alignments = check_output(bowtie_R1_alignments_command,shell=True).strip()

        bowtie_R2_alignments_command = """cat ./alignments/bowtie.R2.out.sam | awk '{count++; if(count > 9 && $3 != "*") print $0}' | wc -l | awk '{print $1}'"""
        self.bowtie_R2_alignments = check_output(bowtie_R2_alignments_command,shell=True).strip()

        Candidate_Reads_count_command = "cat Candidate_Reads.txt | wc -l | awk '{print $1}'"
        self.Candidate_Reads_count = check_output(Candidate_Reads_count_command,shell=True).strip()
        
        R1_slim_clean_count_command = "cat R1.slim.clean | wc -l | awk '{print $1}'"
        self.R1_slim_clean_count = check_output(R1_slim_clean_count_command,shell=True).strip()

        R2_slim_clean_count_command = "cat R2.slim.clean | wc -l | awk '{print $1}'"
        self.R2_slim_clean_count = check_output(R2_slim_clean_count_command,shell=True).strip()

        R1R2_join_count_command = "cat R1R2 | wc -l | awk '{print $1}'"
        self.R1R2_join_count = check_output(R1R2_join_count_command,shell=True).strip()

        # ------------------------- Output ---------------------------------- #
        print("Writing RESULTS to disk.")

        with open("RESULTS.txt","w") as results:

            results.write("Results:\n")
            results.write("\nOverall Number of reads: %s" % self.count_of_reads)
            results.write("\nCount of R1 reads that aligned to genome: %s" % self.bowtie_R1_alignments)
            results.write("\nCount of R2 reads that aligned to genome: %s" % self.bowtie_R2_alignments)
            results.write("\nCount of R1 reads after clone removal:    %s" % self.R1_slim_clean_count)
            results.write("\nCount of R2 reads after clone removal:    %s" % self.R2_slim_clean_count)
            results.write("\nCount of reads after joining R1 and R2:   %s" % self.R1R2_join_count)
            results.write("\nCount of Candidate Reads                  %s" % self.Candidate_Reads_count)

            results.write("\n\nNumber of predicted genes: %s" % len(predictedGenes))
            results.write("\nNumber of predicted genes in reads: %s" % len(predictedGenes.difference(alignedGenes)))
            results.write("\nGenes not found: ")

            genesNotFound = [x for x in predictedGenes.difference(alignedGenes)]
            results.write(" ".join(genesNotFound))

            results.write("\n\nNumber of predicted gene pairs: %s" % len(predictedPairs))
            results.write("\nNumber of predicted gene pairs found in reads: %s" % len(predictedGenesInPredictedPair))

            # Overal Print outs
            results.write("\n\nPredicted Genes in a Predicted Pair\n")
            for tup in predictedGenesInPredictedPairList:
                results.write(" ".join([tup[0],tup[1],str(tup[2]) + "\n"]))

            results.write("\n\nPredicted Gene Paired with unpredicted Gene:\n")

            for tup in predictedGenePairedWithUnpredictedGeneList:
                check = self.checkAcessionNumber(tup[1],predictedGenes)

                results.write(" ".join([tup[0],tup[1],str(tup[2]),check + "\n"]))

            results.write("\n\nPredicted Genes but the Pair itself is not Predicted:\n")
            for tup in predictedGenesPairedWithPredcitedGeneButPairIsNotPredictedList:
                results.write(" ".join([tup[0],tup[1],str(tup[2]) + "\n"]))

    def checkAcessionNumber(self,gene2check,predictedGenes,acession_range=10):

        if gene2check != "Intergenic":
            acession_number = int(gene2check.split("G")[1])
            head            = gene2check.split("G")[0] + "G"

            # Doesn't work!
            to_check_above = [head + str(x) for x in range(acession_number + 10,acession_number + 10 + (acession_range - 10))]
            to_check_below = [head + str(x) for x in range(acession_number - 10 + (-acession_range + 10),acession_number-10)]
            to_check = to_check_below + to_check_above

            for ac in to_check:
                print head + ac

                if head + ac in predictedGenes:
                    print "poop"
                    return "!"
            else:
                return ""

        else:
            return ""


class smartDict(dict):

    def __str__(self):

        stringToReturn = []

        for key in self.keys():
            stringToReturn.append(str(key))
            stringToReturn.append("  ")
            stringToReturn.append(str(self[key]))
            stringToReturn.append("\n")
        
        return "".join(stringToReturn[:len(stringToReturn)-1])

    def add(self,key,value):
        """
        key is added as a key to the dictionary and value is added to the value dictionary.

        If Value is not in the Value dictionary associated with the key it is added and it's frequency is set to one. Else the frequency of that Value is updated
        """
        if self.get(key,0) == 0:
            self[key] = {}

        self[key][value] = self[key].get(value,0) + 1

    def sort(self,highToLow = True):
        """
        Sorts the dict by how many values a key has

        Returns that sorted list
        """
        countOfKeysAndTheirValues = self.count()

        sortedList = []

        for key in sorted(countOfKeysAndTheirValues,key = itemgetter(1), reverse = highToLow):
            sortedList.append(str(key[0]) + "  " + str(key[1]))
            sortedList.append("\n")

        return "".join(sortedList[:len(sortedList)-1])

    def sortAndDisplayValues(self,highToLow = True):
        """
        Displays keys sorted by how values it has and shows that value

        Returns: list of tuples
        """

        countOfKeysAndTheirValues = self.count()

        sortedList = []

        for key in sorted(countOfKeysAndTheirValues, key = itemgetter(1), reverse = highToLow):
            # sortedList.append(str(key[0]) + ":" + " " + str(key[1]) + " " + str(self[key[0]]))

            sortedList.append(str(key[0]))
            sortedList.append(":")
            sortedList.append(" ")
            sortedList.append(str(key[1]))
            sortedList.append(" ")

            valueList = self[key[0]].items()

            for value in sorted(valueList, key = itemgetter(1), reverse = True):
                sortedList.append("(")
                sortedList.append(str(value[0]))
                sortedList.append(":")
                sortedList.append(str(value[1]))
                sortedList.append(")")
                sortedList.append(" ")


            sortedList.append("\n")

        return "".join(sortedList[:len(sortedList)-1])

    def count(self):
        """
        Counts how many items (and their frequencies if more than one) a key has

        returns a list of tuples with the keys and their counts
        """
        countOfKeysAndTheirValues = []

        for key in self.keys():
            count = 0

            for valueKey in self[key].keys():
                count += self[key][valueKey]

            countOfKeysAndTheirValues.append((key,count))

        return countOfKeysAndTheirValues

def startupChecks():

    # Check for Bowtie2
    if not check_PATH_for_program("bowtie2"):
        print("\nBowtie2 is not in your $PATH. Can you modify your .bash_profile or .bashrc, please?\n")
        sys.exit(1)

    # Check for a .CSV file for a Candidate (if any) reads

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
    # path    = startupChecks()

    # loxData = loxData(path)

    # loxData.bowtie()
    # loxData.slim_and_clean_sam_files()
    # loxData.align2gff()
    # loxData.getCandidateReads()
    # loxData.compareCandidateReads2Predicted(predictedPairsCSV="PRSPosGenePairs.csv")

    # --------- Testing ------------#
    path    = startupChecks()

    loxData = loxData(path)
    loxData.compareCandidateReads2Predicted(predictedPairsCSV="PRSPosGenePairs.csv")