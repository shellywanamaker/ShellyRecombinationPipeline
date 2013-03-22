#!/usr/bin/env python

import os
import sys
import signal
from subprocess import call
from subprocess import check_output


class filterData(object):

    def __init__(self,csv_file,complement=False):

        self.csv        = csv_file
        self.complement = False

    def compareCandidateReads2Predicted(self):
        """
        """
        predictedPairsCSV = self.csv

        predictedGenes = set()
        predictedPairs = set()

        predictedGenesInReads = set()

        alignedGenes   = set()
        alignedPairs   = smartDict()

        predictedGenePairedWithUnpredictedGene                     = smartDict()
        predictedGenesInPredictedPair                              = smartDict()
        predictedGenesPairedWithPredcitedGeneButPairIsNotPredicted = smartDict()
        unpredictedGenePairs                                       = smartDict()
        unpredictedGenes                                           = set()

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

        if not predictedPairsCSV:
            print("\nYou didn't give me a csv file!\n")
            sys.exit(1)

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

                    unpredictedGenes.add(geneA)
                    unpredictedGenes.add(geneB)
                    unpredictedGenePairs.add(geneA,geneB)
                
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

        # Pull out unpredicted Pairs
        unpredictedGenePairsList = []

        for g1 in unpredictedGenePairs:

            for g2 in unpredictedGenePairs[g1]:

                l = [g1,g2]
                l.sort()

                tup = (l[0],l[1],unpredictedGenePairs[g1][g2])

                unpredictedGenePairsList.append(tup)

        unpredictedGenePairsList.sort()


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

            results.write("# Results:\n")
            results.write("\n# Overall Number of reads: %s" % self.count_of_reads)
            results.write("\n# Count of R1 reads that aligned to genome: %s" % self.bowtie_R1_alignments)
            results.write("\n# Count of R2 reads that aligned to genome: %s" % self.bowtie_R2_alignments)
            results.write("\n# Count of R1 reads after clone removal:    %s" % self.R1_slim_clean_count)
            results.write("\n# Count of R2 reads after clone removal:    %s" % self.R2_slim_clean_count)
            results.write("\n# Count of reads after joining R1 and R2:   %s" % self.R1R2_join_count)
            results.write("\n# Count of Candidate Reads                  %s" % self.Candidate_Reads_count)

            results.write("\n\n# Number of predicted genes: %s" % len(predictedGenes))
            results.write("\n# Number of predicted genes in reads: %s" % len(predictedGenes.difference(alignedGenes)))
            results.write("\n# Genes not found: \n")

            genesNotFound = [x for x in predictedGenes.difference(predictedGenes.difference(alignedGenes))]
            results.write(" ".join(genesNotFound))

            results.write("\n\n# Number of predicted gene pairs: %s" % len(predictedPairs))
            results.write("\n# Number of predicted gene pairs found in reads: %s" % len(predictedGenesInPredictedPair))

            # Overal Print outs
            results.write("\n\n# Predicted Genes in a Predicted Pair\n")
            for tup in predictedGenesInPredictedPairList:
                results.write(" ".join([tup[0],tup[1],str(tup[2]) + "\n"]))

            results.write("\n\n# Predicted Gene Paired with unpredicted Gene:\n")

            for tup in predictedGenePairedWithUnpredictedGeneList:
                check = self.checkAcessionNumber(tup[1],predictedGenes)

                results.write(" ".join([tup[0],tup[1],str(tup[2]),check + "\n"]))

            results.write("\n\n# Predicted Genes but the Pair itself is not Predicted:\n")
            for tup in predictedGenesPairedWithPredcitedGeneButPairIsNotPredictedList:
                results.write(" ".join([tup[0],tup[1],str(tup[2]) + "\n"]))

            results.write("\n\n# Genes aligned but not predicted\n")
            genes_aligned_but_not_predicted = list(predictedGenes.difference(alignedGenes))
            genes_aligned_but_not_predicted.sort()

            # for gene in genes_aligned_but_not_predicted:
            #     check = self.checkAcessionNumber(gene,predictedGenes)
            #     results.write(gene + check + "\n")

            results.write("\n\n# Genes pairs aligned but not Predicted\n")
            for tup in unpredictedGenePairsList:
                check1 = self.checkAcessionNumber(tup[0],predictedGenes)
                check2 = self.checkAcessionNumber(tup[1],predictedGenes)
                
                results.write(" ".join([tup[0],tup[1],str(tup[2]),check1 + check2 + "\n"]))

    @staticmethod
    def checkAcessionNumber(gene2check,predictedGenes,acession_range=10):

        if gene2check != "Intergenic":
            acession_number = int(gene2check.split("G")[1])
            head            = gene2check.split("G")[0] + "G"


            to_check = [head + "%05d" % x for x in range(acession_number-10, acession_number +10 + 1)]
            for gene in to_check:

                if gene in predictedGenes:
                    return "! %s" % (gene)

            else:
                return ""

        else:
            return ""

    def cleanUp(self):
        print("Cleaning up Current Directory")
        command = "rm R1.slim* R2.slim*"
        call(command,shell=True)


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


if __name__=="__main__":
    """
    Script Expects:

    1) There is a Candidate_Reads.txt file in the directory where the script is Run
    2) This is one and only one .csv file in the diretcory where the script is run and it is the GenePairs.csv
    """

    # ---- Prepping for Run
    candidate_reads = os.path.join(os.getcwd(),"Candidate_Reads.txt")

    if not os.path.isfile(candidate_reads):
        print("\nCould not find Candidate_Reads.txt in your current folder\nBye!\n")
        sys.exit(1)

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

    # ---- Running
    f = filterData(csv_file)

    f.compareCandidateReads2Predicted()

