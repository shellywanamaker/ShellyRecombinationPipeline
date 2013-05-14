ShellyRecombinationPipeline
===========================

Shelly Recombination Pipeline

# Prerequisites:
- Bowtie2 is available via command line (ie it's in your bash profile or bashrc)
- You have a CSV of gene pairs that you expect to see

## To Run
Each module can actually run on it's own. The recombinationPipeline is actually just a wrapper script that calls the modules in order.

## NOTE
The chromosome annotations that the script uses are not built by default. However, the script checks the first time it is run and if it cannot find the annotations it will build them from the included TAIR10 GFF3. From then on the script can run without having to build th annotations every time. The script, parseGFF.py, can actually be used on any arbitrary genome's GFF as long as it is in the same format as the TAIR10 GFF3.

### bowtieFastqs.py
- usage -> bowtieFastqs.py /Folder/With/Fastqs/
- You can also optionally run bowtieFastqs without a specified folder and it will assume that the current working directory is the folder with fastqs.
- Although the INDEXES are hardcoded into the module you can change them easily in the __main__ section of either bowtieFastqs or in recombinationPipeline.py
- NOTE: all reads whether they are unaligned or not are kept in the SAM files. This means the SAM files can be quite big.

### getCandidates.py
- The majority of the processing and logic happens here.
- The sam file outputs are slimmed down by removing multi-aligned, unaligned, and low quality reads.
- R1 and R2 files are then full outer joined. If either R1 or R2 didn't align then a "NA" will be in the columns that should have had information.
- Clones are collapsed. Clones are defined as multiple reads with the same R1 and R2 alignment positions.
- Gff alignment. The resulting R1R2.ggf.out file is aligned to a TAIR10 annotation that contains the genes on the positive and negative strand at every position in the chromosome.
- Finally Candidate Reads are called. The strand that has the most likely gene is taken for both R1 and R2. If the gene on R1 != R2 we output that to a Candidate Reads file.

### filterCandidates.py
- With a CSV that contains all the possible expected interactions between Genes this module returns a RESULTS file that contains the number of expected pairs, expected genes pair with unexpected genes, predicted genes in unpredicted pairs, and unpredicted genes.
- Additionally, the script reports the drop off of reads that make it through each step of the pipeline
