TODO: Pick a better title for this program
==========================================

Table of contents
-----------------

NAME
SYNOPSIS
DESCRIPTION
OPTIONS
INVOCATION
SEE ALSO
FILES
AUTHORS
BUG REPORTS

## 1. NAME


TODO - Synthetic read generator for T-cell receptor data

## 2. SYNOPSIS


  ./lib/main.py [options] [allele fasta files]

## 3. DESCRIPTION


TODO



## 4. OPTIONS

```

usage: main.py [-h] [--chr7-filename FILE] [--chr14-filename FILE]
               [--tcell-data FILE] [--output BASENAME] [--repertoire-size N]
               [--population-size N]
               [--population-distribution {gaussian,chisquare}]
               [--population-gaussian-parameters N | --population-chisquare-parameters k:cutoff]
               [--read-type {paired,single}] [--sequence-type {dna,rna}]
               [--sequence-count SEQ] [--read-length-mean READ_LENGTH_MEAN]
               [--read-length-sd READ_LENGTH_SD] [--read-length-sd-cutoff N]
               [--insert-length-mean INSERT_LENGTH_MEAN]
               [--insert-length-sd INSERT_LENGTH_SD]
               [--insert-length-sd-cutoff N]
               [--degrade-logistic B:L:k:mid | --degrade-phred PHRED_STRING]
               [--degrade-variability FLOAT] [--display-degradation]
               [--receptor-ratio RATIO]
               [--log-level {debug,info,warning,error,critical}]
               [ALLELE_FILES [ALLELE_FILES ...]]

Generate synthetic TCR read data

positional arguments:
  ALLELE_FILES          Read TCR allele data from one or more IMGT-formatted
                        fasta files

optional arguments:
  -h, --help
			show this help message and exit
  --chr7-filename FILE
			Filename of a FASTA formatted file with chromosome 7
                        reference data
  --chr14-filename FILE
                        Filename of a FASTA formatted file with chromosome 7
                        reference data
  --tcell-data FILE
			Filename of a tab-separated file with T cell receptor
                        segments and reference gene coordinates
  --output BASENAME
			Name of output fastq file. This should be a basename,
                        e.g. '--output=foo' will write to 'foo.fastq',
                        'foo.statistics.csv', etc. Default is 'tcr_synth'
  --repertoire-size N
			Size of the TCR repertoire (i.e. the number of unique
                        TCRs that are generated). Default is 10
  --population-size N
			The number of T-cells in the repertoire (e.g. if
                        repertoire-sze=5 and population-size=15, then there
                        are 3 clones of each unique TCR on average). Default
                        is 100
  --population-distribution {gaussian,chisquare}
                        Population distribution function. This defines the
                        function used to distribute the population among the
                        repertoire. Default is the (normalized) gaussian. See
                        --population-gaussian-parameters, --population-
                        chisquare-parameters
  --population-gaussian-parameters N
                        Parameter for the normalized gaussian distribution.
                        The number of standard deviations to include in our
                        population distribution. Decimal value. Default is 3
  --population-chisquare-parameters k:cutoff
                        Parameters for the chi-square distribution. Takes an
                        argument formatted as 'k:cutoff', where k - degrees of
                        freedom. Default is 3. cutoff - X-axis maximum.
                        Default is 8
  --read-type {paired,single}
                        Generate either single or paired-end reads
  --sequence-type {dna,rna}
                        Generate sequences from simulated DNA or RNA
  --sequence-count SEQ  Number of sequences (reads) to generate. Default is
                        1000
  --read-length-mean READ_LENGTH_MEAN
                        The average length of reads in nucleotides. Default is
                        48
  --read-length-sd READ_LENGTH_SD
                        The SD of read length variation in nucleotides. Set to
                        zero for fixed-length reads. Default is 4
  --read-length-sd-cutoff N
                        Read lengths are restricted to less than N standard
                        deviations from the mean. Default is 4
  --insert-length-mean INSERT_LENGTH_MEAN
                        The average length of the insert for paired end reads.
                        Default is 48
  --insert-length-sd INSERT_LENGTH_SD
                        The standard deviation of insert length variation in
                        nucleotides. Set to zero for fixed-length inserts.
                        Default is 4
  --insert-length-sd-cutoff N
                        Insert lengths are restricted to less than N standard
                        dviations from the mean. Default is 4
  --degrade-logistic B:L:k:mid
                        Simulate non-optimal quality using the logstic
                        (sigmoid) function. Takes an argument formatted as
                        'B:L:k:mid'. B - Base error rate probability. L -
                        Maximum error rate. k - Steepness factor. mid -
                        Midpoint, this is the base position where error rate
                        is equal to 1/2 of L. Default is off. This option is
                        mutually exclusive to --degrade-phred. See: --display-
                        degradation, --degrade-variability
  --degrade-phred PHRED_STRING
                        Simulate non-optimal quality using a Phred+33 string
                        to specify quality on a per-nucleotide basis. If a
                        generated read is longer than the given phred string,
                        then the last character in the phred string is used.
                        Default is off. This option is mutually exclusive to
                        --degrade-logistic. See: --degrade-variability
  --degrade-variability FLOAT
                        Applies a relative variability in the per-nucleotide
                        error applied by the --degrade option. If a given base
                        were to have an error rate of 0.1 (10%), then a
                        degrade-variability of 0.5 (50%) would result in an
                        error rate in the range of 0.1 +/- 0.1 * 0.5. Default
                        is 0
  --display-degradation
                        Display the error rate per base pair for a given
                        B:L:k:mid value and exit. The number of positions
                        displayed is adjustable through the --read-length-mean
                        option. This is mostly useful in adjusting these
                        parameters to be passed to the --degrade option. Note
                        that no reads or repertoire will be generated when
                        this option is given
  --receptor-ratio RATIO
                        Ratio of alpha/beta vs gamma/delta sequences. Default
                        is 0.9 (9 alpha/beta per 1 gamma/delta TCR)
  --log-level {debug,info,warning,error,critical}
                        Logging level. Default is warning and above

```

### 4.1 File I/O options

    --output BASENAME

Write output files using `BASENAME` as a prefix.  Can be a directory followed by a basename.

`--chr7-filename FILE`
`--chr14-filename FILE`

Read reference chromosome data from FASTA-formatted files.  Required when invoking to generate read data



| Option                | Defaults  | Description                               |
|-----------------------|-----------|-------------------------------------------|
| `--output BASENAME`   |`tcr_synth`| Write output files using BASENAME as a prefix.  Can be a directory followed by a basename.|
| `--chr7-filename FILE`| None      | Location of the FASTA-formatted file with the reference chromosome 7 |
| `--chr14-filename FILE`| None     | Location of the FASTA-formatted file with the reference chromosome 14 |
|                       |           |                                           |
|                       |           |                                           |
|                       |           |                                           |
|                       |           |                                           |
|                       |           |                                           |
|                       |           |                                           |
|                       |           |                                           |
|                       |           |                                           |
|                       |           |                                           |
|                       |           |                                           |
|                       |           |                                           |
|                       |           |                                           |
|                       |           |                                           |
|                       |           |                                           |
|                       |           |                                           |
|                       |           |                                           |



 ALLELE_FILES          Read TCR allele data from one or more IMGT-formatted
                        fasta files

optional arguments:
  -h, --help            show this help message and exit
  --chr7-filename FILE  Filename of a FASTA formatted file with chromosome 7
  --chr14-filename FILE
  --tcell-data FILE     Filename of a tab-separated file with T cell receptor
  --output BASENAME     Name of output fastq file. This should be a basename,
  --repertoire-size N   Size of the TCR repertoire (i.e. the number of unique
  --population-size N   The number of T-cells in the repertoire (e.g. if
  --population-distribution {gaussian,chisquare}
  --population-gaussian-parameters N
  --population-chisquare-parameters k:cutoff
  --read-type {paired,single}
  --sequence-type {dna,rna}
  --sequence-count SEQ  Number of sequences (reads) to generate. Default is
  --read-length-mean READ_LENGTH_MEAN
  --read-length-sd READ_LENGTH_SD
  --read-length-sd-cutoff N
  --insert-length-mean INSERT_LENGTH_MEAN
  --insert-length-sd INSERT_LENGTH_SD
  --insert-length-sd-cutoff N
  --degrade-logistic B:L:k:mid
  --degrade-phred PHRED_STRING
  --degrade-variability FLOAT
  --display-degradation
  --receptor-ratio RATIO
  --log-level {debug,info,warning,error,critical}









## 5. INVOCATION

### 5.1 Quick usage



## 6. SEE ALSO

## 7. FILES

### 7.1 Input files
1. TCR FASTA files were downloaded from the ImMunoGeneTics GENE-DB site
(www.imgt.org/genedb).  The datasets were chosen to include "nucleotide 
sequences for F+ORF+all P".  As of this writing, a description of the available
datasets can be found at http://www.imgt.org/genedb/directlinks

    An example of these files can be found at:
http://www.imgt.org/genedb/GENElect?query=7.2+TRAV&species=Homo+sapiens

2. TCR segment descriptors (e.g. TRBV22-1 L & V-gene unit is on forward strand of 7q34 at nucleotides 142641746->142642235) were copied from the "Localization in genome assemblies" tables for the various TCR loci.  As of this writing, these can be found at http://www.imgt.org/genedb/

    An example of this table can be found at:
http://www.imgt.org/genedb/resultPage_localizations.action?mode=localizations&collectionId=4&regionId=1

3. Reference genome data for chromosomes 7 and 14.

#### 7.1.1 Modifications to input files
N.b. that the versions of these files were hand-curated to exclude various logical inconsistencies that would interfere with the operation of this software.  These exclusions and edits fall into several broad categories:

1. V-REGIONs without L-V-GENE-UNIT regions
2. Orphons and pseudogenes
3. V-REGIONs without L-PART1+L-PART2 data (which is needed to generate RNA)
4. Addition of explicit L-PART1-L-PART2 reference coordinates
5. Addition of D-REGION reference coordinates
6. Correction of NCBI-format reference coordinates


## 8. AUTHORS


## 9. BUGS
### 9.1 Reporting Bugs
### 9.2 Known issues/limitations