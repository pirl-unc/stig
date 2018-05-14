STIG: Synthetic TCR Informatics Generator
=========================================

Current for v1.2.1

Table of contents
-----------------

1. NAME
2. SYNOPSIS
3. DESCRIPTION
4. OPTIONS
5. INVOCATION
6. SEE ALSO
7. FILES
8. AUTHORS
9. BUG REPORTS

## 1. NAME


STIG - Synthetic TCR Informatics Generator

## 2. SYNOPSIS


  ./lib/stig [options] [allele fasta files]

## 3. DESCRIPTION

STIG is a tool for creating artificial T-cell repertoires and producing simulated sequencing data from them.  Many characteristics of the repertoires and the sequencing output can be customized.  Reads can be generated in both RNA and DNA space.  Applications include evaluating and optimizing tools for performing analysis of T-cell receptors.

### 3.1 Conventions

STIG utilizes a number of underlying assumptions when reading input and output data:

1. Quality strings are Phred+33 (Sanger) format: (0,40) or (!, I), with I being highest quality
2. DNA and RNA inputs and outputs are 5' --> 3' direction, even when generating paired end reads or specififying amplicon probes


## 4. OPTIONS

```
usage: stig [-h] [--chr7-filename FILE] [--chr14-filename FILE]
            [--tcell-data FILE] [--output BASENAME] [--load-population FILE]
            [--repertoire-size N] [--repertoire-unique]
            [--repertoire-chain-unique] [--repertoire-cdr3-unique]
            [--population-size N]
            [--population-distribution {gaussian,chisquare}]
            [--population-gaussian-parameters N | --population-chisquare-parameters k:cutoff]
            [--read-type {paired,single,amplicon}] [--sequence-type {dna,rna}]
            [--sequence-count N] [--read-length-mean READ_LENGTH_MEAN]
            [--read-length-sd READ_LENGTH_SD] [--read-length-sd-cutoff N]
            [--insert-length-mean INSERT_LENGTH_MEAN]
            [--insert-length-sd INSERT_LENGTH_SD]
            [--insert-length-sd-cutoff N] [--amplicon-probe STR]
            [--degrade-logistic B:L:k:mid | --degrade-phred PHRED_STRING | --degrade-fastq FILE[,FILE2]
            | --degrade-fastq-random FILE[,FILE2]]
            [--degrade-variability FLOAT] [--display-degradation]
            [--receptor-ratio RATIO]
            [--log-level {debug,info,warning,error,critical}]
            [ALLELE_FILES [ALLELE_FILES ...]]

Generate synthetic TCR read data

positional arguments:
  ALLELE_FILES          Read TCR allele data from one or more IMGT-formatted
                        fasta files

optional arguments:
  -h, --help            show this help message and exit
  --chr7-filename FILE  Filename of a FASTA formatted file with chromosome 7
                        reference data
  --chr14-filename FILE
                        Filename of a FASTA formatted file with chromosome 7
                        reference data
  --tcell-data FILE     Filename of a tab-separated file with T cell receptor
                        segments and reference gene coordinates
  --output BASENAME     Basename for output files, e.g. '--output=foo' will
                        write to 'foo.fastq', 'foo.statistics.csv', etc.
                        Default is 'stig.out'
  --load-population FILE
                        Load TCR population and repertoire data from FILE,
                        rather than generating from scratch
  --repertoire-size N   Size of the TCR repertoire (i.e. the number of unique
                        TCRs that are generated). Default is 10
  --repertoire-unique   Force each TCR to be unique on the RNA level. Default
                        is to allow collisions
  --repertoire-chain-unique
                        Force each TCR chain (e.g. alpha) to be unique on the
                        RNA level. Implies unique TCRs as per --repertoire-
                        unique. Default is to allow collisons
  --repertoire-cdr3-unique
                        Force each CDR3 of each chain to be unique on the
                        nucleotide level. Implies unique TCRs as per
                        --repertoire-unique and unique chains as per
                        --repertoire-chain-unique. Note this may cause
                        performance issues as repertoire size increases.
                        Default is to allow collisons
  --population-size N   The number of T-cells in the repertoire (e.g. if
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
  --read-type {paired,single,amplicon}
                        Generate either single, paired-end, or amplicon reads.
                        Default is single
  --sequence-type {dna,rna}
                        Generate sequences from simulated DNA or RNA. Default
                        is DNA
  --sequence-count N    Number of sequences (reads) to generate. Default is
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
                        deviations from the mean. Default is 4
  --amplicon-probe STR  Anchoring/priming sequence for generating amplicon
                        reads. This should align with some RNA or DNA
                        sequence, either sense or anti-sense. Read 1 will have
                        length given by --read-length-* options. Read 2 will
                        be complementary to read 1 and of an identical length.
                        The default value is a 27-mer that anchors on the
                        reverse strand in EX1 of the beta chain C-region
  --degrade-logistic B:L:k:mid
                        Simulate non-optimal quality using the logistic
                        (sigmoid) function. Takes an argument formatted as
                        'B:L:k:mid'. B - Base error rate probability. L -
                        Maximum error rate. k - Steepness factor. mid -
                        Midpoint, this is the base position where error rate
                        is equal to 1/2 of L. Default is off. This option is
                        mutually exclusive to --degrade-phred, --degrade-
                        fastq, and --degrade-fastq-random. See: --degrade-
                        variability
  --degrade-phred PHRED_STRING
                        Simulate non-optimal quality using a Phred+33 string
                        to specify quality on a per-nucleotide basis. If a
                        generated read is longer than the given phred string,
                        then the last character in the phred string is used.
                        Default is off. This option is mutually exclusive to
                        --degrade-logistic, --degrade-fastq and --degrade-
                        fastq-random. See: --degrade-variability
  --degrade-fastq FILE[,FILE2]
                        Simulate non-optimal quality by degrading reads based
                        on Phred+33 quality strings from the given fastq FILE,
                        or files FILE1,FILE2. Two files required when
                        generating paired or amplicon reads. Output quality
                        strings are assigned from FILE in a stepwise fashion
  --degrade-fastq-random FILE[,FILE2]
                        Simulate non-optimal quality by degrading reads based
                        on Phred+33 quality strings from the given fastq FILE,
                        or files FILE1,FILE2. Two files required when
                        generating paired or amplicon reads. Output quality
                        strings are assigned from FILE randomly
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

Please see manual or README for further details
```

## 4.1 Options overview

The majority of options fall into a few broad categories:

1. File input/output options, which determine how alleles, reference chromosomes, and T cell receptor reference data is stored
2. Generation of the T cell receptor repertoire.  Determines what types and how many T cell receptors are created
3. Data generation, which controls if reads are DNA or RNA space, and the read depth, read length, etc. of those reads.  Options also exist to 'degrade' the data in a manner that simlulates less-than-perfect read quality


## 5. INVOCATION

### 5.1 Quick usage



### 5.2 Degrading reads


By default, STIG generates outputs with 100% accuracy against the underlying simulated TCR.  These reads are placed in the output fastq file(s) (by default 'stig.out.fastq', but can be overridden by the `--output` option).  n.b. that because the Phred+33 standard does not have a "0% chance of error" score, STIG will label these reads as of the highest quality available ('I', by Phred+33).

It is frequently useful to have less-than-perfect quality for testing purposes, as this simulates real-world data.  STIG has several mechanisms for degrading the output to match certain criteria.  Degraded outputs are given the suffix `.degraded` in the fastq output filename.  All reads are labeled with a `readnum` value that corresponds between degraded and non-degraded fastq files, so one can see exactly what changes were made to the output based on the quality string.

Note that just because a given position has a low quality score does not necessarily mean that the read doesn't match the underlying DNA/RNA.  STIG will randomly assign bases at each position based on the probabilities given by the quality string.  For example, a location with a degraded quality score of `I` will have an error probability of 0.001, or 0.1%.  If STIG determines a location should contain a simulated error, then the location is replaced with a random nucleotide.  Of course, there is a 25% chance it will substitute the same base value at that location.

#### 5.2.1 Displaying degredation

For tuning `--degrade-logistic` and `--degrade-variability` parameters, the `--display-degradation` option is available.  When given, STIG will not read input or generate outputs, but instead will display an example quality string and display the error rates per position.

#### 5.2.2 Degradation with a logistic function

Specified by the `--degrade-logistic=B:L:k:mid` option.  This uses a logistic (sigmoid) curve to degrade outputs based on their length.  The function parameters B:L:k:mid refer to the base error rate, maximum error rate, steepness factor, and midpoint length, respectively.  The `--display-degradation` option is often useful for tuning these options.

The `--degrade-variability` option can also be applied to introduce additional error.

This option is mutually exclusive to `--degrade-phred`, `--degrade-fastq` and `--degrade-fastq-random`.

#### 5.2.3 Degradation specified by a Phred+33 string

Specified by the `--degrade-phred=PHRED` option.  This uses the error probabilites speficied by the input string and applies them to the generated reads, where the error probability of the nth read in an output is given by the nth quality score of `PHRED`.  If the generated read is longer than `PHRED`, then the last quality score in `PHRED` is used for the remaining bases.

The `--degrade-variability` option can also be applied to introduce additional error.

This option is mutually exclusive to `--degrade-logistic`, `--degrade-fastq` and `--degrade-fastq-random`.

#### 5.2.4 Degradation specified by fastq quality scores

Specified by the `--degrade-fastq=FILE[,FILE2]` option.  This uses the error probabilities specified by reads in a fastq formatted file `FILE` and applies them to the generated reads.  This is done in a stepwise fashion, such that the nth read generated has the quality score of the nth read in `FILE`.  In the event that more reads are generated than reads in `FILE`, the quality strings will be recycled starting from the start of `FILE`.

For paired-end or amplicon reads one may specify a second file, `FILE2` from which the quality scores for paired reads are pulled.

The `--degrade-variability` option can also be applied to introduce additional error.

This option is mutually exclusive to `--degrade-logistic`, `--degrade-phred`, and `--degrade-fastq-random`.

#### 5.2.5 Degradation specified by fastq quality scores, randomized

Specified by the `--degrade-fastq-random=FILE[,FILE2]` option.  This uses the error probabilities specified by reads in the fastq formatted file `FILE` and applies them to the generated reads.  This is done in a random fashion, such that each generated read uses the quality score of a random read from `FILE`.  The random order is without replacement, so fastq quality strings are not used more than once unless more reads are requested than there are lines in the fastq file.

For paired-end or amplicon reads one may specify a second file, `FILE2` from which the quality scores for paired reads are pulled.  These are randomized and assigned exclusively to the paired read.

The `--degrade-variability` option can also be applied to introduce additional error.

This option is mutually exclusive to `--degrade-logistic`, `--degrade-phred`, and `--degrade-fastq`.

#### 5.2.6 Degradation variability

Specificed by the `--degrade-variability=N` option.  `N` should be in range (0,1).  Default is 0 (off).  This introduces variability into underlying error specified by `--degrade-logistic`, `--degrade-phred`, `--degrade-fastq`, and `--degrade-fastq-random`.  `N` is interpreted as a maximum value by which the error rate should fluctuate at each position, relative to the value of the error rate at each position.  e.g. with N = 0.01, a position with an error rate `e` will have an effective error rate in the range `e +/- e * 0.01`.

Note that quality scores in degraded output reflect the additional variability introduced from this option on a per-position base.  Some positions may (randomly) have no extra error induced, whereas others may have up to N*100% additional error added (or subtracted).


### 5.3 Amplicon data

Amplicon sequencing uses a "priming" string to anchor reads near an area of interest, thus enriching the read depth at that location (versus standard single or pair-end sequencing).  This technique is commonly used in T-cell repertoire analysis to target the V/D/J recombination portion.

Amplicon probes are always given in a 5' --> 3' direction.  STIG will generate reads in the 3' direction from the amplicon probe, and supports reverse-strand matching.  The read length for amplicon sequencing can be specified by the --read-length parameters -- although note that in practice, most amplicon reads will be long enough to cover the entire mRNA from the probe location into the UTR.  In STIG's amplicon sequencing, a second "paired" read is generated which is simply the complement of the original strand.

#### 5.3.1 Forward strand example 

Let us say this is our forward strand (either DNA or RNA, the choice is not relevant for this example):

```
5' <--- AAAAAAAAAAAATTGTCCCCCCCCCCCCATAA ---> 3' 
```

If we call STIG with `--amplicon-probe=TTGT`, then this will "anchor" at the expected area and generate a read in the 3' direction (to the specified read-length).

```
5' <--- AAAAAAAAAAAATTGTCCCCCCCCCCCCATAA ---> 3' 
Primary read:       TTGTCCCCCCCCCCCCATAA...
Probe               ^^^^
```
Do not forget that we get a reverse complement read as well:

```
 Paired read:    ...TTATGGGGGGGGGGGGACAA
 Probe                              ^^^^

```
(n.b. STIG's output is always 5' --> 3', which is why the above string is reversed.)

#### 5.3.2 Reverse strand example

The C-region is a nice "constant" area to anchor probes in, but if the amplicon sequencing proceeds in a 3' direction this doesn't help capture the CDR3 region.  The solution is to use a probe which matches the reverse strand, so that reads capture the V, D, J portions 'upstream' (i.e. in the 5' direction) from the C-region probe.

Our DNA strand again, now with the complementary strand given:

```
5' <--- AAAAAAAAAAAATTGTCCCCCCCCCCCCATAA ---> 3'
3' <--- TTTTTTTTTTTTAACAGGGGGGGGGGGGTATT ---> 5'
```

If we call STIG with `--amplicon-probe=TTAT`, then this will anchor on the reverse strand.  The read will proceed in a 3' direction on the reverse strand, which means the reads come from the 5' direction on the forward strand:

```
5' <--- AAAAAAAAAAAATTGTCCCCCCCCCCCCATAA ---> 3'
3' <--- TTTTTTTTTTTTAACAGGGGGGGGGGGGTATT ---> 5'
Primary read: ...AAATTGTCCCCCCCCCCCCATAA
Probe                               ^^^^

Paired read:    TTATGGGGGGGGGGGGACAATTT...
Probe           ^^^^
```

(Again, the paired read is given 5' --> 3' )

#### 5.3.3 Amplicon probes

The default value for --amplicon-probe option is a 27-mer that anchors on the reverse strand of EX1 of the C-region on beta chains, about ~480nt from the V-region start codon: `GATCTCTGCTTCTGATGGCTCAAACAC`

Here is a slightly longer probe that anchors on the reverse strand of EX1 of the C-region on alpha chains, about ~510nt from the V-region start codon: `AGAATCCTTACTTTGTGACACATTTGTTTGAGA`

These values have only undergone limited testing in STIG, and may not function as intended in any wet-lab/in-vitro setting.


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

1. Loading a repertoire/population will attempt to re-access the chromosome files it was originally saved with.  There currently is not a way to change/update this, and modifying the chromosome files may result in unintended behavior when generating reads with a previously-saved repertoire
2. Loading a repertoire/population will contain the allele data from when the repertoire was generated.  As the alleles are used to generate the repertoires, this is intended behavior.
