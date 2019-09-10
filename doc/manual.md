STIG: Synthetic TCR Informatics Generator
=========================================

Current for v0.5.4

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
10. COPYRIGHT 

## 1. NAME


STIG - Synthetic TCR Informatics Generator

## 2. SYNOPSIS


  ./lib/stig [options] working_dir

## 3. DESCRIPTION

STIG is a tool for creating artificial T-cell repertoires and producing simulated sequencing data from them.  Many characteristics of the repertoires and the sequencing output can be customized.  Reads can be generated in both RNA and DNA space.  Applications include evaluating and optimizing tools for performing analysis of T-cell receptors.

### 3.1 Conventions

STIG utilizes a number of underlying assumptions when reading input and output data:

1. Quality strings are Phred+33 (Illumina 1.8+) format: [0,41] or [!, J], with J being highest quality
2. DNA and RNA inputs and outputs are 5' --> 3' direction, even when generating paired end reads or specififying amplicon probes


## 4. OPTIONS

```
usage: stig [-h] [--output BASENAME] [--load-population FILE]
            [--repertoire-size N] [--repertoire-unique]
            [--repertoire-chain-unique] [--repertoire-cdr3-unique]
            [--population-size N]
            [--population-distribution {unimodal,chisquare,stripe,equal,logisticcdf}]
            [--population-unimodal-parameters N | --population-chisquare-parameters k:cutoff | --population-logisticcdf-parameters s:cutoff]
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
            WORKING_DIR

Generate synthetic TCR read data

positional arguments:
  WORKING_DIR           Directory with tcell_receptor.tsv,
                        tcell_recombination.yaml, reference chromosome(s), &
                        allele subdir. Try STIG's H. sapiens directory named
                        'data'

optional arguments:
  -h, --help            show this help message and exit
  --output BASENAME     Basename for output files, e.g. '--output=foo' will
                        write to 'foo.fastq', 'foo.statistics.csv', etc.
                        Default is 'stig.out'
  --load-population FILE
                        Load TCR population and repertoire data from FILE,
                        rather than generating from scratch
  --repertoire-size N   Size of the TCR repertoire (i.e. the number of unique
                        TCR clonotypes that are generated). Default is 10
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
  --population-size N   The approximate number of T-cells in the repertoire
                        (e.g. if repertoire-size=5 and population-size=15,
                        then there are, on average, 3 clones of each unique
                        TCR clonotype). Note that some population distribution
                        options may choose slightly fewer or more "cells"
                        depending on the particulars of the distribution.
                        Default is 100
  --population-distribution {unimodal,chisquare,stripe,equal,logisticcdf}
                        Population distribution function. This defines the
                        function used to distribute the population among the
                        repertoire. Default is the logistic CDF, approximating
                        a normalized distribution of TCR subclone population
                        sizes. 'stripe' will assign the Nth cell in the
                        population to the (N % repertoire-size) clonotype.
                        'equal' assigns cells in the population to each
                        clonotype with equal probability. 'unimodal' produces
                        a small set of clones with high population sizes
                        relative to the others. See --population-unimodal-
                        parameters, --population-chisquare-parameters,
                        --population-logisticcdf-parameters
  --population-unimodal-parameters N
                        Parameter for the unimodal population. The width of
                        the peak is defined by number of standard deviations
                        to include in our population distribution. Decimal
                        value. Default is 3
  --population-chisquare-parameters k:cutoff
                        Parameters for the chi-square distribution. Takes an
                        argument formatted as 'k:cutoff', where k - degrees of
                        freedom. Default is 3. cutoff - X-axis +/- maximum.
                        Default is 8
  --population-logisticcdf-parameters s:cutoff
                        Parameter for the logistic cumulative distribution
                        function. Takes an argument formatted as 's:cutoff',
                        where s - logistic scale. Default is 1. cutoff -
                        X-axis +/- maximum. Default is 3
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
                        Simulate non-optimal quality using a Phred+33
                        (Illumina 1.8+) string to specify quality on a per-
                        nucleotide basis. If a generated read is longer than
                        the given phred string, then the last character in the
                        phred string is used. Default is off. This option is
                        mutually exclusive to --degrade-logistic, --degrade-
                        fastq and --degrade-fastq-random. See: --degrade-
                        variability
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

The majority of options fall into two broad categories:

1. Generation of the T cell receptor repertoire.  Determines what types and how many T cell receptors are created
2. Data generation, which controls if reads are DNA or RNA space, and the read depth, read length, etc. of those reads.  Options also exist to 'degrade' the data in a manner that simlulates less-than-perfect read quality


## 5. INVOCATION

### 5.1 Quick usage

	./lib/stig --repertoire-size=100 --population-size 10000 --output=devel --sequence-count=0 ./data
Will create a STIG population of 100 subclones consisting of 10000 virtual cells in total, based on alleles and reference chromosomes in the working directory (`./data`, by default).  Population and statistics will be saved in `devel.population.bin` & `devel.statistics.tsv` files (`--output=devel`).  No output sequences are generated (`--sequence-count=0`).

	./lib/stig --load-population=devel.population.bin --read-type=paired --sequence-type=rna --sequence-count 50000000 ./data
Generate 50 million paired-end reads in RNA-space from the previously saved population(`--load-population=devel.population.bin`).  The default paired end options will define the characteristics of the paired reads (e.g. average read and insert length).

	./lib/stig --repertoire-size=100 --population-size 10000 --output=devel --read-type=paired --sequence-type=rna --sequence-count=50000000 ./data
This combines the two previous commands into a single step.



### 5.2 Degrading reads

By default, STIG generates outputs with 100% accuracy against the underlying simulated TCR.  These reads are placed in the output fastq file(s) (by default 'stig.out.fastq', but can be overridden by the `--output` option).  Because the Phred+33 standard does not have a "0% chance of error" score, STIG will label these reads as of the highest quality available ('J', by Phred+33 as implemented for Illumina 1.8+).

It is frequently useful to have less-than-perfect quality for testing purposes, as this simulates real-world data.  STIG has several mechanisms for degrading the output to match certain criteria.  Degraded outputs are given the suffix `.degraded` in the fastq output filename, and are generated in parallel with perfect-quality reads.  Each read is labeled with a `readnum` value that corresponds between degraded and non-degraded fastq files, so one can see exactly what changes were made to the output based on the quality string.

Note that just because a given position has a low quality score does not necessarily mean that the read doesn't match the underlying DNA/RNA.  STIG will randomly assign bases at each position based on the probabilities given by the quality string.  For example, a location with a degraded quality score of `I` will have an error probability of 0.001, or 0.1%.  If STIG determines a location should contain a simulated error, then the location is replaced with a random nucleotide.  Of course, there is a 25% chance it will substitute the same base value at that location.

#### 5.2.1 Displaying degredation

For tuning `--degrade-logistic` and `--degrade-variability` parameters, the `--display-degradation` option is available.  When given, STIG will not read input or generate outputs, but instead will display an example quality string and display the error rates per position.

#### 5.2.2 Degradation with a logistic function

Specified by the `--degrade-logistic=B:L:k:mid` option.  This uses a logistic (sigmoid) curve to degrade outputs based on their length.  The function parameters B:L:k:mid refer to the base error rate, maximum error rate, steepness factor, and midpoint length, respectively.  The `--display-degradation` option is often useful for tuning these options.

The `--degrade-variability` option can also be applied to introduce additional error.

This option is mutually exclusive to `--degrade-phred`, `--degrade-fastq` and `--degrade-fastq-random`.

#### 5.2.3 Degradation specified by a Phred+33 string

Specified by the `--degrade-phred=PHRED` option.  This uses the error probabilites speficied by the input string and applies them to the generated reads, where the error probability of the nth read in an output is given by the nth quality score of `PHRED`.  If the generated read is longer than `PHRED`, then the last quality score in `PHRED` is used for the remaining bases.  STIG uses the Illumina 1.8+ Phred+33 standard [i.e. each position must be one of: !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ ]

The `--degrade-variability` option can also be applied to introduce additional error, see below for more on this.

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

Specificed by the `--degrade-variability=N` option.  `N` should be in range (0,1).  Default is 0 (off).  This introduces variability into underlying error specified by `--degrade-logistic`, `--degrade-phred`, `--degrade-fastq`, and `--degrade-fastq-random`.  `N` is interpreted as a maximum value by which the error rate should fluctuate at each position, relative to the value of the error rate at each position.  e.g. with N = 0.01, a position with an error rate `e` will have an effective error rate in the range `[e - e * 0.01, e + e * 0.01]`.

Note that quality scores in degraded output reflect the additional variability introduced from this option on a per-position base.  Some positions may (randomly) have no extra error induced, whereas others may have up to N*100% additional error added or subtracted.


### 5.3 Amplicon data

Amplicon sequencing uses a "priming" string to anchor reads near an area of interest, thus enriching the read depth at that location (versus standard single or pair-end sequencing).  This technique is commonly used in T-cell repertoire analysis to target the V/D/J recombination portion.

Amplicon probes are always given in a 5' --> 3' direction.  STIG will generate reads in the 3' direction from the amplicon probe, and supports reverse-strand matching.  The read length for amplicon sequencing can be specified by the `--read-length parameters` -- although note that in practice, most amplicon reads will be long enough to cover the entire mRNA from the probe location into the UTR.  In STIG's amplicon sequencing, a second "paired" read is generated which is simply the complement of the original strand.

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

These values have only undergone limited testing within STIG, and may not function as intended outside of this setting.


### 5.4 Defining T-cell repertoires

( Terminology note: We use the terms *clonotype* or *subclone* to refer to a particular T-cell clone: a collection of cells all with identical receptors.  *Repertoire* refers to the set of all clonotypes within a population of cells.  When referring to individual cells, each will be an instance of a particular clonotype.  Of course, any two cells may be from the same, or different, clonotypes.)

There are two command line options that define the number of subclones and the number of cells shared between them: `--repertoire-size` defines the number of subclones, and `--population-size` defines the number of cells.

#### 5.4.1 Subclone distributions

STIG supports multiple options for distributing virtual cells between the underlying repertoire clonotypes, defined by the `--population-distribution` option.  The default value is `logisticcdf`.


##### logisticcdf 

This uses the logistic (or 'sigmoid') cumulative distribution function to distribute cells amongs the various clonotypes.  Thus, there will be a few subclones with very few cells (corresponding to the bottom tail of the sigmoid curve), there will be a larger collection of subclones with increasingly higher levels of clones (the middle of the curve), and few subclones with a very large number of cells.  The parameters of the logistic curve can be adjusted with the `--population-logisticcdf-parameters` argument, which can adjust the scale (or 'steepness') of the curve, as well as +/- X-axis cutoff (since the logistic is defined over an infinte range of negative and positive values).

n.b. That, in some cases, rounding losses when trying to fit the given number of cells into the given number of clonotypes can lead to cells that don't fit within the distribution.  In this case, STIG will first retry the fit multiple times by re-rolling the distribution with identical parameters.  In the event that it hits a maximum retry value (currently 500 times), STIG will instead distribute these losses by either subtracting clones from the left side of the distribution, or adding clones to the right side of the distribution.  Except for very small population values, this is unlikely to meaningfully impact the desired distribution.  STIG will print a warning message to alert the user of the adjustments for the rounding losses, to the effect of: `logisticcdf distribution encountered a rounding error when assigning X out of Y requested cells.  This may affect the distribution in rare cases, please verify acceptable subclone populations in your STIG population file. See the manual for more details`.

##### unimodal

Uses a gaussian-like CDF to distribute cells, thus creating a moderate number of subclones with lots of cells, and a tail of subclones with progressively fewer cells.  The width of this peak/tail can be adjusted with the `--population-unimodal-parameters` argument.

##### chisquare

Uses a chi-square CDF to distribute cells.  There are parameters to adjust the chi-square curve (`--population-chisquare-parameters`), which adjust the degrees of freedom and X-axis cutoff.  As the chi-square function has varied shapes over differing degrees of freedom, this allows one a lot of flexibility in generating population distributions that existing STIG options do not provide.

#### equal

Uses a flat function to distribte cells, thus providing equal odds of a given cell being assigned to each clonotype.  This approximates *stripe* below, but allows for some random variation between the sizes of each subclone.

#### stripe

This assigns cells by a round robin approach, where the Nth cell will belong to clonetype N % (repertoire-size).  If your `--population-size` is wholely divisible by your `--repertoire-size`, then each subclone will contain the same number of cells (e.g. 100 cells striped across 20 subclones will give 5 cells in each subclone).


## 6. SEE ALSO
* IMGT's overview of V(D)J recombination: http://www.imgt.org/IMGTeducation/Tutorials/index.php?article=IGandBcells&chapter=VariableRegion&lang=UK&nbr=article

## 7. FILES

### 7.1 The Working Directory

A directory which contains four key components:
1. `tcell_receptor.tsv`: A T-cell receptor component definition file
2. Some number of chromosome reference files, formatted as `chrN.fa` for chromosome N
3. `allele` directory: A subdirectory with FASTA files with IMGT-formatted headers which provides the nucleotide sequences of various T-cell receptor component alleles (e.g. the V, D, J alleles)
4. `tcell_recombination.yaml`: A YAML-formatted file with probabilities for gene segment recombination and chewback/nucleotide addition parameters

### 7.2 Input files
1. TCR FASTA files were downloaded from the ImMunoGeneTics GENE-DB site
(www.imgt.org/genedb).  The datasets were chosen to include "nucleotide 
sequences for F+ORF+all P".  As of this writing, a description of the available
datasets can be found at http://www.imgt.org/genedb/directlinks

    An example of these files can be found at:
http://www.imgt.org/genedb/GENElect?query=7.2+TRAV&species=Homo+sapiens

2. TCR segment descriptors (e.g. TRBV22-1 L & V-gene unit is on forward strand of 7q34 at nucleotides 142641746->142642235) were copied from the "Localization in genome assemblies" tables for the various TCR loci.  As of this writing, these can be found at http://www.imgt.org/genedb/

    An example of this table can be found at:
http://www.imgt.org/genedb/resultPage_localizations.action?mode=localizations&collectionId=4&regionId=1

3. Reference genome data for chromosomes 7 and 14 (assuming human TCRs).

#### 7.2.1 TCR Component Definitions

The T-cell receptor component file (`tcell_receptor.tsv` in the working_dir) specifies the (NCBI) coordinates of each component of the T-cell receptor that will be used to create a rearranged receptor.  In STIG, the D and J region formats are fairly straightforward, the C region is somewhat complicated, and the V region is the most complicated.

A good overview of how V(D)J recombination is done using IMGT-named blocks can be found at: http://www.imgt.org/IMGTeducation/Tutorials/index.php?article=IGandBcells&chapter=VariableRegion&lang=UK&nbr=article

##### V-Region Definition

Each V-region consists of three major components, defined by coordinates in the reference chromosome:

* L-V-GENE-UNIT: The entire gene, this includes the L-PART1, an intron, L-PART2, the V-region portion of the TCR, and finally some trailing bases that are trimmed during TCR rearrangement.  The L-PART1 contains the start codon (the stop codon spans the J-C junction).  This coordinate is needed so STIG can know where, in the reference chromosome, the V-region gene starts and thus where the 5' UTR data ends.
* V-REGION: The V-region itself.  This portion will contain the CYS residue that marks the start of the CDR3 sequence.  This coordinate is needed so STIG can splice in V-region alleles of varying lengths.
* L-PART1+L-PART2: The start and end coordinates of segment spanning L-PART1, the V-region intron, and L-PART2.  This coordinate is needed so STIG can splice in L-PART1+L-PART2 alleles of varying lengths when generating RNA reads.

##### C-Region Definition

The C-region is composed of a large segment with multiple (3 or 4, in humans) exons.  Because the C-region spans multiple exons, STIG requires the start and end coordinates of each exon, defined as segments EX1, EX2, EX3, and (sometimes) EX4.  These coordinates are used to splice in exon allele data.

###### D and J-Region Definition

The D and J regions don't have intronic data, and thus only a D-REGION and J-REGION coordinates are required so STIG can splice in alleles of varying length.


#### 7.2.2 Modifications to input files

N.b. that the versions of these files were hand-curated to exclude various logical inconsistencies that would interfere with the operation of this software.  These exclusions and edits fall into several broad categories:

1. V-REGIONs without L-V-GENE-UNIT regions
2. Orphons and pseudogenes
3. V-REGIONs without L-PART1+L-PART2 data (which is needed to generate RNA)
4. Addition of explicit L-PART1-L-PART2 reference coordinates
5. Addition of D-REGION reference coordinates
6. Correction of NCBI-format reference coordinates


## 8. AUTHORS

Initial code: Mark Woodcock, University of North Carolina at Chapel Hill
	mark.woodcock@unchealth.unc.edu
	


## 9. BUGS
### 9.1 Reporting Bugs

If you find a bug in STIG, please report it.  To ensure the problem exists in the most current version of STIG, you may download the latest version at https://github.com/vincentlaboratories/stig

Bug report can be filed under the project page on GitHub.


### 9.2 Known issues/limitations

1. Loading a repertoire/population will attempt to re-access the chromosome files it was originally saved with. There currently is not a way to change/update this, and modifying the chromosome files may result in unintended behavior when generating reads with a previously-saved repertoire

2. Loading a repertoire/population will contain the allele data from when the repertoire was generated (i.e. it will not attempt to re-load allele data from the working directory).  As the initial alleles are used to generate the TCR repertoires, this is intended behavior.

3. The handling of C-region alleles is semi-functional.  In RNA, an allele for each exon are randomly pulled when requested (e.g. Requesting EX1 of TRBC1*01 will return either EX1 of TRBC1*01 or EX1 of TRBC1*02).  In DNA, the alleles are not used at all (data is instead pulled from the chromosome file).  The complexity of code needed to splice in multiple exon alleles is nontrivial, and the C-regions aren't highly utilized in TCR reconstruction and analysis, in any event.  A fix can be implemented if there is sufficient demand.  A workaround for this is to only have a single allele for each C region, and ensure this allele matches the reference chromosomes: this ensures that all C-regions will have the same exonic sequences in DNA and RNA.

4. Similar to the above C-region alleles, there are no alleles for L-PART1, L-PART2, or the DNA-space of the V-region intron.  Thus, these nucleotides are pulled directly from the reference chromosome when simulating DNA sequencing data.  There's a fair bit of complexity in providing this functionality, and this region likely does not contribute much to TCR diversity or functionality.  A fix can be implemented if there is sufficient demand.

## 10. COPYRIGHT

Copyright (C) 2019 The University of North Carolina at Chapel Hill.  

License: GPLv3+: GNU  GPL  version  3  or  later http://gnu.org/licenses/gpl.html  This is free software: you are free to change and redistribute it.  There is NO WARRANTY, to the extent permitted by law.
