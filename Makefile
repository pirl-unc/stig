BASE_OPTS=\
--chr7-filename=./data/chr7.fa \
--chr14-filename=./data/chr14.fa \
--tcell-data=./data/tcell_receptor.tsv \
data/*.fasta

LOAD_DEVEL_OPTS=\
--load-population=devel.population.bin

TCR_OPTS=\
--chr7-filename=./data/chr7.fa \
--chr14-filename=./data/chr14.fa \
--tcell-data=./data/tcell_receptor.tsv \
--output=test_output \
--repertoire-size=100 \
--population-size=10000 \
--read-type=single \
--sequence-type=dna \
--sequence-count=500 \
--read-length-mean=80 \
--read-length-sd=0 \
data/*.fasta


AMPLICON_OPTS=\
--read-type=amplicon \
--sequence-type=rna \
--sequence-count=5 \
--read-length-mean=48 \
--read-length-sd=0

DEGRADE_LOGISTIC_OPTS=\
--degrade-logistic=

DEGRADE_PHRED_OPTS=\
--degrade-phred='5555555IJ5555555555555' \
--degrade-variability=0.5

DEGRADE_DUAL_FASTQ_OPTS=\
--degrade-fastq=./data/test_degrade.fastq,./data/test_degrade.fastq
DEGRADE_DUAL_FASTQRAND_OPTS=\
--degrade-fastq=./data/test_degrade.fastq,./data/test_degrade.fastq


TR_BIN=/usr/bin/tr
GREP_BIN=/bin/grep


all:
	@echo "Valid targets:";
	@echo "work  - Main workflow";
	@echo "devel - Development workflow";

work:
	./lib/stig $(TCR_OPTS)

devel: distribution_debug

degrade_fastq: devel.population.bin
	./lib/stig --log-level=debug --degrade-fastq=./data/test_degrade.fastq $(LOAD_DEVEL_OPTS) $(BASE_OPTS)
	./lib/stig --degrade-fastq-rand=./data/test_degrade.fastq $(LOAD_DEVEL_OPTS) $(BASE_OPTS)
	./lib/stig --read-type=paired --degrade-fastq=./data/test_degrade.fastq,./data/test_degrade.fastq $(LOAD_DEVEL_OPTS) $(BASE_OPTS)
	./lib/stig --read-type=amplicon --degrade-fastq=./data/test_degrade.fastq,./data/test_degrade.fastq $(LOAD_DEVEL_OPTS) $(BASE_OPTS)
	./lib/stig --read-type=paired --degrade-fastq-rand=./data/test_degrade.fastq,./data/test_degrade.fastq $(LOAD_DEVEL_OPTS) $(BASE_OPTS)
	./lib/stig --read-type=amplicon --degrade-fastq-rand=./data/test_degrade.fastq,./data/test_degrade.fastq $(LOAD_DEVEL_OPTS) $(BASE_OPTS)

degrade_phred: devel.population.bin
	./lib/stig $(LOAD_DEVEL_OPTS) $(DEGRADE_PHRED_OPTS) $(BASE_OPTS)

devel.population.bin:
	./lib/stig --output=devel --repertoire-size=100 --population-size=100000 --sequence-count=0 $(BASE_OPTS)

amplicon: devel.population.bin
	./lib/stig --amplicon-probe=GATCTCTGCTTCTGATGGCTCAAACAC --output=devel --load-population=devel.population.bin $(DEGRADE_PHRED_OPTS) $(AMPLICON_OPTS) $(BASE_OPTS)
	./lib/stig --amplicon-probe=AGAATCCTTACTTTGTGACACATTTGTTTGAGA --output=devel --load-population=devel.population.bin $(DEGRADE_PHRED_OPTS) $(AMPLICON_OPTS) $(BASE_OPTS)

distribution_debug:
	./lib/stig --population-distribution=logisticcdf --output=distribution --repertoire-size=40 --population-size=1000 --sequence-count=0 $(BASE_OPTS)


# Test targets
test: test_help test_degrade test_amplicon
test_degrade: display_degradation degrade_fastq
test_amplicon: amplicon

test_help:
	./lib/stig --help 2>&1 >/dev/null

display_degradation:
	./lib/stig --degrade-phred='555555555555' --degrade-variability=0.5 --display-degradation > /dev/null

clean:
	rm -f stig.out*.fastq devel*.fastq devel.statistics.csv devel.population.bin
