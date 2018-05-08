BASE_OPTS=\
--chr7-filename=./data/chr7.fa \
--chr14-filename=./data/chr14.fa \
--tcell-data=./data/tcell_receptor.tsv \
data/*.fasta



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

DEGRADE_OPTS=\
--degrade-phred='55555555555555555555' \
--degrade-variability=0.5

TR_BIN=/usr/bin/tr
GREP_BIN=/bin/grep


all:
	@echo "Valid targets:";
	@echo "work  - Main workflow";
	@echo "devel - Development workflow";

work:
	./lib/stig $(TCR_OPTS)

devel: amplicon

devel.population.bin:
	./lib/stig --output=devel --repertoire-size=100 --population-size=100000 --sequence-count=0 $(BASE_OPTS)

amplicon: devel.population.bin
	./lib/stig --amplicon-probe=GATCTCTGCTTCTGATGGCTCAAACAC --log-level=debug --output=devel --load-population=devel.population.bin $(DEGRADE_OPTS) $(AMPLICON_OPTS) $(BASE_OPTS)
	./lib/stig --amplicon-probe=AGAATCCTTACTTTGTGACACATTTGTTTGAGA --log-level=debug --output=devel --load-population=devel.population.bin $(DEGRADE_OPTS) $(AMPLICON_OPTS) $(BASE_OPTS)

# Test targets
test: test_help test_display_degradation
	./lib/test.py $(TEST_OPTS)

test_help:
	./lib/main.py --help

test_display_degradation:
	./lib/main.py --degrade-phred='555555555555' --degrade-variability=0.5 --display-degradation

