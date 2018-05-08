

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


DEVEL_OPTS=\
--chr7-filename=./data/chr7.fa \
--chr14-filename=./data/chr14.fa \
--tcell-data=./data/tcell_receptor.tsv \
--output=test_output \
--repertoire-size=10 \
--population-size=1000 \
--read-type=paired \
--sequence-type=rna \
--sequence-count=5 \
--read-length-mean=48 \
--read-length-sd=0 \
--degrade-phred='555555555555776869382849@@@AAABBCCC' \
--degrade-variability=0.25 \
--log-level=debug \
data/*.fasta

TEST_OPTS=\
--chr7-filename=./data/chr7.fa \
--chr14-filename=./data/chr14.fa \
--tcell-data=./data/tcell_receptor.tsv \
--log-level=debug \
data/*.fasta

TR_BIN=/usr/bin/tr
GREP_BIN=/bin/grep


all:
	@echo "Valid targets:";
	@echo "work  - Main workflow";
	@echo "devel - Development workflow";

work:
	./lib/stig $(TCR_OPTS)

devel:
	./lib/stig $(DEVEL_OPTS)



# Test targets
test: test_help test_display_degradation
	./lib/test.py $(TEST_OPTS)

test_help:
	./lib/stig --help

test_display_degradation:
	./lib/stig --degrade-phred='555555555555' --degrade-variability=0.5 --display-degradation

