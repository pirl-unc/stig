

TCR_OPTS=\
--output=test_output \
--repertoire-size=10 \
--population-size=100 \
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
--population-size=100 \
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
--log-level=warning \
data/*.fasta

TR_BIN=/usr/bin/tr
GREP_BIN=/bin/grep


all:
	@echo "Valid targets:";
	@echo "work  - Main workflow";
	@echo "devel - Development workflow";

work:
	./lib/main.py $(TCR_OPTS)

devel:
	./lib/main.py $(DEVEL_OPTS)



# Test targets
test: test_help test_display_degradation test_chisquare

test_help:
	./lib/main.py --help

test_display_degradation:
	./lib/main.py --degrade-phred='555555555555' --degrade-variability=0.5 --display-degradation

test_chisquare:
	./lib/main.py --population-distribution=chisquare $(TEST_OPTS)
