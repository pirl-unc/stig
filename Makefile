

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
--output=test_output \
--repertoire-size=1 \
--population-size=10 \
--read-type=paired \
--sequence-type=rna \
--sequence-count=5 \
--read-length-mean=48 \
--read-length-sd=0 \
--degrade=0.01:0.2:0.15:45 \
--log-level=debug \
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
