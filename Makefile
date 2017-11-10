

TCR_OPTS=\
--output=outputs.fastq \
--repertoire-size=5000 \
--population-size=100000 \
--read-type=single \
--sequence-type=dna \
--sequences=500 ./data/*.fasta \
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
	./lib/main.py --log-level=debug $(TCR_OPTS)
