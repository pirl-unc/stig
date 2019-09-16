all:
	@echo "Valid make targets:";
	@echo "fetch - Retrieve and unpack hg38 reference chromosomes 7 and 14";

fetch:
	wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr7.fa.gz
	wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr14.fa.gz
	gunzip chr7.fa.gz
	gunzip chr14.fa.gz
