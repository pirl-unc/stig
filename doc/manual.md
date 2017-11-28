TODO: Pick a better title for this program
==========================================

## Table of contents
1.
NAME
SYNOPSIS
DESCRIPTION
OPTIONS
INVOCATION
SEE ALSO
FILES
AUTHORS
BUG REPORTS

# 1. NAME
TODO - Synthetic read generator for T-cell receptor data

# 2. SYNOPSIS

  ./lib/main.py [options] [allele fasta files]

# 3. DESCRIPTION

TODO

# 4. OPTIONS


# 5. INVOCATION

## 5.1 Quick usage

To 

# 2. DATA FILES

## 2.1 Input files
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

