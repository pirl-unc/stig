STIG: Synthetic TCR Informatics Generator
=========================================


INTRODUCTION
------------

STIG is a tool for creating artificial T-cell repertoires and producing simulated sequencing data from them.  Many characteristics of the repertoires can be customized (e.g. alpha/beta & gamma/delta distribution; allele substitution for V, D, J, and C regions).  Simulated sequencing data can be generated in DNA or RNA space, with a variety of customization options (e.g. paired vs. single-end reads, user-defined distribution of read & insert lengths).  Applications include evaluating and optimizing tools for performing analysis of T-cell receptors.

## Getting started
### Prerequisites

1. A working Python 2 installation.  This has been been tested with 2.7.12 and 2.7.5.
2. Two python package requirements: 'numpy' for numeric distributions required for making repertoires, and 'pyyaml' for unpacking the YAML-formatted TCR recombination data
3. Reference chromosomes described below.

#### Reference chromosomes

STIG is distributed with allele data from IMGT, as well as reference coordinates for that allele data compatible with hg38.  You will need a copy of the hg38 chromosome 7 and 14 reference files, which may be found at <http://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/>  Look for chr7.fa.gz and chr14.fa.gz, and unpack these into the working directory (defaults to ./data).

### Installation

No installation is required and STIG can be run directly out of its extracted package directory.


### Running STIG

STIG can be called directly from its package directory.  Use `./lib/stig --help` to view the command-line help menu.

LATEST VERSION
--------------

The latest release can be found at <https://github.com/vincentlaboratories/stig>


DOCUMENTATION
-------------

Documentation for the most recent version is available on the project website.  A [copy of the manual](doc/manual.md) is included in each version as Markdown-formatted text.


LICENSE
-------
STIG is licensed for non-commercial research purposes only - see LICENSE.txt for full license

T-cell segment allele files are copyright IMGT(R), and for academic research only.  Please note that some of these files have been modified for the purposes of STIG: See the manual for details.
- IMGT(R), the international ImMunoGeneTics information system(R) http://www.imgt.org (founder and director: Marie-Paule Lefranc, Montpellier, France)

CONTACT
-------

* Maintainer: Mark Woodcock <mark.woodcock@unchealth.unc.edu>
* GitHub: https://github.com/vincentlaboratories/stig
