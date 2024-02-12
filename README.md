STIG: Synthetic TCR Informatics Generator
=========================================

* Note original repository was at: https://github.com/Benjamin-Vincent-Lab/stig . Moved on 2/12/2024

INTRODUCTION
------------

STIG is a tool for creating artificial T-cell repertoires and producing simulated sequencing data from them.  Many characteristics of the repertoires can be customized (e.g. alpha/beta & gamma/delta distribution; allele substitution for V, D, J, and C regions).  Simulated sequencing data can be generated in DNA or RNA space, with a variety of customization options (e.g. paired vs. single-end reads, user-defined distribution of read & insert lengths).  Applications include evaluating and optimizing tools for performing analysis of T-cell receptors.

## Getting started
### Prerequisites

1. A working Python 3 installation.  This has been been tested with 3.5.2.
2. Two python package requirements: 'numpy' for numeric distributions required for making repertoires, and 'pyyaml' for unpacking the YAML-formatted TCR recombination data
3. Reference chromosomes placed in the working directory (`data`, by default).  See "Reference chromosomes" below.

#### Reference chromosomes

STIG is distributed with allele data from IMGT, as well as reference coordinates for that allele data compatible with hg38.  You will need a copy of the hg38 chromosome 7 and 14 reference files, which may be found at <http://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/>  Look for `chr7.fa.gz` and `chr14.fa.gz`, and unpack these into your working directory (defaults to `data`).  There's also a Makefile in the working directory included with STIG, running `cd data && make fetch` will fetch and unpack hg38 chromosomes 7 and 14.


### Installation

No installation is required. STIG can be run directly out of its extracted package directory, once the reference chromosomes are installed in the working directory.


### Running STIG

STIG can be called directly from its package directory.  Use `./lib/stig --help` to view the command-line help menu.

LATEST VERSION
--------------

The latest release can be found at <https://github.com/vincentlaboratories/stig>


DOCUMENTATION
-------------

Documentation for the most recent version is available on the project website.  A [copy of the manual](doc/manual.md) is included in each version as Markdown-formatted text.


COPYRIGHT AND LICENSE
-------
This software is copyright (C) 2019 The University of North Carolina at Chapel Hill.

STIG is free sofware licensed under the GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.

T-cell segment allele files are copyright IMGT(R), and for academic research only.  Please note that some of these files have been modified for the purposes of STIG: See the manual for details.
- IMGT(R), the international ImMunoGeneTics information system(R) http://www.imgt.org (founder and director: Marie-Paule Lefranc, Montpellier, France)

CONTACT
-------

* Maintainer: Mark Woodcock <mark.woodcock@unchealth.unc.edu>
* GitHub: https://github.com/vincentlaboratories/stig
