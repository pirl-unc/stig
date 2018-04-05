#! /usr/bin/python


# TCR Synth - Generate synthetic T-cell receptor reads in DNA and RNA
# v0.0.1
# Mark Woodcock (mark.woodcock@unchealth.unc.edu)
# Based on makeVDJfastq by David Marron (dmarron@email.unc.edu)


import sys
import re
import random
import argparse
from string import maketrans
from itertools import ifilter
from itertools import ifilterfalse
import logging
import math
import cPickle

import tcrFOOBAR


# Configure our logging
log = logging.getLogger('main')
log.setLevel(logging.DEBUG)

# Stream handler to log warning and higher messages
sh = logging.StreamHandler()
sh.setFormatter(logging.Formatter(fmt='%(asctime)s.%(msecs)03d [%(levelname)s] %(name)s %(message)s',
                                  datefmt='%Y%m%d%H%M%S'))
log.addHandler(sh);

# Seed our internal random number generator
random.seed()


parser = argparse.ArgumentParser(description = "Generate synthetic TCR read data",
																 epilog = "Please see manual or README for further details" )

parser.add_argument('allele_files', metavar='ALLELE_FILES', nargs='*',
										help='Read TCR allele data from one or more IMGT-formatted fasta files')

parser.add_argument('--chr7-filename', metavar='FILE',
										help = "Filename of a FASTA formatted file with chromosome 7 reference data")
parser.add_argument('--chr14-filename', metavar='FILE',
										help = "Filename of a FASTA formatted file with chromosome 7 reference data")

parser.add_argument('--tcell-data', metavar='FILE',
										help = "Filename of a tab-separated file with T cell receptor segments and reference gene coordinates")

parser.add_argument("--output", metavar='BASENAME', default='tcr_synth',
										help='Basename for output files, e.g. \'--output=foo\' will write to \'foo.fastq\', \'foo.statistics.csv\', etc.  Default is \'tcr_synth\'')

parser.add_argument("--load-population", metavar='FILE', type=str,
										help='Load TCR population and repertoire data from FILE, rather than generating from scratch')

parser.add_argument('--repertoire-size', metavar='N', type=int, default=10,
										help='Size of the TCR repertoire (i.e. the number of unique TCRs that are generated).  Default is 10')
parser.add_argument('--population-size', metavar='N', type=int, default=100,
										help='The number of T-cells in the repertoire (e.g. if repertoire-sze=5 and population-size=15, then there are 3 clones of each unique TCR on average).  Default is 100')
parser.add_argument('--population-distribution', choices = ['gaussian', 'chisquare'], default='gaussian',
										help = 'Population distribution function.  This defines the function used to distribute the population among the repertoire.  Default is the (normalized) gaussian.  See --population-gaussian-parameters, --population-chisquare-parameters')

parserGroup1 = parser.add_mutually_exclusive_group()
parserGroup1.add_argument('--population-gaussian-parameters', metavar='N', type=float, default=3.0,
													help='Parameter for the normalized gaussian distribution.  The number of standard deviations to include in our population distribution.  Decimal value.  Default is 3')
parserGroup1.add_argument('--population-chisquare-parameters', metavar='k:cutoff', default='2:8',
													help='Parameters for the chi-square distribution.  Takes an argument formatted as \'k:cutoff\', where k - degrees of freedom.  Default is 3. cutoff - X-axis maximum.  Default is 8')

parser.add_argument('--read-type', choices = ['paired', 'single'], default = 'single',
										help='Generate either single or paired-end reads.  Default is single')
parser.add_argument('--sequence-type', choices = ['dna', 'rna'], default = 'dna',
										help='Generate sequences from simulated DNA or RNA. Default is DNA')
parser.add_argument("--sequence-count", metavar="N", type=int, default=1000,
										help='Number of sequences (reads) to generate.  Default is 1000')
parser.add_argument("--read-length-mean", type=int, default=48,
										help='The average length of reads in nucleotides. Default is 48')
parser.add_argument("--read-length-sd", type=int, default=4,
										help='The SD of read length variation in nucleotides. Set to zero for fixed-length reads.  Default is 4')
parser.add_argument("--read-length-sd-cutoff", type=int, default=4, metavar='N',
										help='Read lengths are restricted to less than N standard deviations from the mean.  Default is 4')
parser.add_argument("--insert-length-mean", type=int, default=48,
										help='The average length of the insert for paired end reads.  Default is 48')
parser.add_argument("--insert-length-sd", type=int, default=4,
										help='The standard deviation of insert length variation in nucleotides. Set to zero for fixed-length inserts.  Default is 4')
parser.add_argument("--insert-length-sd-cutoff", type=int, default=4, metavar='N',
										help='Insert lengths are restricted to less than N standard dviations from the mean.  Default is 4')

parserGroup2 = parser.add_mutually_exclusive_group()
parserGroup2.add_argument("--degrade-logistic", default=None, metavar="B:L:k:mid",
										help='Simulate non-optimal quality using the logstic (sigmoid) function.  Takes an argument formatted as \'B:L:k:mid\'.  B - Base error rate probability.  L - Maximum error rate. k - Steepness factor. mid - Midpoint, this is the base position where error rate is equal to 1/2 of L. Default is off.  This option is mutually exclusive to --degrade-phred.  See: --display-degradation, --degrade-variability')
parserGroup2.add_argument("--degrade-phred", metavar="PHRED_STRING", default=None,
										help='Simulate non-optimal quality using a Phred+33 string to specify quality on a per-nucleotide basis.  If a generated read is longer than the given phred string, then the last character in the phred string is used.  Default is off.  This option is mutually exclusive to --degrade-logistic.  See: --degrade-variability')

parser.add_argument("--degrade-variability", default=0, metavar='FLOAT', type=float,
										help='Applies a relative variability in the per-nucleotide error applied by the --degrade option.  If a given base were to have an error rate of 0.1 (10%%), then a degrade-variability of 0.5 (50%%) would result in an error rate in the range of 0.1 +/- 0.1 * 0.5.  Default is 0')

parser.add_argument("--display-degradation", action = 'store_true',
										help='Display the error rate per base pair for a given B:L:k:mid value and exit.  The number of positions displayed is adjustable through the --read-length-mean option.  This is mostly useful in adjusting these parameters to be passed to the --degrade option.  Note that no reads or repertoire will be generated when this option is given')

parser.add_argument("--receptor-ratio", metavar="RATIO", type=float, default=0.9,
										help='Ratio of alpha/beta vs gamma/delta sequences.  Default is 0.9 (9 alpha/beta per 1 gamma/delta TCR)')
parser.add_argument('--log-level', choices=['debug', 'info', 'warning', 'error', 'critical'], default='warning',
										help='Logging level.  Default is warning and above')

args = parser.parse_args()

# Process our logging level arguments
if( args.log_level == 'debug'):
  sh.setLevel(logging.DEBUG)
elif( args.log_level =='info' ):
  sh.setLevel(logging.INFO)
elif( args.log_level =='warning' ):
  sh.setLevel(logging.WARNING)
elif( args.log_level =='error' ): 
  sh.setLevel(logging.ERROR)
elif( args.log_level =='critical' ):
  sh.setLevel(logging.CRITICAL)
else:
  log.error("Error: Unknown log level %s", args.log_level)


# Process options 'display-degradation', 'degrade-logistic' and 'degrade-phred'
degradeOptions = None
if args.degrade_logistic is not None or args.degrade_phred is not None :
		method, baseError, L, k, midpoint, phred = [0] * 6 # Initialize to zero
		
		if args.degrade_logistic is not None:
				if re.match('^((?:\d+)|(?:\d*.\d+)):((?:\d+)|(?:\d*.\d+)):((?:\d+)|(?:\d*.\d+)):((?:\d+)|(?:\d*.\d+))$', args.degrade_logistic):
						log.info("Using logistic function for degradation")
						method = 'logistic'
						baseError, L, k, midpoint = args.degrade.split(':')
				else:
						log.critical("Invalid string for --degrade-logistic: \"%s\".  Valid example: 0.005:0.2:0.25:15", args.degrade_logistic)
						exit(-1)
		elif args.degrade_phred is not None:
				if re.match(r'^[!\"#\$%&\'\(\)\*\+,-./0123456789:;<=>?@ABCDEFGHI]+$', args.degrade_phred):
						log.info("Using Phred string for degradation")
						method = 'phred'
						matches = re.match('^(.+)$', args.degrade_phred)
						phred = matches.groups()[0]
				else:
						log.critical("Invalid argument for --degrade-phred: \"%s\".  Valid example: IIIIIIII444433", args.degrade_phred)
						exit(-1)
		else:
				raise ValueError("Fallen through to an invalid location")

		degradeOptions = {
				'method': method,
				'baseError': baseError,
				'L': L,
				'k': k,
				'midpoint': midpoint,
				'phred': phred
				}

# Display degradation output, if --display-degradation given
if( args.display_degradation is True and
		degradeOptions is not None ):
		displayString = "A" * args.read_length_mean
		tempConfig = tcrFOOBAR.tcrConfig()
		tempConfig.getDegradedFastq(displayString, method, 'ident',  variability=args.degrade_variability,
																phred=degradeOptions['phred'],
																baseError=degradeOptions['baseError'], L=degradeOptions['L'],
																k=degradeOptions['k'], midpoint=degradeOptions['midpoint'],
																display=True)
		exit(0)
elif args.display_degradation is True:
		raise ValueError("--display-degradation requires a degradation method.  See --degrade-logistic, --degrade-phred under help")


#  Make sure we have the necessary options given to Generate synthetic reads
if( args.chr7_filename  is None or
		args.chr14_filename is None or
		args.tcell_data     is None or
		args.allele_files   is None ):
		log.critical("Required command-line options missing (chr7-filename, chr14-filename, tcell-data or allele files).  See --help")
		exit(-10)

# Create our configuration object
my_configuration = tcrFOOBAR.tcrConfig(log=log.getChild('tcrConfig'))
my_configuration.readTCRConfig(args.tcell_data)
my_configuration.readAlleles( args.allele_files )
my_configuration.setChromosomeFiles(chr7=args.chr7_filename, chr14=args.chr14_filename)





##################				
#
log.info("Performing segment and allele sequence testing")

vGene = filter(lambda x:x['gene'] == 'TRAV20' and
								 x['region'] == 'L-V-GENE-UNIT', my_configuration.receptorSegment)
if len(vGene) == 0:
		log.error("Failed to identify L-V-GENE-UNIT for TRAV20")
		exit(-10)
elif len(vGene) > 1:
		log.error("Duplicate L-V-GENE-UNIT locations for TRAV20")
		exit(-10)
vGeneSeq = my_configuration.readChromosome(14, vGene[0]['start_position'], vGene[0]['end_position'], vGene[0]['strand'])
trav20GeneSeq = 'ATGGAGAAAATGTTGGAGTGTGCATTCATAGTCTTGTGGCTTCAGCTTGGCTGTAAGTTGGAGGTTAAGAAATGAGAACCATTCGTGCAAAGCTGAGGAAGAGAGATGAAGATTTAATCTCAGCTCATTCAGTTTTCCTGTTGGAATTTTTCAAAGCTGCTAAACTAAAAGCTTTGTATAATTTCTGATTTTGTTTTCCTGAACAGGGTTGAGTGGAGAAGACCAGGTGACGCAGAGTCCCGAGGCCCTGAGACTCCAGGAGGGAGAGAGTAGCAGTCTCAACTGCAGTTACACAGTCAGCGGTTTAAGAGGGCTGTTCTGGTATAGGCAAGATCCTGGGAAAGGCCCTGAATTCCTCTTCACCCTGTATTCAGCTGGGGAAGAAAAGGAGAAAGAAAGGCTAAAAGCCACATTAACAAAGAAGGAAAGCTTTCTGCACATCACAGCCCCTAAACCTGAAGACTCAGCCACTTATCTCTGTGCTGTGCAGGCACAGCGTTCCCCAGGCACCTGCAACTTGTATCAAAACC'
if vGeneSeq != trav20GeneSeq:
		log.error("Unexpected L-V-GENE region (expected %s, got %s)", trav20GeneSeq, vGeneSeq)
		exit(-10)

vRegion = filter(lambda x:x['gene'] == 'TRAV20' and
								 x['region'] == 'V-REGION', my_configuration.receptorSegment)
if len(vRegion) == 0:
		log.error("Failed to identify V-REGION for TRAV20")
		exit(-10)
elif len(vRegion) > 1:
		log.error("Duplicate V-REGION locations for TRAV20")
		exit(-10)
vRegionIndex = my_configuration.receptorSegment.index(vRegion[0])
vAlleleSeq = my_configuration.getSegmentSequences((vRegionIndex, '01'))
trav20LVGeneSeq='atggagaaaatgttggagtgtgcattcatagtcttgtggcttcagcttggctggttgagtggagaagaccaggtgacgcagagtcccgaggccctgagactccaggagggagagagtagcagtcttaactgcagttacacagtcagcggtttaagagggctgttctggtataggcaagatcctgggaaaggccctgaattcctcttcaccctgtattcagctggggaagaaaaggagaaagaaaggctaaaagccacattaacaaagaaggaaagctttctgcacatcacagcccctaaacctgaagactcagccacttatctctgtgctgtgcagg'.upper()

if vAlleleSeq[1] != trav20LVGeneSeq:
		log.error("Unexpected V-region RNA (expected %s, got %s)", trav20LVGeneSeq, vAlleleSeq[1])
		exit(-10)


vRegion = filter(lambda x:x['gene'] == 'TRAJ22' and
								 x['region'] == 'J-REGION', my_configuration.receptorSegment)
if len(vRegion) == 0:
		log.error("Failed to identify J-REGION for TRAJ22")
		exit(-10)
elif len(vRegion) > 1:
		log.error("Duplicate J-REGION locations for TRAJ22")
		exit(-10)
vRegionIndex = my_configuration.receptorSegment.index(vRegion[0])
traj22AlleleSeq = my_configuration.getSegmentSequences((vRegionIndex, '01'))
targetSeq='TTTCTTCTGGTTCTGCAAGGCAACTGACCTTTGGATCTGGGACACAATTGACTGTTTTACCTG'

if traj22AlleleSeq[1] != targetSeq:
		log.error("Unexpected J-region RNA (expected \n%s, got \n%s)", targetSeq, traj22AlleleSeq[1])
		exit(-10)

VJSeq = trav20LVGeneSeq + traj22AlleleSeq[1]
VJCSeq = 'ATGGAGAAAATGTTGGAGTGTGCATTCATAGTCTTGTGGCTTCAGCTTGGCTGGTTGAGTGGAGAAGACCAGGTGACGCAGAGTCCCGAGGCCCTGAGACTCCAGGAGGGAGAGAGTAGCAGTCTCAACTGCAGTTACACAGTCAGCGGTTTAAGAGGGCTGTTCTGGTATAGGCAAGATCCTGGGAAAGGCCCTGAATTCCTCTTCACCCTGTATTCAGCTGGGGAAGAAAAGGAGAAAGAAAGGCTAAAAGCCACATTAACAAAGAAGGAAAGCTTTCTGCACATCACAGCCCCTAAACCTGAAGACTCAGCCACTTATCTCTGTGCTTAGTCGTTTCTTCTGGTTCTGCAAGGCAACTGACCTTTGGATCTGGGACACAATTGACTGTTTTACCTGTCTGTCTGCCTATTCACCGATTTTGATTCTCAAACAAATGTGTCACAAAGTAAGGATTCTGATGTGTATATCACAGACAAAACTGTGCTAGACATGAGGTCTATGGACTTCAAGAGCAACAGTGCTGTGGCCTGGAGCAACAAATCTGACTTTGCATGTGCAAACGCCTTCAACAACAGCATTATTCCAGAAGACACCTTCTTCCCCAGCCCAGAAAGTTCCTGTGATGTCAAGCTGGTCGAGAAAAGCTTTGAAACAGATACGAACCTAAACTTTCAAAACCTGTCAGTGATTGGGTTCCGAATCCTCCTCCTGAAAGTGGCCGGGTTTAATCTGCTCATGACGCTGCGGCTGTGGTCCAGCTGA'
log.info("%s", VJSeq)
prettySeq = []
for i in range(0, len(VJSeq), 3):
		prettySeq.append(VJSeq[0+i:3+i])
log.info("%s", ','.join(prettySeq))



vGene = filter(lambda x:x['gene'] == 'TRBV20-1' and
								 x['region'] == 'L-V-GENE-UNIT', my_configuration.receptorSegment)
if len(vGene) == 0:
		log.error("Failed to identify L-V-GENE-UNIT for TRBV20-1")
		exit(-10)
elif len(vGene) > 1:
		log.error("Duplicate L-V-GENE-UNIT locations for TRBV20-1")
		exit(-10)
extractedSeq = my_configuration.readChromosome(7, vGene[0]['start_position'], vGene[0]['end_position'], vGene[0]['strand'])
targetSeq = 'ATGCTGCTGCTTCTGCTGCTTCTGGGGCCAGGTATAAGCCTCCTTCTACCTGGGAGCTTGGGTGGGCATGTGCGTGTGTTGGCATGGTCAAGTGGTGGCCAGCAGGGTTGCAATGTGGATTGTTTATGCTCATCGAAGGGAGAGGGAGAGGCCCTGCTCTCTAGAGGTGTAAATGGTAAGGTGAAAGCCGCCGGTCAGAGGAGATGGGGGATTATGGCCCTAGGGAGATGACGGGAAGATTGCACAAAACAAACAGGACTCTCCAGGAGCTGGGAGCACAGGGAGGGAGTGAGGCTCAGCTCTGCCTGGCGTCCTGTCTGACTCGGCTCCCACTGGGCTCTCCTCTCTCTCTGGCTTCTGTCTCAGCAGGCTCCGGGCTTGGTGCTGTCGTCTCTCAACATCCGAGCAGGGTTATCTGTAAGAGTGGAACCTCTGTGAAGATCGAGTGCCGTTCCCTGGACTTTCAGGCCACAACTATGTTTTGGTATCGTCAGTTCCCGAAACAGAGTCTCATGCTGATGGCAACTTCCAATGAGGGCTCCAAGGCCACATACGAGCAAGGCGTCGAGAAGGACAAGTTTCTCATCAACCATGCAAGCCTGACCTTGTCCACTCTGACAGTGACCAGTGCCCATCCTGAAGACAGCAGCTTCTACATCTGCAGTGCTAGAGACACAGCGCCAGGAGGGGATCAGACACCGCGGCAAGAACC'

if extractedSeq != targetSeq:
		log.error("Unexpected L-V-GENE region (expected %s, got %s)", targetSeq, extractedSeq)
		exit(-10)

vRegion = filter(lambda x:x['gene'] == 'TRBJ2-1' and
								 x['region'] == 'J-REGION', my_configuration.receptorSegment)
if len(vRegion) == 0:
		log.error("Failed to identify J-REGION for TRBJ2-1")
		exit(-10)
elif len(vRegion) > 1:
		log.error("Duplicate J-REGION locations for TRBJ2-1")
		exit(-10)
jRegionIndex = my_configuration.receptorSegment.index(vRegion[0])
vAlleleSeq = my_configuration.getSegmentSequences((jRegionIndex, '01'))
targetSeq=['atgctgctgcttctgctgcttctggggccaggtataagcctccttctacctgggagcttggcaggctccgggcttggtgctgtcgtctctcaacatccgagctgggttatctgtaagagtggaacctctgtgaagatcgagtgccgttccctggactttcaggccacaactatgttttggtatcgtcagttcccgaaacagagtctcatgctgatggcaacttccaatgagggctccaaggccacatacgagcaaggcgtcgagaaggacaagtttctcatcaaccatgcaagcctgaccttgtccactctgacagtgaccagtgcccatcctgaagacagcagcttctacatctgcagtgctagaga'.upper(), 'ATGCTGCTGCTTCTGCTGCTTCTGGGGCCAGCAGGCTCCGGGCTTGGTGCTGTCGTCTCTCAACATCCGAGCTGGGTTATCTGTAAGAGTGGAACCTCTGTGAAGATCGAGTGCCGTTCCCTGGACTTTCAGGCCACAACTATGTTTTGGTATCGTCAGTTCCCGAAACAGAGTCTCATGCTGATGGCAACTTCCAATGAGGGCTCCAAGGCCACATACGAGCAAGGCGTCGAGAAGGACAAGTTTCTCATCAACCATGCAAGCCTGACCTTGTCCACTCTGACAGTGACCAGTGCCCATCCTGAAGACAGCAGCTTCTACATCTGCAGTGCTAGAGA']

if not targetSeq.index(vAlleleSeq[1]):
		log.error("Unexpected V-region RNA (expected \n%s, got \n%s)", targetSeq, jAlleleSeq[1])
		exit(-10)

seq = ##
prettySeq = []
for i in range(0, len(seq), 3):
		prettySeq.append(seq[0+i:3+i])
log.info("%s", ','.join(prettySeq))

quit()



##################				
# TEST
# Verify out 5' UTR and 3' UTR reads
my_repertoire = tcrFOOBAR.tcrRepertoire(my_configuration, 1, log=log.getChild('tcrRepertoire'))
my_repertoire.populate(1, 'gaussian')
outputSequences = my_repertoire.simulateRead(1, 'rna',
																						 read_length_mean = 100,
																						 read_length_sd   = 0,
																						 read_length_sd_cutoff = 1,
																						 paired_end = False )

for read in outputSequences:
		log.info("Read: %s", read)

quit()

# Obtain our simulated reads
pairedReadOpt = True if args.read_type == 'paired' else False
outputSequences = my_repertoire.simulateRead(args.sequence_count, args.sequence_type,
																						 read_length_mean      = args.read_length_mean,
																						 read_length_sd        = args.read_length_sd,
																						 read_length_sd_cutoff = args.read_length_sd_cutoff,
																						 inner_mate_length_mean      = args.insert_length_mean,
																						 inner_mate_length_sd        = args.insert_length_sd,
																						 inner_mate_length_sd_cutoff = args.insert_length_sd_cutoff,
																						 paired_end = pairedReadOpt)

quit()

# Write the read sequences to output file(s)
if args.read_type == 'single':
		outputFilename = args.output + '.fastq'
		with open(outputFilename, 'w') as fp:
				i = 0
				for read in outputSequences:
						qualStr = 'I'*len(read)
						fp.write("@TCR_SIMULATED_%d\n" %(i))
						fp.write("%s\n" % (read))
						fp.write("+\n")
						fp.write("%s\n" % qualStr)
						i += 1
elif args.read_type == 'paired':
		output1Filename = args.output + '1.fastq'
		output2Filename = args.output + '2.fastq'
		with open(output1Filename, 'w') as output1:
				with open(output2Filename, 'w') as output2:
						i = 0
						for readPair in outputSequences:
								read1, read2 = readPair
								qualStr1 = 'I'*len(read1)
								qualStr2 = 'I'*len(read2)

								output1.write("@TCR_SIMULATED_%d\n" %(i))
								output1.write("%s\n" % (read1))
								output1.write("+\n")
								output1.write("%s\n" % qualStr1)

								output2.write("@TCR_SIMULATED_%d\n" %(i))
								output2.write("%s\n" % (read2))
								output2.write("+\n")
								output2.write("%s\n" % qualStr2)

								i += 1
else:
		raise ValueError("Unknown read_type encountered " + args.read_type)


# Write degraded-quality reads, if requested by the user
# n.b. the cmd line options were parsed previously and placed in degradeOptions dict
if degradeOptions is not None:

		method = degradeOptions['method']
		baseError = float(degradeOptions['baseError'])
		L = float(degradeOptions['L'])
		k = float(degradeOptions['k'])
		midpoint = float(degradeOptions['midpoint'])
		phred = degradeOptions['phred']
		
		if args.read_type == 'single':
				outputFilename = args.output + '.degraded.fastq'
				with open(outputFilename, 'w') as fp:
						i = 0
						for read in outputSequences:
								ident = "@TCR_SIMULATED_%d" % (i)
								fp.write(my_configuration.getDegradedFastq(read, method, ident, variability=args.degrade_variability, phred=phred, baseError=baseError, L=L, k=k, midpoint=midpoint))
								i += 1
		elif args.read_type == 'paired':
				output1Filename = args.output + '1.degraded.fastq'
				output2Filename = args.output + '2.degraded.fastq'
				with open(output1Filename, 'w') as output1:
						with open(output2Filename, 'w') as output2:
								i = 0
								for readPair in outputSequences:
										read1, read2 = readPair
										ident = "@TCR_SIMULATED_%d" % (i)
										output1.write(my_configuration.getDegradedFastq(read1, method, ident, variability=args.degrade_variability, phred=phred, baseError=baseError, L=L, k=k, midpoint=midpoint))
										output2.write(my_configuration.getDegradedFastq(read2, method, ident, variability=args.degrade_variability, phred=phred, baseError=baseError, L=L, k=k, midpoint=midpoint))
										i += 1

				
		else:
				raise ValueError("Unknown read_type encountered" + args.read_type)



# Write statistics to an output file
statsFilename = args.output + '.statistics.csv'
with open(statsFilename, 'w') as fp:
		fp.write("CELL_COUNT, CDR3_1, RNA_1, DNA_1, CDR3_2, RNA_2, DNA_2\n")
		for i in my_repertoire.getStatistics():
				fp.write(", ".join(str(e) for e in i) + "\n")

# Write our repertoire object to a file
populationFilename = args.output + '.population.bin'
with open(populationFilename, 'wb') as fp:
		cPickle.dump(my_repertoire.freeze(), fp)


log.info("All actions complete")
exit(0)


