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


parser = argparse.ArgumentParser(description="Generate synthetic TCR read data")
parser.add_argument('input', metavar='INPUT', nargs='*',
										help="Read TCR data from one or more IMGT-formatted fasta files")
parser.add_argument("--output", default='tcr_synth',
										help="Name of output fastq file.  This should be a basename, e.g. '--output=foo' will write to 'foo.fastq'.  Default is 'tcr_synth'")
parser.add_argument('--repertoire-size', type=int, default=10,
										help='Size of the TCR repertoire (i.e. the number of unique TCRs that are generated).  Default is 10')
parser.add_argument('--population-size', type=int, default=100,
										help='The number of T-cells in the repertoire (e.g. if repertoire-sze=5 and population-size=15, then there are 3 clones of each unique TCR on average).  Default is 100')
parser.add_argument('--read-type', choices = ['paired', 'single'], default = 'single',
										help='Generate either single or paired-end reads')
parser.add_argument('--sequence-type', choices = ['dna', 'rna'], default = 'dna',
										help='Generate sequences from simulated DNA or RNA')
parser.add_argument("--sequence-count",
										help="Number of sequences (reads) to generate.  Default is 1000", metavar="SEQ", type=int, default=1000)
parser.add_argument("--read-length-mean", type=int, default=48,
										help="The average length of reads in nucleotides. Default is 48")
parser.add_argument("--read-length-sd", type=int, default=4,
										help="The SD of read length variation in nucleotides. Set to zero for fixed-length reads.  Default is 4")
parser.add_argument("--read-length-sd-cutoff", type=int, default=4, metavar='N',
										help="Read lengths are restricted to less than N standard deviations from the mean.  Default is 4")
parser.add_argument("--insert-length-mean", type=int, default=48,
										help="The average length of the insert for paired end reads.  Default is 48")
parser.add_argument("--insert-length-sd", type=int, default=4,
										help="The SD of insert length variation in nucleotides. Set to zero for fixed-length inserts.  Default is 4")
parser.add_argument("--insert-length-sd-cutoff", type=int, default=4, metavar='N',
										help="Insert lengths are restricted to less than N standard dviations from the mean.  Default is 4")

parser.add_argument("--degrade", default=None, metavar='B:L:k:mid',
										help="Simulate non-optimatal quality.  Value is a colon-separated list to define a logistic function outlining the error rate per base position.  B - Base error rate probability.  L - Maximum error rate. k - Steepness factor. mid - Midpoint, this is the base position where error rate is equal to 1/2 of L")
parser.add_argument("--display-degradation",  metavar='B:L:k:mid', default=None,
										help="Display the error rate per base pair for a given B:L:k:mid value and exit.  The number of positions displayed is adjustable through the --read-length-mean option.  This is mostly useful in adjusting these parameters to be passed to the --degrade option.  Note that no reads or repertoire will be generated when this option is given")

parser.add_argument("--receptor-ratio", metavar="RATIO", type=float, default=0.9,
										help="Ratio of alpha/beta vs gamma/delta sequences.  Default is 0.9 (9 alpha/beta per 1 gamma/delta TCR)")
parser.add_argument('--log-level', choices=['debug', 'info', 'warning', 'error', 'critical'], default='warning',
										help='Logging level.  Default is warning and above')

# Process our command-line arguments
args = parser.parse_args()

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



# Process the 'display-degedation' option
if args.display_degradation is not None:
		components = re.split(':', args.display_degradation)
		if( components == None or
				not isinstance(components, list) or
				len(components) != 4 or
				components[0] <= 0 or components[1] <= 0 or components[2] <= 0 or components[3] <= 0):
				log.error("Invalid string for degradation \"%s\".  Valid example: 0.005:0.2:0.25:15")
				exit(-1)
		baseError, L, k, midpoint = components
		for i in range(0, args.read_length_mean):
				errorRate = (float(L) - float(baseError)) / (1 + math.exp(-float(k)*(i - int(midpoint)))) + float(baseError)
				print("Position %-2d, error rate: %0.4f" % (i, errorRate))
		exit(0)
elif args.input is None:
		print("Error: Must give either input FASTA files, or --display-degradation option")
		exit(-1)


		
my_configuration = tcrFOOBAR.tcrConfig(log=log.getChild('tcrConfig'))
my_configuration.readTCRConfig('./data/tcell_receptor.tsv')
my_configuration.readAlleles( args.input )
my_configuration.setChromosomeFiles(chr7='./data/chr7.fa', chr14='./data/chr14.fa')

my_repertoire = tcrFOOBAR.tcrRepertoire(my_configuration, args.repertoire_size, AB_frequency=args.receptor_ratio, log=log.getChild('tcrRepertoire'))
my_repertoire.populate(args.population_size, 'gaussian')


pairedReadOpt = True if args.read_type == 'paired' else False
outputSequences = my_repertoire.simulateRead(args.sequence_count, args.sequence_type,
																						 read_length_mean      = args.read_length_mean,
																						 read_length_sd        = args.read_length_sd,
																						 read_length_sd_cutoff = args.read_length_sd_cutoff,
																						 inner_mate_length_mean      = args.insert_length_mean,
																						 inner_mate_length_sd        = args.insert_length_sd,
																						 inner_mate_length_sd_cutoff = args.insert_length_sd_cutoff,
																						 paired_end=pairedReadOpt)
# Write the read sequences to output file(s)
if args.read_type == 'single':
		outputFilename = args.output + '.fastq'
		with open(outputFilename, 'w') as fp:
				i = 0
				for read in outputSequences:
						qualStr = 'I'*len(read)
						fp.write("@TCRSIM_%d\n" %(i))
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

								output1.write("@TCRSIM_%d\n" %(i))
								output1.write("%s\n" % (read1))
								output1.write("+\n")
								output1.write("%s\n" % qualStr1)

								output2.write("@TCRSIM_%d\n" %(i))
								output2.write("%s\n" % (read2))
								output2.write("+\n")
								output2.write("%s\n" % qualStr2)

								i += 1
else:
		raise ValueError("Unknown read_type encountered " + args.read_type)


# Write our degraded-quality reads, if requested by the user
if args.degrade is not None:
		components = re.split(':', args.degrade)
		if( components == None or
				not isinstance(components, list) or
				len(components) != 4 or
				components[0] <= 0 or components[1] <= 0 or components[2] <= 0 or components[3] <= 0):
				log.error("Invalid string for degradation \"%s\".  Valid example: 0.005:0.2:0.25:15")
				exit(-1)
		baseError, L, k, midpoint = components

		if args.read_type == 'single':
				outputFilename = args.output + '.degraded.fastq'
				with open(outputFilename, 'w') as fp:
						i = 0
						for read in outputSequences:
								ident = "@TCRIM_%d" % (i)
								fp.write(my_configuration.getDegradedFastq(read, float(baseError), float(L), float(k), int(midpoint), ident))
								i += 1
		elif args.read_type == 'paired':
				output1Filename = args.output + '1.degraded.fastq'
				output2Filename = args.output + '2.degraded.fastq'
				with open(output1Filename, 'w') as output1:
						with open(output2Filename, 'w') as output2:
								i = 0
								for readPair in outputSequences:
										read1, read2 = readPair
										ident = "@TCRIM_%d" % (i)
										output1.write(my_configuration.getDegradedFastq(read1, float(baseError), float(L), float(k), int(midpoint), ident))
										output2.write(my_configuration.getDegradedFastq(read2, float(baseError), float(L), float(k), int(midpoint), ident))
										i += 1

				
		else:
				raise ValueError("Unknown read_type encountered" + args.read_type)






log.info("All actions complete")
exit(0)


