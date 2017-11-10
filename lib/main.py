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

import tcrFOOBAR


# Configure our logging
log = logging.getLogger('main')
log.setLevel(logging.DEBUG)

# File handler to log debug messages
#fh = logging.FileHandler('tcr.py.log');
#fh.setLevel(logging.DEBUG);
#log.addHandler(fh);

# Stream handler to log warning and higher messages
sh = logging.StreamHandler()
sh.setFormatter(logging.Formatter(fmt='%(asctime)s.%(msecs)03d [%(levelname)s] %(name)s %(message)s',
                                  datefmt='%Y%m%d%H%M%S'))
log.addHandler(sh);

# Seed our internal random number generator
random.seed()


parser = argparse.ArgumentParser(description="Generate TCR read data")
parser.add_argument('input',
										help="Read TCR data from one or more IMGT-formatted fasta files", metavar='INPUT', nargs='+')
parser.add_argument("--output",
										help="Name of output fastq file", required=True)
parser.add_argument('--repertoire-size', type=int, required=True,
										help='Size of the TCR repertoire (i.e. the number of unique TCRs that are generated)')
parser.add_argument('--population-size', type=int, required=True,
										help='The number of T-cells in the repertoire (e.g. if repertoire-sze=5 and population-size=15, then there are 3 clones of each unique TCR on average)')
parser.add_argument('--read-type', choices = ['paired', 'single'], default = 'single',
										help='Generate either single or paired-end reads')
parser.add_argument('--sequence-type', choices = ['dna', 'rna'], default = 'dna',
										help='Generate sequences from simulated DNA or RNA')
parser.add_argument("--sequences",
										help="Number of sequences to generate.  Default is 1000", metavar="SEQ", type=int, default=1000)
parser.add_argument("--receptor-ratio",
										help="Ratio of alpha/beta vs gamma/delta sequences.  Default is 0.9 (9 alpha/beta per 1 gamma/delta)", metavar="RATIO", type=float,default=0.9)
parser.add_argument('--log-level',
										help='Logging level.  Default is warning and above', choices=['debug', 'info', 'warning', 'error', 'critical'], default='warning')

args = parser.parse_args()


# Default values
receptor_ratio = 0.9
sequences_requested = 1000

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



log.debug("Processing FASTA inputs: %s", args.input)


if( args.sequences ):
  if( args.sequences <= 0 ):
    log.critical("--sequences requires argument > 0");
    exit(-1);
  sequences_requested = args.sequences
if (args.output):
  fq_out_name = args.output
else:
  log.critical("--output is a required argument");
  exit(-1);
  



my_configuration = tcrFOOBAR.tcrConfig(log=log.getChild('tcrConfig'))
my_configuration.readTCRConfig('./data/tcell_receptor.tsv')
my_configuration.readAlleles( args.input )
my_configuration.setChromosomeFiles(chr7='./data/chr7.fa', chr14='./data/chr14.fa')

log.debug("\n\n\n\n\n\n\n\n\n\n")
my_repertoire = tcrFOOBAR.tcrRepertoire(my_configuration, args.repertoire_size, AB_frequency=args.receptor_ratio, log=log.getChild('tcrRepertoire'))
my_repertoire.populate(args.population_size, 'gaussian')

# Write this repertoire out to a file
with open(fq_out_name, 'w') as fp:
		i = 0
		pairedEndReadOption = (True) if args.read_type == 'paired' else False
		for read in my_repertoire.simulateRead(args.sequences, args.sequence_type, read_length_mean=25, paired_end=pairedEndReadOption):
				qualStr = '~'*len(read)
				fp.write("@TCRSIM_%d\n" %(i))
				fp.write("%s\n" % (read))
				fp.write("+\n")
				fp.write("%s\n" % qualStr)
				i += 1

log.info("All actions complete")
exit(0)


###############################################################################
###############################################################################
###############################################################################
###############################################################################











































def generate_dna_tcr( sequences, receptor_ratio ):
  receptor_type=''
  V_sequence=''
  D_sequence=''
  J_sequence=''
  C_sequence=''
  # First decide if between alpha/beta and gamma/delta
  if random.random() < receptor_ratio:
    if random.random() < 0.5:
      receptor_type = 'alpha'
      
    else:
      receptor_type = 'beta'
  elif random.random() < 0.5:
    receptor_type = 'gamma'
  else:
    receptor_type = 'delta'
    
  log.info('Generating receptor type %s', receptor_type)
  return
    



# read_fasta( inputs )
# Desc: Read FASTA files and return their contents
#
# Args:
# inputs - Single filename (string), or list of filenames (list)
#
# Returns: Array with alternating header (">" line) and sequence elements
#          read from the input files provided
# Notes:
# - This function scoops all file contents in one read, if you need to read in
# large files you may not want to use this function.
# - Large REGEXP: This removes carriage-return+newline between sequences, needed 
# because althout fasta supports them, it is much easier to process if they are 
# removed.  (We also remove the trailing newline from the last line)

def read_fasta( inputs ):
  if not isinstance(inputs, list):
    inputs = [inputs]

  retval = []
  for input in inputs:
    line_num = 1
    content = re.sub('([ctag]{1})(\r?\n)+([ctag]{1}|\Z)(\n\Z)?', r'\1\3', open(input, 'rU').read()).split("\n")
    for line in content:
      if (
        ( (line_num % 2 == 1) and not re.search('^>.+$', line) ) or
        ( (line_num % 2 == 0) and not re.search('^[ctag]+$', line) )
        ):
        log.critical("In file %s, line %d does not appear to be FASTA formatted: \"%s\"",
                     input, line_num, line)
        exit(-1)
      else:
        retval.append(line)
        line_num += 1
  return retval
        


















###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################


# TODO: Review
def reverse_complement(line):
  revcomp_key = maketrans('ACGTRYMKBDHVacgtrymkbdhv', 'TGCAYRKMVHDBtgcayrkmvhdn') # the DNA complementation key
  line = line.rstrip() # remove new line charac2ter from end
  line = line[::-1]  #reverse sequence
  rev_comp_line = line.translate(revcomp_key) #complement sequence
  return rev_comp_line

# TODO: Review
def match_regex(line):
   #regular expression from abmining
   cdr3_pattern = r'(TT[TC]|TA[CT])(TT[CT]|TA[TC]|CA[TC]|GT[AGCT]|TGG)(TG[TC])(([GA][AGCT])|TC)[AGCT]([ACGT]{3}){5,32}TGGG[GCT][GCT]'
   match = re.search(cdr3_pattern, line)
   #match = re.search(cdr3_pattern, reverse_complement(line))
   return match
   #match = re.search(cdr3_pattern, reverse_complement(line)) #check the reverse complement first
   #if match:
   #  return match
   #else:
   #  match = re.search(cdr3_pattern, line)
   #  return match

# TODO: Review
def noStopCodons(sequence,frame):
  #return true if there are no stop codons
  stops = ['TAA','TAG','TGA']
  seqAfter = sequence[frame:]
  k = 0
  while k < len(seqAfter)-2:
    for s in stops:
      if seqAfter[k:k+3] == s:
        return False
    k+=3
  seqBefore = sequence[:frame]
  k = len(seqBefore)
  while k > 2:
    for s in stops:
      if seqBefore[k-3:k] == s:
        return False
    k = k - 3
  return True

# TODO: Review
def getConstant():
  return random.choice(constantSeqs);

# TODO: Review
def mutate(m):
  newN = m
  while (newN == m):
    newN = random.choice(nucleotides)
  return newN

# TODO: Review
def addMutations(s,rate):
  mutated = ""
  for c in s:
    rand = random.random();
    if rand < rate/100.0:
      mutated = mutated + mutate(c)
    else:
      mutated = mutated + c
  return mutated

# TODO: Review
def addIndels(s1,s2):
  #randomly delete 0-10 bases, then add 0-10 bases
  rand = random.random();
  if rand < indel/100.0:
    deletions = random.randint(0,10);
    insertions = random.randint(0,10);
    indelstart = random.randint(0,deletions);
    s1 = s1[:len(s1)-deletions + indelstart];
    s2 = s2[indelstart:];
    for i in range(0,insertions):
      s1 = s1 + random.choice(nucleotides);
  return s1 + s2

# TODO: Review
def writeToFastq(toWrite,header):
  fq_out.write(header);
  fq_out.write('\n');
  fq_out.write(toWrite);
  fq_out.write('\n');
  fq_out.write("+");
  fq_out.write('\n');
  qual = "";
  for c in range(len(toWrite)):
    qual += "C";
  fq_out.write(qual);
  fq_out.write('\n');





###########################################

log = logging.getLogger('main')
log.setLevel(logging.DEBUG);

# File handler to log debug messages
#fh = logging.FileHandler('tcr.py.log');
#fh.setLevel(logging.DEBUG);
#log.addHandler(fh);

# Stream handler to log warning and higher messages
sh = logging.StreamHandler()
#sh.setLevel(logging.WARNING);
sh.setLevel(logging.DEBUG);
sh.setFormatter(logging.Formatter(fmt='%(asctime)s.%(msecs)03d main [%(levelname)s] %(message)s',
                                  datefmt='%Y%m%d%H%M%S'))
log.addHandler(sh);


# TODO: Will clean this up later
parser = argparse.ArgumentParser()
parser.add_argument("--vdj_location", help="path to folder containing idhv.fa, idhd.fa, and idhj.fa")
parser.add_argument("--bed_location", help="path to all.bed")
parser.add_argument("--gtf_location", help="path to gtf")
parser.add_argument("--fa_location", help="path to folder containing chr*.fa")
parser.add_argument("--vdj_mutation", help="percent mutation chance for vdj region",type=float)
parser.add_argument("--constant_mutation", help="percent mutation chance for constant region",type=float)
parser.add_argument("--indel_percent", help="percent chance for a deletion of 0-10 bases followed by insertion of 0-10 bases at each junction",type=float)
parser.add_argument("--num_seq", help="number of vdj sequences to generate",type=int)
parser.add_argument("--constant", help="a single constant region to use")
parser.add_argument("--output", help="name of output fastq file")



args = parser.parse_args()

#vdj_path = "/datastore/nextgenout4/share/labs/bioinformatics/dmarron/simulate_vdj/original_vdj/";
#bed_path = "/datastore/nextgenout4/share/labs/bioinformatics/dmarron/simulate_vdj/original_vdj/all.bed";
#gtf_path = "/datastore/nextgenout4/share/labs/bioinformatics/dmarron/simulate_vdj/original_vdj/hg38.knowngene.chr14.gtf";
#fa_path = "/datastore/nextgenout4/share/labs/bioinformatics/dmarron/simulate_vdj/original_vdj/";
#fq_out_name = "/datastore/nextgenout4/share/labs/bioinformatics/dmarron/simulate_vdj/python/vdj_simulated.fastq"

# Default values for our command line arguments
vdj_path = ""
bed_path = "./data/all.bed"
gtf_path = "./data/hg38.knowngene.chr14.gtf"
fa_path = "./data"
fq_out_name = ""


if (args.vdj_location):
  log.info("Using %s as vdj_path", args.vdj_location)
  vdj_path = args.vdj_location;
if (args.bed_location):
  bed_path = args.bed_location;
if (args.gtf_location):
  gtf_path = gtf.bed_location;
if (args.fa_location):
  fa_path = args.fa_location;
if (args.constant):
  constant = args.constant;
if (args.output):
  fq_out_name = args.output;
else:
  log.critical("--output is a required argument");
  exit(-1);
  

# Igheavy regions from these files
v_in = vdj_path + "ighv.fa"
d_in = vdj_path + "ighd.fa"
j_in = vdj_path + "ighj.fa"


log.debug("opening files");

nucleotides = ['A','T','C','G'];
#constants = ['IGHM','IGHD','IGHG1','IGHA1'];
#constants = ['IGHG1'];
constants = [constant];

# Data files
# BED: chromosomes, start & end loci for various IGH genes
# GTF: Chromosomes, start & end loci, classification (exon, start_codon, end_codon, CDS), and gene and transcript ID using UC naming
bed_file = open(bed_path, 'rU');
bed_list = list(bed_file);
gtf_file = open(gtf_path, 'rU');
gtf_list = list(gtf_file);
fq_out = open(fq_out_name, 'w');



# Grab only lines that contain 'exon' from out gtf file
gtf_exons = filter(lambda e:'exon' in e, gtf_list);

# Loop description
# 1. For each given constant region C:
# 2. For each occurance of C in the BED file:
#    (Although there should only be one, it appears)
#    - Open the chromosome file associated with this C in the BED file
#    - Read all the sequences in...
# 3. For each GTF line matching current chromosome, e:
#    - Grab the start, end loci, as well as the UC transcript ID from GTF file
# TODO: Trim these comments

cchrom = "chr14"
cstart = 0;
cend = 0;
strand = "";
constantSeqs = [];
for c in constants:
  # Grab only lines from bed list with our constant c in it
  for b in filter(lambda x:c in x, bed_list):
    if b.split()[3] == c:
      log.info("Pulling details for constant %s from .bed file", c);
      cchrom = b.split()[0];
      cstart = int(b.split()[1]);
      cend = int(b.split()[2]);
      strand = b.split()[5];
      log.info("- Chromosome: %s", cchrom);
      log.info("- Start: %s", cstart);
      log.info("- End: %s", cend);
      log.info("- Strand: %s", strand);
      #extract constant sequences from fa file
      fa_file = open(fa_path + "/" + cchrom + ".fa", 'rU'); #FYI this is always chr14 with our current data
      fa_list = list(ifilterfalse(lambda line: '>' in line,fa_file));

      faLineLength = len(fa_list[1].rstrip());
      log.debug("Setting faLineLength based on line 1: %d", faLineLength);
      #check splices between cstart and cend
      prevend = 0;
      prevtrans = "";
      spliceStarts = [];
      spliceEnds = [];
      foundTranscript = False;
      log.info("Looking for transcripts from chromosome %s in our gtf exons", cchrom);
      for e in filter(lambda ch:cchrom in ch, gtf_exons):
        estart = int(e.split("\t")[3]);
        eend = int(e.split("\t")[4]);
        trans = e.split("\t")[8].split("\"")[3];
        log.info("Examining transcript %s [%s-%s]", trans, estart, eend);
        if (trans == prevtrans):
          #splice found in gtf file
          if (prevend > cstart and prevend < cend):            
            log.info("- Splice detection triggered by previous transcript!\n                                              (DNA coordinates [%s-%s]\n                                              (splice coords   [%s-%s]", cstart, cend, prevend, estart);
            spliceStarts.append(prevend);
            spliceEnds.append(estart);
            foundTranscript = True;
        else:
          if foundTranscript:
            foundTranscript = False;
            break;
        prevtrans = trans;
        prevend = eend;


      nstart = cstart;
      nend = cend;
      constantSequence = "";
      log.info("Processing our splice list N=%d...", len(spliceStarts));
      for s in range(0,len(spliceStarts)):
        
        nend = spliceStarts[s];
        startLineNum = nstart/faLineLength;
        startBase = nstart - faLineLength * startLineNum;
        endLineNum = nend/faLineLength;
        endBase = nend - faLineLength * endLineNum;
        
        constantSequence += fa_list[startLineNum][startBase:].rstrip();
        for i in range(startLineNum+1,endLineNum):
          constantSequence += fa_list[i].rstrip();
        constantSequence += fa_list[endLineNum][:endBase].rstrip();
        nstart = spliceEnds[s];
      nend = cend;
      startLineNum = nstart/faLineLength;
      startBase = nstart - faLineLength * startLineNum
      endLineNum = nend/faLineLength;
      endBase = nend - faLineLength * endLineNum;
      constantSequence += fa_list[startLineNum][startBase:].rstrip();
      for i in range(startLineNum+1,endLineNum):
        constantSequence += fa_list[i].rstrip();
      constantSequence += fa_list[endLineNum][:endBase].rstrip();
      if strand == "-":
        revSequence = reverse_complement(constantSequence);
        constantSeqs.append(revSequence);
      else:
        constantSeqs.append(constantSequence);
      fa_file.close();

log.info("We dropped down here...");
exit(0);
r1 = 1
r2 = 0.2
if (args.vdj_mutation):
  r1 = args.vdj_mutation
if (args.constant_mutation):
  r2 = args.constant_mutation
indel = 50
if (args.indel_percent):
  indel = args.indel_percent
v_file = open(v_in, 'rU')
#list only sequences from fa and filter out headers
#v_list = list(ifilterfalse(lambda line: '>' in line,v_file))
v_lines = v_file.readlines();
v_list,v_headers = [],[]
for line in v_lines:
  if '>' in line:
    v_headers.append(line);
  else:
    v_list.append(line);
d_file = open(d_in, 'rU')
#d_list = list(ifilterfalse(lambda line: '>' in line,d_file))
d_lines = d_file.readlines();
d_list,d_headers = [],[]
for line in d_lines:
  if '>' in line:
    d_headers.append(line);
  else:
    d_list.append(line);
j_file = open(j_in, 'rU')
#j_list = list(ifilterfalse(lambda line: '>' in line,j_file))
j_lines = j_file.readlines();
j_list,j_headers = [],[]
for line in j_lines:
  if '>' in line:
    j_headers.append(line);
  else:
    j_list.append(line);
sequenceCount = 1000
if (args.num_seq):
  sequenceCount = args.num_seq;
numSequences = 0
numIterations = 0
regexRejected = 0
stopRejected = 0
while (numSequences < sequenceCount and numIterations < 100*sequenceCount):
  #randomV = random.choice(v_list);
  #randomD = random.choice(d_list);
  #randomJ = random.choice(j_list);
  rand = random.randint(0, len(v_list) - 1);
  randomV = v_list[rand];
  vHeader = v_headers[rand];
  rand = random.randint(0, len(d_list) - 1);
  randomD = d_list[rand];
  dHeader = d_headers[rand];
  rand = random.randint(0, len(j_list) - 1);
  randomJ = j_list[rand];
  jHeader = j_headers[rand];
  vdjun = addIndels(randomV.rstrip(),randomD.rstrip());
  combined = addIndels(vdjun,randomJ.rstrip());
  #combined = randomV.rstrip() + randomD.rstrip() + randomJ.rstrip();
  constant = getConstant();
  combined = addMutations(combined,r1)
  constant = addMutations(constant,r2)
  finalSeq = combined + constant[:400]
  m = match_regex(finalSeq)
  if (m) :
    if (noStopCodons(finalSeq,m.start())) :
    #if (noStopCodons(combined,m.start())):
      finalHeader = "@VDJSIM_" + str(numSequences) + "_" + vHeader.split("|")[1] + "_" + dHeader.split("|")[1] + "_" + jHeader.split("|")[1];
      writeToFastq(finalSeq,finalHeader)
      numSequences+=1
    else:
      stopRejected+=1
  else:
    regexRejected+=1
  numIterations+=1
#for y in constantSeqs:
#  print y[:10]
print 'total iterations',numIterations
print 'rejected by regex',regexRejected
print 'rejected by stop codon',stopRejected
v_file.close();
d_file.close();
j_file.close();
bed_file.close();
fq_out.close();
