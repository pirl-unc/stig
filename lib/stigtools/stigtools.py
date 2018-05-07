import re
import string
import logging
import random
import math
import numpy
import time


# TCR configuration class
#
# Self variables:
# log - A logging object.  Can be passed as an option to our constructor via
#       log=<someobject>.  If a pre-configured logging object is not passed
#       then a non-functioning logging object will be created and used
#       throughout this class, effectively disabling logging.
# 
# receptorSegment - An array of dicts containing entries for each of our TCR
#   components.
# = [
#     {
#       gene: String.  The full name of the receptor segment defined here.
#         (e.g. TRGC2)
#
#       region: The type of information stored at this segment (e.g.
#         J-GENE-UNIT, J-REGION, V-REGION, EX1)
#
#       segment_type: One of 'V', 'D', 'J', or 'C'
#
#       segment_number: String, the segment number.  (e.g. '13-1' for TRAV13-1)
#
#       receptor_type: TCR type.  One of 'A', 'B', 'G', 'D', representing
#         alpha, beta, gamma, and delta receptors, respectively.
#
#       start_position:  Integer, coordinates in NCBI format indicating the
#         start position of the region of this TCR segment
#
#       end_position:  Integer, coordinates in NCBI format
#
#       chromosome: String, e.g. '14q11.2'
#
#       strand: String.  Strand orientation.  One of 'forward' or 'reverse'
#
#       allele: Dictionary, containing allele_name:sequence key:value pairs
#         n.b. the sequence does not have to be the same length as denoted
#         by the (start_position, end_position) coordinates
#          e.g. {
#                 '01': 'gctactgatctcaattacagt',
#                 '02': 'gctactgatctcaattacagttatt', ...
#               }
#     }, ...
#   ]
#



class tcrConfig:



		# Define our probabilities for VDJ recombination likelihoods.
		# The numbers for beta-chain likelihoods are adapted from figures in the work of J Freeman and R Warren as described below
		# See: JDFreeman, RLWarren, et al "Profiling the T-cell receptor beta-chain repertoire by massively parallel sequencing",
		#      Genome Research, 2009. https://doi.org/10.1101/gr.092924.109
		#
		# Examples
		#
		# 15% chance of choosing TRAV1-1 when selecting an alpha chain:
		#				( 'TRAV1-1', 0.15 ),
		# 
		#
		# 25% chance of choosing TRAJ22 if TRAV20 was the chosen V-segment for an alpha chain: 
		#				('TRAV20', 'TRAJ22', 0.25),
		#
		# Other examples for gamma and delta regions:
		#				( 'TRGV1', 0.30 ),
		#				( 'TRGV4', 0.18 ),
		#				( 'TRGV8', 0.17 ),
		#				( 'TRDV1', 0.21 ),

		VDJprobability = [
				( 'TRBV20-1', 0.24 ),
				( 'TRBV5-1', 0.125 ),
				( 'TRBV29-1', 0.105 ),
				( 'TRBV28', 0.051 ),
				( 'TRBV10-3', 0.045 ),
				( 'TRBV4-2', 0.04 ),
				
				('TRBV20-1', 'TRBJ2-1', 0.17 ),
				('TRBV20-1', 'TRBJ1-1', 0.12 ),
				('TRBV20-1', 'TRBJ2-7', 0.11 ),
				('TRBV20-1', 'TRBJ1-5', 0.10 ),
				('TRBV20-1', 'TRBJ2-3', 0.09 ),
				('TRBV20-1', 'TRBJ2-2', 0.08 ),

		]

		# Define our probabilities for chewback and NT addition at VDJ junction sites
		# Probabilities are defined by index of the number of nucleotides chewed or added
		# e.g. Vchewback[3] is the probability that 3 nucleotides will be chewed from V segments
		#
		# Vchewback   - Chewback for V region
		# VJadditions - Nucleotide additions during VJ fusion (alpha/gamma chains)
		# D5chewback  - Chewback for D region adjacent to V region
		# D3chewback  - Chewback for D region adjacent to J region
		# DJadditions - Nucleotide additions during DJ fusion (beta/delta chains)
		# Jchewback   - Chewback for J region
		#
		# Any/all credit for these numbers go to J Freeman and R Warren for their analysis of the beta-chain TCR
		# See: JDFreeman, RLWarren, et al "Profiling the T-cell receptor beta-chain repertoire by massively parallel sequencing",
		#      Genome Research, 2009. https://doi.org/10.1101/gr.092924.109
		#
		junctionProbability = {
				'Vchewback' : [ 0.2384, 0.1629, 0.1008, 0.1159, 0.1359, 0.1044, 0.0743, 0.0391, 0.0152, 0.0061, 0.0043, 0.0022, 0.0004, 0.000044 ],
				'VJaddition': [ 0.1625, 0.1179, 0.1466, 0.1403, 0.1330, 0.0817, 0.0635, 0.048, 0.0357, 0.0233, 0.0184, 0.0073, 0.0073, 0.0071, 0.0033, 0.0013, 0.0009, 0.0005, 0.0002, 0.0004, 0.0002, 0, 0.0009 ],
				'VDaddition': [ 0.1625, 0.1179, 0.1466, 0.1403, 0.1330, 0.0817, 0.0635, 0.048, 0.0357, 0.0233, 0.0184, 0.0073, 0.0073, 0.0071, 0.0033, 0.0013, 0.0009, 0.0005, 0.0002, 0.0004, 0.0002, 0, 0.0009 ],
				'D5chewback': [ 0.3231, 0.1683, 0.1299, 0.0746, 0.1077, 0.0555, 0.0504, 0.0444, 0.0462 ],
				'D3chewback': [ 0.2401, 0.1621, 0.2431, 0.2072, 0.0631, 0.0355, 0.0202, 0.0133, 0.0153 ],
				'DJaddition': [ 0.1246, 0.1035, 0.1457, 0.1466, 0.1257, 0.0975, 0.0793, 0.0389, 0.0369, 0.0227, 0.0133, 0.0164, 0.0076, 0.0135, 0.0053, 0.0060, 0.0045, 0.0045, 0.0029, 0.0020, 0.0009, 0.0005, 0.0002 ],
				'Jchewback' : [ 0.1197, 0.0867, 0.1037, 0.0930, 0.1253, 0.1143, 0.0897, 0.0734, 0.0548, 0.0327, 0.0207, 0.0131, 0.0101, 0.0095, 0.0072, 0.0139, 0.0059, 0.0077, 0.0065, 0.0032, 0.0027, 0.0019, 0.0013, 0.0008, 0.0005 ],
		}
				
		def __init__( self, log=None ):
				# Initialize our instance variables
				self.receptorSegment = []
				self.geneName = []
				self.chr7Filename = None
				self.chr7Offset = 0
				self.chr7LineLength = 0
				self.chr14Filename = None
				self.chr14Offset = 0
				self.chr14LineLength = 0

				self.setLog(log)

				return

		# setLog - Configure our logging object
		# 
		# args:
		# log - If a logging object, we will use this for our logging
		#       If None, we will configure a new, non-functioning logging object
		#       
		# Returns:
		#  nothing
		# 
		def setLog( self, log ):
				if( isinstance(log, logging.Logger) ):
						self.log = log
				elif log is None:
						self.log = logging.getLogger(__name__)
						self.log.setLevel(99) # A high level, effectively disabling logging
				else:
						raise ValueError("Log object for tcrConfig must be a logging.Logger (or None)")


				
		# rmLog - Remove our logging object
		#         This is usually called prior to serializing this object, as 
		#         logger objects cannot be serialized. see setLog()
		# 
		# Arguments: none
		# Returns: nothing
		# 
		def rmLog( self ):
				self.log = None

				
		def readTCRConfig( self, filename ):
				self.log.debug("Processing config file %s", filename)
				open( filename, 'rU' )
				line_number = 0
				for line in list(open( filename, 'rU' )):
						line_number += 1
						new_segment = {}
						
						line = re.sub('^(.+)#.*$', r'\1', line ) # Strip off trailing comments
						if ( re.match('^\s*$', line) or re.match('^#.*$', line) ):
								self.log.debug("Ignoring comment/blank line")
								continue

						if( len(line.split("\t")) == 15 ):
								component, chromosome, strand, x, x, x, x, x, region, x, x, x, x, coordinates, x = line.split("\t")
						else:
								raise ValueError("Unexpected line #%din file %s" % (line_number, filename))
						

						# Validate the TCR component field (e.g. TRAV13-1)
						if( not re.match('^TR([ABGD])(?:([VDJ])(\d+(?:-\d+)?)|(?:(C)(\d*)))$', component) ):
								self.log.info("Invalid receptor description %s on line %d, ignoring", component, line_number)
								continue
						else:
								matches = re.match('^TR([ABGD])(?:([VDJ])(\d+(?:-\d+)?)|(?:(C)(\d*)))$', component)

								new_segment['gene'] = component
								new_segment['receptor_type'] = matches.groups()[0]
								new_segment['segment_type'] = matches.groups()[1] if matches.groups()[1] else matches.groups()[3]
								new_segment['segment_number'] = matches.groups()[2] if matches.groups()[2] else matches.groups()[4]

						if( not re.match('^(14q11.2|7q34|7p14)$', chromosome) ):
								self.log.info("Invalid chromosome %s on line %d, ignoring", chromosome, line_number)
								continue
						else:
								new_segment['chromosome'] = chromosome

						if( not re.match('^(forward|reverse)$', strand) ):
								self.log.info("Invalid strand %s on line %d, ignoring", strand, line_number)
								continue
						else:
								new_segment['strand'] = strand

						if( not re.match('^([VDJ]-REGION)|([VDJ]-GENE-UNIT)|(L-V-GENE-UNIT)|(EX\d+)|(L-PART1\+L-PART2)$', region) ):
								self.log.info("Invalid region: %s on line %d, ignoring", region, line_number)
								continue
						else:
								new_segment['region'] = region

						if( not re.match('^(\d+)\.\.(\d+)$', coordinates) ):
								self.log.info("Invalid coordinates %s on line %d, ignoring", coordinates, line_number)
								continue
						else:
								matches = re.match('^(\d+)\.\.(\d+)$', coordinates)
								new_segment['start_position'] = int(matches.groups()[0])
								new_segment['end_position'] = int(matches.groups()[1])

						self.log.info("Adding component: %s:%s [%s-%s]",
													new_segment['gene'], new_segment['region'],
													new_segment['start_position'], new_segment['end_position'])

						# Enforce uniqueness of our genes & regions (e.g. J-REGION of TRAJ24)
						for i in self.receptorSegment:
								if(i['gene'] == new_segment['gene'] and i['region'] == new_segment['region'] ):
										raise ValueError("Entry for gene and region pair was previously defined!", new_segment)
						self.receptorSegment.append(new_segment)
						#self.receptorSegment[component] = new_segment

				# Store a unique set of all gene names
				geneName = []
				for i in self.receptorSegment:
						geneName.append(i['gene'])
				self.geneName = set(geneName)

				return


		# readAlleles - Read allele info from fasta files with IMGT/GENE-DB compliant header lines
    #
    # See examples of these files at http://www.imgt.org/vquest/refseqh.html
    # This function populates the "allele" dict inside our receptorSegment
    # variable as described above
		#
		# Arguments:
		# filenames - A string (or array of strings) containing the file name of a fasta-formatted
    #             file with IMGT/GENE-DB headers describing gene allele sequences


		def readAlleles( self, filenames ):

				if isinstance(filenames, str):
						filenames = [filenames]

				for file in filenames:
						line_num = 0
						self.log.info("Processing file %s", file)
						content = re.sub('([ctag]{1})(\r?\n)+([ctag]{1}|\Z)(\n\Z)?', r'\1\3', open(file, 'rU').read()).split("\n")
						new_allele = None
						for line in content:
								line_num += 1
								if (
												( (line_num % 2 == 1) and re.search('^>.+$', line) ) or
												( (line_num % 2 == 0) and re.search('^[ctag]+$', line) )
								):
										if( line_num % 2 == 1 and len(line.split("|")) == 16 ):
												x,allele,x,functionality,region,x,x,x,x,x,x,x,x,x,x,x = line.split("|",)
												matches = re.match('^(TR[ABGD](?:[VDJ]\d+(?:-\d+)?|(?:C\d*)))\*(\d+)$', allele)
												if( matches is None ):
														self.log.warning("Skipping unsupported gene allele name %s", allele)		
														continue

												if re.match('^(V-REGION|J-REGION|EX\d|D-REGION|L-PART1\+L-PART2)$', region):
														gene, allele_name = matches.groups()
														
														segments = filter(lambda e:e['gene'] == gene and e['region'] == region,
																							self.receptorSegment)
														if( len(segments) is 1 ):
																self.log.debug("Found %s for %s*%s", region, gene, allele_name)

																# Populate new allele info, will complete addition after the sequence is read from file
																new_allele = {}
																new_allele['name'] = allele_name
																new_allele['index'] = self.receptorSegment.index(segments[0])
														elif( len(segments) is 0 ):
																self.log.warning("No corresponding receptor localization data for %s of %s*%s", region, gene, allele_name)
														else:
																self.log.critical("Multiple receptor localization data for %s of %s*%s", region, gene, allele_name)
												else:
														self.log.warning("Skipping unsupported gene region \"%s\" in allele %s", region, allele)
														continue
												
														
										elif( line_num % 2 == 1 ):
												self.log.error("In file %s, header on line %d does not appear to be in IMGT/GENE-DB format, ignoring...",
															 		file, line_num)
												continue
										elif( line_num %2 == 0 ):
												# Add any found new allele to our struct
												if( new_allele is not None):
														if( new_allele.get('index', None) is not None and
																new_allele.get('name', None) is not None ):
																
																if( self.receptorSegment[new_allele['index']].get('allele', None) is None ):
																		self.receptorSegment[new_allele['index']]['allele'] = {}
																#self.log.debug("Assigning %s sequence: %s", new_allele['name'], line)
																self.receptorSegment[new_allele['index']]['allele'][new_allele['name']] = line
																new_allele = None
																exit
										
								else:
										self.log.critical("In file %s, line %d does not appear to be FASTA formatted: \"%s\"",
																			file, line_num, line)
										raise ValueError("In file, line does not appear to be FASTA formatted: ",
																			file, line_num, line)
								
				self.log.info("readAlleles(): Processing complete")


		# setChromosomeFiles - Identify the location of necessary chromosome reference files
		# 
		# Arguments:
		# chr7  - File name of chromsome 7 reference file
		# chr14 - File name of chromosome 14 reference file
		# 
		# Returns:
		# nothing

		def setChromosomeFiles(self, chr7=None, chr14=None):

				if not chr7 is None:
						self.chr7Filename = chr7
						with open(self.chr7Filename) as fp:
								self.chr7Offset = len(fp.readline())
								self.chr7LineLength = len(fp.readline()) - 1
								self.log.debug("Chr14: Offset %d, line length %d", self.chr7Offset, self.chr7LineLength)
								if( self.chr7Offset > 10 ):
										self.log.warning("Chromosome file %s does not appear to be in proper format, reading may fail", chr7)
								
				if not chr14 is None:
						self.chr14Filename = chr14
						with open(self.chr14Filename) as fp:
								self.chr14Offset = len(fp.readline())
								self.chr14LineLength = len(fp.readline()) - 1
								self.log.debug("Chr14: Offset %d, line length %d", self.chr14Offset, self.chr14LineLength)
								if( self.chr14Offset > 10 ):
										self.log.warning("Chromosome file %s does not appear to be in proper format, reading may fail", chr14)
								
		# readChromosome - Request data from a chromosome reference
		# 
		# Arguments:
		# chromosome - Chromsome reference number.  Currently supported values are 7 and 14.
		# start - Start coordinates (IMGT)
		# end   - End coordinates (IMGT)
		# strand - Strand to read from.  Can be one of: forward, reverse.  Default is forward.
		
		def readChromosome(self, chromosome, start, end, strand):
				self.log.info("readChromosome() starting")
				self.log.debug("Arguments: %s %s %s %s", chromosome, start, end, strand)

				if start < 0 or end < 0:
						raise ValueError("Values must be non-zero integers")
				
				if not chromosome in (7, 14):
						raise ValueError("Chromosome must be either 7 or 14")
				elif chromosome == 7:
						filename = self.chr7Filename
						offset = self.chr7Offset
						lineLength = self.chr7LineLength
				elif chromosome == 14:
						filename = self.chr14Filename
						offset = self.chr14Offset
						lineLength = self.chr14LineLength
				else:
						raise ValueError("Unknown type passed for chromosome number")
				
				if start <= 0 or end <= 0:
						raise ValueError("Start and end values must be positive integers")
				
				with open(filename) as fp:
						self.log.info("Bytes requested %d, seek %d, offset %d, reading +%d",
													end-start,
													(end-start+int(math.floor(start/lineLength)) -1),
													offset,
													int(math.floor((end-start)/lineLength)))
						
						fp.seek(offset + start + int(math.floor(start/lineLength)) - 1)
						data = fp.read(end - start + int(math.floor((end-start)/lineLength)) + 2)
						data = data.replace("\n", '').upper()
						data = data[:(end - start + 1)]
						if strand == 'reverse':
								data = data.translate(string.maketrans('CATG', 'GTAC'))
								data = data[::-1]

				#self.log.debug("Read: %s (%db)", data, len(data))
				return data


				
		# chooseRandomSegment - Pick an (appropriately) random V, D, J or C segment for a provided receptor type
		#
		# Arguments:
		# receptorType -  A, B, G, or D. For the receptor type (A = alpha, B = beta, etc.)
    # componentName - V, D, J, or C. For the requested segment type (V - Variable, D - Diversity, etc.) 
		# V, D, J -       A 2-tuple consisting of an index to the V/D/J-REGION CDR3 component, and an allele name
		#
		# Returns:
		# A 2-tuple with an index and allele name to the requested component type
		def chooseRandomSegment(self, receptorType, componentName, V=None, D=None, J=None):
				self.log.debug("chooseRandomSegment() starting")
				self.log.debug("Arguments: %s, %s, %s, %s, %s", receptorType, componentName, V, D, J)

				if( receptorType not in ('A', 'B', 'G', 'D') ):
						raise ValueError("Receptor type must be either A, B, G, or D (alpha, beta, gamma or delta, respectively)")
				elif( receptorType in ('A', 'G') and componentName == 'D' ):
						return None
				elif( componentName not in ('V', 'D', 'J', 'C') ):
						raise ValueError("componentName must be one of V, D, J or C")
				elif( componentName is 'D' and V is None ):
						raise ValueError("Must define your V segment when choosing D segments")
				elif( componentName is 'J' and receptorType in ('A', 'G') and V is None ):
						raise ValueError("Must define your V segment when choosing alpha or gamma J segments")
				elif( componentName is 'J' and receptorType in ('B', 'D') and D is None ):
						raise ValueError("Must define your D segment when choosing beta or delta J segments")

				if V is not None:
						Vindex, Vallele = V
				if D is not None:
						Dindex, Dallele = D
				if J is not None:
						Jindex, Jallele = J

				# Generate a list of all valid segments (n.b. we pick CDR3 components here, not gene units)
				segmentChoices = []
				for i in range(0, len(self.receptorSegment)):
						if( re.match('^TR'+receptorType+componentName, self.receptorSegment[i]['gene'] ) and
								re.match('^[VDJ]-REGION|EX1', self.receptorSegment[i]['region']) ):

								if componentName == 'V':
										self.log.debug("Valid segment choice %s (%d)", self.receptorSegment[i]['gene'], i)
										segmentChoices.append(i)
										
								if( componentName == 'J' and
										self.receptorSegment[i]['chromosome'] == self.receptorSegment[Vindex]['chromosome'] ):
										if not D is None:
												if( self.receptorSegment[i]['start_position'] < self.receptorSegment[Dindex]['start_position'] and
														self.receptorSegment[i]['strand'] == 'forward' ):
														self.log.debug("This is not a valid choice: %s (%d)", self.receptorSegment[i]['gene'], i)
														continue
										self.log.debug("Valid segment choice %s (%d)", self.receptorSegment[i]['gene'], i)
										segmentChoices.append(i)

								if( componentName == 'D' ):
										self.log.debug("Valid segment choice %s (%d)", self.receptorSegment[i]['gene'], i)
										segmentChoices.append(i)

								if( componentName == 'C' and
										self.receptorSegment[i]['chromosome'] == self.receptorSegment[Vindex]['chromosome'] and
										( (self.receptorSegment[i]['start_position'] > self.receptorSegment[Jindex]['start_position'] and
											 self.receptorSegment[i]['strand'] == 'forward' ) or
											(self.receptorSegment[i]['start_position'] < self.receptorSegment[Jindex]['start_position'] and
											 self.receptorSegment[i]['strand'] == 'reverse' ) ) ):
										self.log.debug("Valid segment choice %s (%d)", self.receptorSegment[i]['gene'], i)
										if len(segmentChoices) > 0:
												a = self.receptorSegment[i]
												b = self.receptorSegment[segmentChoices[0]]
												if( ( a['strand'] == 'forward' and a['start_position'] < b['start_position']) or
														( a['strand'] == 'reverse' and a['start_position'] > b['start_position']) ):
														segmentChoices = [ i ]
												self.log.debug("C segment choices are now: %s", segmentChoices)
										else:
												segmentChoices.append(i)
												
				# Throw an error here if there are no possible join candidates.
				# This should never happen
				if len(segmentChoices) == 0:
						self.log.critical("No possible segments to join")
						exit(-10)

										
				segmentProbabilities = []						
		    # Pull out our predefined probabilities, moving them from segmentChoices to segmentProbabilities as we find them
				for i in self.VDJprobability:
						for j in segmentChoices:

								if ( componentName == 'V' and
										 len(i) == 2 and
										 self.receptorSegment[j]['gene'] in i ):
										segmentProbabilities.append((j, i[len(i)-1]))
										segmentChoices.remove(j)
										break

								if( componentName == 'D' and
										len(i) == 3 and
										self.receptorSegment[j]['gene'] in i and
										self.receptorSegment[Vindex]['gene'] in i ):
										segmentProbabilities.append((j, i[len(i)-1]))
										segmentChoices.remove(j)
										break

								if( componentName == 'J' and
										( len(i) == 3 and
											self.receptorSegment[j]['gene'] in i and
											self.receptorSegment[Vindex]['gene'] in i ) or
										( len(i) == 4 and
											self.receptorSegment[j]['gene'] in i and
											self.receptorSegment[Vindex]['gene'] in i and
											self.receptorSegment[Dindex]['gene'] in i ) ):
										segmentProbabilities.append((j, i[len(i)-1]))
										segmentChoices.remove(j)
										break


				# Fill in the undefined probabilities
				probabilityTotal = 0
				for i in segmentProbabilities:
						probabilityTotal += i[1]
				defaultProbability = float(1 - probabilityTotal) / len(segmentChoices)
				for i in segmentChoices:
						segmentProbabilities.append((i, defaultProbability))
						
				
				# Randomly choose a segment
				rand = random.random()
				self.log.debug("Roll: %0.3f", rand)

				cumulativeProbability = 0
				for i in segmentProbabilities:
						segmentIndex, probability = i
						cumulativeProbability += probability
						self.log.debug("Examining %s, cumulative: %0.3f", i, cumulativeProbability)
						if rand < cumulativeProbability:
								alleles = self.receptorSegment[segmentIndex]['allele'].keys()
								self.log.debug("Allele choices: %s", ', '.join(alleles))
								random.shuffle(alleles)
								allele = alleles[0]
								self.log.info("Choosing %d(%s) allele %s", segmentIndex, self.receptorSegment[segmentIndex]['gene'], allele)
								return (segmentIndex, allele)

				raise ValueError("We fell through the rabbit hole")



		# Perform the recombination of V, D, J and C segments.
    # This includes chewback and nucletide addition, as well as validating the new TCR
		# to ensure it is structurally correct (e.g. no errant stop codons, etc.)
    #
		# Return value:
    # 1. An array of two 6-tuples, defining our DNA and RNA sequences:
		#   ( chromosome, 5' end coordinates (NCBI), 5' strand (forward/reverse), DNA/RNA sequence, 3' start coordinates (NCBI), 3' strand (forward/reverse) )
    # OR
		# 2. None: if this VDJ recombination failed (i.e. early stop codon or invalid CDR3 AA sequence)
		
		def recombinate( self, V, D, J, C):
				self.log.info("Recombinate() called...")
				self.log.debug("Arguments: %s", (V, D, J, C))

				DNA = None
				RNA = None

				vIndex, vAllele = V
				jIndex, jAllele = J
				cIndex, cAllele = C
				
				chromosome = None
				if re.match('^7', self.receptorSegment[jIndex]['chromosome']):
						chromosome = 7
				elif re.match('^14', self.receptorSegment[jIndex]['chromosome']):
						chromosome = 14
						
				self.log.debug("Calculating V segment...")
				vSegmentDNA, vSegmentRNA = self.getSegmentSequences(V)
				vChewback = self.roll(self.junctionProbability['Vchewback'])
				if vChewback > 0:
						vSegmentDNA = vSegmentDNA[0:-vChewback]
						vSegmentRNA = vSegmentRNA[0:-vChewback]

				if D is not None:
						self.log.debug("Calculating D segment...")
						dSegmentDNA, dSegmentRNA = self.getSegmentSequences(D)
						d5Chewback = self.roll(self.junctionProbability['D5chewback'])
						d3Chewback = self.roll(self.junctionProbability['D3chewback'])
						vdAdditions = self.getRandomNucleotides(self.roll(self.junctionProbability['VDaddition']))
						djAdditions = self.getRandomNucleotides(self.roll(self.junctionProbability['DJaddition']))
						if d3Chewback > 0:
								dSegmentDNA = dSegmentDNA[d3Chewback:]
								dSegmentRNA = dSegmentRNA[d3Chewback:]
						if d5Chewback > 0:
								dSegmentDNA = dSegmentDNA[0:-d5Chewback]
								dSegmentRNA = dSegmentRNA[0:-d5Chewback]
						dSegmentDNA = vdAdditions + dSegmentDNA + djAdditions
						dSegmentRNA = vdAdditions + dSegmentRNA + djAdditions
				else:
						dSegmentDNA = self.getRandomNucleotides(self.roll(self.junctionProbability['VJaddition']))
						dSegmentRNA = dSegmentDNA
						
				self.log.debug("Calculating J segment...")
				jChewback = self.roll(self.junctionProbability['Jchewback'])
				jSegmentDNA, jSegmentRNA = self.getSegmentSequences(J)
				if jChewback > 0:
						jSegmentDNA = jSegmentDNA[jChewback:]
						jSegmentRNA = jSegmentRNA[jChewback:]

				self.log.debug("Calculating C segment...")
				cSegmentDNA, cSegmentRNA = self.getSegmentSequences(C)
				
				self.log.debug("Calculating JC segment DNA...")
				jcSegmentDNA = None
				if self.receptorSegment[jIndex]['strand'] == 'forward':
						jcSegmentDNA = self.readChromosome(chromosome, self.receptorSegment[jIndex]['end_position'] + 1, self.receptorSegment[cIndex]['start_position'] - 1, self.receptorSegment[cIndex]['strand'])
				else:
						jcSegmentDNA = self.readChromosome(chromosome, self.receptorSegment[cIndex]['end_position'] + 1, self.receptorSegment[jIndex]['start_position'] - 1, self.receptorSegment[cIndex]['strand'])

				# Assemble DNA/RNA sequences and validate it for early stops and functional CDR amino acid sequence
				dnaSequence = vSegmentDNA + dSegmentDNA + jSegmentDNA + jcSegmentDNA + cSegmentDNA
				rnaSequence = vSegmentRNA + dSegmentRNA + jSegmentRNA + cSegmentRNA

				self.log.debug("Validating RNA: %s", rnaSequence)

				# Ensure string is in-frame first...
				matches = re.match('^ATG((?:[CTAG]{3})+)$', rnaSequence)
				if matches is None:
						self.log.critical("Invalid CDR3: Frame shifted (%d)", len(rnaSequence)%3)
						return None
				else:
						self.log.critical("NOT Frame shifted")
				# Check for early stop codons
				matches = re.match('^((?:[CTAG]{3})*)(TAA|TAG|TGA)((?:[CTAG]{3})+)$', rnaSequence)
				if matches is not None:
						self.log.critical("Invalid CDR3: Stop codon found in sequence...")
						self.log.critical("%s", '.'.join(matches.groups()))
						self.log.critical("%s", '_'.join([vSegmentRNA, dSegmentRNA, jSegmentRNA, cSegmentRNA]))
						return None

				# Continue only if our CDR3 sequence is valid
				if not self.validateCDR3Sequence(rnaSequence):
						self.log.critical("Invalid CDR3: Amino acid sequence incorrect")
						return None
										
				# Calculate the location of our DNA/RNA 5' and 3' UTR areas
		    # We need to locate the corresponding L-V-GENE-UNIT to locate our start codon
				dnaStartPosition = None
				dnaEndPosition = None
				geneUnits = filter(lambda x:x['gene'] == self.receptorSegment[vIndex]['gene'] and
																	 re.match('^L-V-GENE-UNIT$', x['region']), self.receptorSegment)
				if len(geneUnits) == 0:
						self.log.critical("No corresponding GENE-UNIT found for gene %s", self.receptorSegment[segmentIndex]['gene'])
						exit(-10)
				elif len(geneUnits) > 1:
						self.log.critical("Too many GENE-UNITs found for gene %s", self.receptorSegment[segmentIndex]['gene'])
						exit(-10)

				geneUnit = geneUnits[0]
				geneCoordinates = (geneUnit['start_position'], geneUnit['end_position'], geneUnit['strand'])
				if self.receptorSegment[vIndex]['strand'] == 'forward':
						dnaStartPosition = geneUnit['start_position']
						rnaStartPosition = dnaStartPosition
						dnaEndPosition = self.receptorSegment[cIndex]['start_position'] + len(cSegmentDNA)
						rnaEndPosition = dnaEndPosition
				else:
						dnaStartPosition = geneUnit['end_position']
						rnaStartPosition = dnaStartPosition
						dnaEndPosition = self.receptorSegment[cIndex]['end_position'] - len(cSegmentDNA)
						rnaEndPosition = dnaEndPosition
						
				startStrand = self.receptorSegment[vIndex]['strand']
				endStrand = self.receptorSegment[cIndex]['strand']
				
				DNA = (chromosome, dnaStartPosition, startStrand, dnaSequence, dnaEndPosition, endStrand)
				RNA = (chromosome, rnaStartPosition, startStrand, rnaSequence, rnaEndPosition, endStrand)
				self.log.debug("recombinate() returning DNA: %s \nRNA %s", DNA, RNA)
				self.log.info("recombinate() returning...")
				return [ DNA, RNA ]

		

		# validateCDR3Sequence - Determine if an RNA sequence represents a valid CDR3
    # 
    # Arguments:
    # rnaSeq - RNA sequence to examine
    # 
    # Returns:
    # Boolean - True if there is a valid CDR3 sequence, in frame, in rnaSeq
    #           False, if otherwise
    # 
		
		def validateCDR3Sequence(self, rnaSeq):
				
        # Enforce CxxxxxFGxG in AA space
				cdr3Pattern = '^((?:[CTAG]{3})+)(TG[TC])((?:[CTAG]{3}){5,32})(TT[TC]GG[CTAG][CTAG]{3}GG[CTAG])'
				## Check for the AA sequence following CxxxxxF
				#cdr3Consolation = '^((?:[CTAG]{3})+)(TG[TC])((?:[CTAG]{3}){5,32})(TT[TC](?:[CTAG]{3}){3})'

				matches = re.match(cdr3Pattern, rnaSeq)
				if matches is not None:
						self.log.info("Valid CDR3")
						self.log.debug(matches.groups())
				else:
						#matches = re.match(cdr3Consolation, rnaSequence)
						self.log.info("Invalid CDR3")
						#if matches is not None:
						#		self.log.debug(matches.groups())
						return False

				return True

		# getCDR3Sequence - Obtain the CDR3 RNA sequence from an RNA sequence
    # 
    # Arguments:
    # rnaSeq - RNA sequence to examine
    # 
    # Returns:
    # String, or None - String containing the nucletides of the CDR3 sequence,
		#                   if a valid CDR3 sequence was found.
		#                   None, if otherwise
    # 
		def getCDR3Sequence(self, rnaSeq):

        # Enforce CxxxxxFGxG in AA space
				cdr3Pattern = '^((?:[CTAG]{3})+)(TG[TC])((?:[CTAG]{3}){5,32})(TT[TC]GG[CTAG][CTAG]{3}GG[CTAG])'

				matches = re.match(cdr3Pattern, rnaSeq)
				if matches is not None:
						return(''.join(matches.groups()[1:4]))
				return None
				

		# getRandomNucleotides - Generate random nucleotide strings
    # 
    # Arguments:
    # count - Number of random bases
    # 
    # Returns:
    # String - String of count random C, T, A, or G characters
    # 
    # 
    # 

		def getRandomNucleotides(self, count):
				if count < 0:
						raise ValueError("Count must be a non-negative integer (zero is permissible)")
				return ''.join(random.choice('CATG') for i in range(count))


		
		# getSegmentSequences - Return the DNA & RNA of a gene segment
		#
		# N.b. The segment argument will likely be a 2-tuple pointing to a CDR3
		#      component and an allele name. This function will search out the
		#      GENE-UNIT component and return the DNA sequence for that with the
		#      appropriate allele DNA spliced in.
		#
		# Arguments:
		# segment - 2-tuple of 1) index into self.receptorSegment and 2) an allele name
		#
		# Returns:
		# Array, with strings representing the DNA and RNA (respectively) of this segment
		#
		
		def getSegmentSequences( self, segment ):
				self.log.info("getSegmentSequences() called...")
				self.log.debug("Arguments: %s", segment)
				if not len(segment) == 2:
						raise ValueError("Argument must be a 2-tuple")
				
				segmentIndex, segmentAllele = segment
				chromosome = None
				if(re.match('^7', self.receptorSegment[segmentIndex]['chromosome']) ):
						chromosome = 7
				elif(re.match('^14', self.receptorSegment[segmentIndex]['chromosome']) ):
						chromosome = 14
				else:
						raise ValueError("Unknown chromosome " + self.receptorSegment[segmentIndex]['chromosome'])

				self.log.debug("Segment sequence requested: %s", self.receptorSegment[segmentIndex])
				
				# If a [VDJ]-REGION provided, find the GENE-UNIT for the given segment
				if re.match('^[VDJ]-REGION', self.receptorSegment[segmentIndex]['region'] ):
						# Replace the V-REGION and L-PART1+L-PART2 sequences within the GENE-UNIT sequence
						if self.receptorSegment[segmentIndex]['region'] == 'V-REGION':
								geneUnits = filter(lambda x:x['gene'] == self.receptorSegment[segmentIndex]['gene'] and
																	 re.match('^(L-V|D|J)-GENE-UNIT$', x['region']), self.receptorSegment)
								if len(geneUnits) == 0:
										self.log.critical("No corresponding GENE-UNIT found for gene %s", self.receptorSegment[segmentIndex]['gene'])
										exit(-10)
								elif len(geneUnits) > 1:
										self.log.critical("Too many GENE-UNITs found for gene %s", self.receptorSegment[segmentIndex]['gene'])
										exit(-10)

								geneUnit = geneUnits[0]

								geneCoordinates = (geneUnit['start_position'], geneUnit['end_position'], geneUnit['strand'])
								alleleCoordinates = (self.receptorSegment[segmentIndex]['start_position'], self.receptorSegment[segmentIndex]['end_position'], self.receptorSegment[segmentIndex]['strand'])


								geneData = self.readChromosome(chromosome, geneCoordinates[0], geneCoordinates[1], geneCoordinates[2])

								geneHeaderLength = alleleCoordinates[0] - geneCoordinates[0]
								geneAlleleLength = alleleCoordinates[1] - alleleCoordinates[0] + 1
								if self.receptorSegment[segmentIndex]['strand'] == 'reverse':
										geneHeaderLength = abs(geneCoordinates[1] - alleleCoordinates[1])
										geneAlleleLength = abs(alleleCoordinates[1] - alleleCoordinates[0]) + 1
										
								self.log.debug("Header %d, Allele %d, total %d",
															 geneHeaderLength,
															 geneAlleleLength,
															 len(geneData))
						
								self.log.debug("Head:   %s", geneData[0:geneHeaderLength])
								self.log.debug("Allele: %s", self.receptorSegment[segmentIndex]['allele'][segmentAllele])
								self.log.debug("Tail:   %s", geneData[geneHeaderLength+geneAlleleLength:])
								dnaData = geneData[0:geneHeaderLength] + self.receptorSegment[segmentIndex]['allele'][segmentAllele]
								# Now substitute the L-PART allele in RNA (as the geneHeader portion has an intron in it)
								lPartSegments = filter(lambda x:x['gene'] == self.receptorSegment[segmentIndex]['gene'] and
																			 x['region'] == 'L-PART1+L-PART2', self.receptorSegment)
								if len(lPartSegments) == 1:
										if segmentAllele in lPartSegments[0]['allele']:
												rnaData = lPartSegments[0]['allele'][segmentAllele] + self.receptorSegment[segmentIndex]['allele'][segmentAllele]
										else:
												rnaData = lPartSegments[0]['allele'][random.choice(lPartSegments[0]['allele'].keys())] + self.receptorSegment[segmentIndex]['allele'][segmentAllele]
										
								else:
										self.log.error("Did not find matching L-PART segment for this V-REGION")
										exit(-10)

								self.log.debug("Returning data for V segment")
								return [ dnaData.upper(), rnaData.upper() ]
						
						elif self.receptorSegment[segmentIndex]['region'] == 'D-REGION':
								dnaData = self.receptorSegment[segmentIndex]['allele'][segmentAllele]
								rnaData = dnaData
								self.log.debug("Returning data for D segment")
								return [ dnaData.upper(), rnaData.upper() ]

						elif self.receptorSegment[segmentIndex]['region'] == 'J-REGION':
								dnaData = self.receptorSegment[segmentIndex]['allele'][segmentAllele]
								rnaData = dnaData
								self.log.debug("Returning data for J segment")
								return [ dnaData.upper(), rnaData.upper() ]
								
								
				# If an EX1 provided, pull ALL of the exons for that C-segment
		    # RNA is just sum of the alleles, but DNA includes intronic segments as well
				elif re.match('^EX1', self.receptorSegment[segmentIndex]['region']):
						cStartPosition = None
						cEndPosition = None
						rnaData = ['', '', '', '']
						for i in filter(lambda j:re.match('^EX\d$', j['region']) and j['gene'] == self.receptorSegment[segmentIndex]['gene'], self.receptorSegment):
								if i['start_position'] < cStartPosition or cStartPosition is None:
										cStartPosition = i['start_position']
								if i['end_position'] > cEndPosition or cEndPosition is None:
										cEndPosition = i['end_position']
								if i['region'] == 'EX1':
										rnaData[0] = i['allele'][random.choice(i['allele'].keys())]
								if i['region'] == 'EX2':
										rnaData[1] = i['allele'][random.choice(i['allele'].keys())]
								if i['region'] == 'EX3':
										rnaData[2] = i['allele'][random.choice(i['allele'].keys())]
								if i['region'] == 'EX4' and 'allele' in i:
										rnaData[3] = i['allele'][random.choice(i['allele'].keys())]
						rnaData = ''.join(rnaData)
						#self.log.debug("Returning RNA: %s", rnaData)
						dnaData = self.readChromosome(chromosome, cStartPosition, cEndPosition, self.receptorSegment[segmentIndex]['strand'])
						#self.log.debug("Returning DNA: %s", dnaData)
						self.log.debug("Returning data for C segment")
						return [ dnaData.upper(), rnaData.upper() ]

				self.log.critical("We shouldn't be here")
				exit(-10)
				return

				
		# Choose an index from a probability array
		# Arguments:
		# probability - array of probabilities (e.g. [0.5, 0.25, 0.125, 0.125] )
		#               These values SHOULD total to 100%, as any remaining probability is implicitly assigned to the last index
		# Returns:
		# Integer representing the index to the probability array randomly chosen based on the provided probabilities
		# (e.g. 0, 1, 2, or 3 in the above example)
		#
		def roll( self, probability ):
				rand = random.random()
				cumulativeProbability = 0
				for i in range(0, len(probability)):
						cumulativeProbability += probability[i]
						index = i
						if rand < cumulativeProbability:
								break
				return index


		# Degrade a sequence read, based on some parameters
    #
		# Arguments:
    # read   - String.  The read to be degraded
		# method - String.  The method for degradation to be performed. Valid values: "logistic" or "phred"
    # ident  - String. Label for this entry in a FASTQ formatted file
		# variability   - Float.  Optional.  This adds additional per-nucleotide error
    #          to the methods described below.  This is a percentage
		#          relative error e.g. For some base, if the method determines the
    #          error to be 0.05 or 5% and fuzz was given as 0.5 or 50%, then 
		#          the actual error is 0.05 +/- 0.025, or range [2.5%, 7.5%]
    #          Default is zero.
		# display - Boolean.  If True, then will print details of degradation to
    #           stdout.  Default is False.
    # [ phred ] are parameters for using a Phred score to determine the per-nucleotide error rate
		# phred - Error rate probability per nucleotide.  This is Phred+33 format
    #         (i.e. range [!, I] as acceptable characters)
    #         If the read being degraded is longer than the given Phred string,
		#         then the last base of the phred string is used for subsequent bases
    #         (e.g. "IIIIJJ55" is equivalent to "IIIIJJ5555555555555")
    # [ baseError, L, k, midpoint] are parameters for a logistic function defining the error rate:
    # baseError - Base error rate probability
		# L - Maximum error rate
    # k - Steepness factor
		# midpoint - Base position where the error rate is equal to 1/2 of L
		#
    # Returns:
    # FASTQ string of the "degraded" read, with sequence label and quality score
		#
		def getDegradedFastq(self, read, method, ident, variability=0, phred='', baseError=0, L=0, k=0, midpoint=0, display=False):
				self.log.info("getDegradedFastq() called")
				self.log.debug("Arguments: %s", (read, method, ident, variability, phred, baseError, L, k, midpoint))

				if display == True:
						print("Displaying degradation output with method %s, variability %0.5f" % (method, variability))
				
				phred33Reference = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI"
        #                   0 1234567890123456789012345678901234567890
				if method not in ('phred', 'logistic'):
						raise ValueError("Method %s must be either logistic or phred" % method)
						exit(-10)

				if method == 'logistic':
						readStr = ''
						qualStr = ''
						for i in range(0, len(read)):
								errorRate = (L - baseError) / (1 + math.exp(-k*(i - midpoint))) + baseError

								if variability != 0:
										errorRate += random.random() * 2 * errorRate * variability - errorRate * variability

								phredScore = int(-10 * math.log(errorRate))
								if phredScore > 40:
										phredScore = 40
										if random.random() < errorRate:
												readStr += self.getRandomNucleotides(1)
										else:
												readStr += read[i]
										qualStr += phred33Reference[phredScore]
										if display == True:
												print("Position %02d: error rate %0.4f, Phred+33 %s" % (i, errorRate, phred33Reference[phredScore]))

						fastqOutput = "%s\n%s\n+\n%s\n" % (ident, readStr, qualStr)
				
						self.log.debug("Orig: %s", read)
						self.log.debug("Read: %s", readStr)
						self.log.debug("Qual: %s", qualStr)

						return fastqOutput

				elif method == 'phred':
						readStr = ''
						qualStr = ''
						for i in range(0, len(read)):
								errorRate = 0
								if i >= len(phred):
										errorRate = 10 ** ( phred33Reference.index(phred[-1]) / float(-10) )
								else:
										errorRate = 10 ** ( phred33Reference.index(phred[i]) / float(-10) )
								if errorRate > 1:
										errorRate = 1
								elif errorRate < 0:
										errorRate = 0
								
								if variability != 0:
										errorRate += random.random() * 2 * errorRate * variability - errorRate * variability
								
								if random.random() < errorRate:
										readStr += self.getRandomNucleotides(1)
								else:
										readStr += read[i]
										
								phredScore = int(-10 * math.log10(errorRate))
								if phredScore > 40:
										phredScore = 40
										
								qualStr += phred33Reference[phredScore]

								if display is True:
										print("Position %02d: error rate %0.4f, Phred+33 %s" % (i, errorRate, phred33Reference[phredScore]))

								
						fastqOutput = "%s\n%s\n+\n%s\n" % (ident, readStr, qualStr)
						return fastqOutput
				else:
						raise ValueError("Invalid method \"%s\". We should not be here" % method)





				
class tcr:

		def __init__( self, AB_frequency, config, log=None ):
				if( isinstance(config, tcrConfig) ):
						self.config = config
				else:
						raise ValueError("config object must be a tcrConfig")

				self.setLog(log)

				if( AB_frequency >= 0 and
						AB_frequency <= 1):
						self.AB_frequency = AB_frequency
				else:
						raise ValueError("AB_frequency must be in range [0, 1]")
						
				self.type1 = None
				self.V1 = None
				self.D1 = None
				self.J1 = None
				self.C1 = None
				self.type2 = None
				self.V2 = None
				self.D2 = None
				self.J2 = None
				self.C2 = None

				# These tuples define our DNA and RNA sequences
				# Consist of a 6-tuple with the following definition:
				# ( chromosome,  5' end coordinates (NCBI), 5' strand (forward/reverse), DNA/RNA sequence, 3' start coordinates (NCBI), 3' strand (forward/reverse) )
				self.DNA1 = ()
				self.RNA1 = ()
				self.DNA2 = ()
				self.RNA2 = ()

				return

		# setLog - Configure our logging object
		# 
		# args:
		# log - If a logging object, we will use this for our logging
		#       If None, we will configure a new, non-functioning logging object
		#       
		# Returns:
		#  nothing
		# 
		def setLog( self, log ):
				if( isinstance(log, logging.Logger) ):
						self.log = log
				elif log is None:
						self.log = logging.getLogger(__name__)
						self.log.setLevel(99) # A high level, effectively disabling logging
				else:
						raise ValueError("Log object for tcr must be a logging.Logger (or None)")


				
		# rmLog - Remove our logging object
		#         This is usually called prior to serializing this object, as 
		#         logger objects cannot be serialized. see setLog()
		# 
		# Arguments: none
		# Returns: nothing
		# 
		def rmLog( self ):
				self.log = None

				
		# freeze - Render this self object suitable for pickling (with pickle or
    #          cPickle) Mostly this just involves discarding our self.log 
		#          logger objects, as these contain filehandles that cannot
		#          be serialized.
		# 
		# Arguments:
		# none
		# 
		# Returns:
		# A copy of this object
    #
		def freeze( self ):
				self.rmLog()
				self.config = None
				return self

		# thaw - Recover this object after being serialized
		#        Currently, this involves re-establishing the logger and
		#        config objects
    #
		# Arguments:
		# log - Optional.  A logger object to log to.  If empty or None, will 
		#       disable logging
		# config - Mandatory.  A tcrConfig object to use
    #
		# Returns: Nothing
		#
		def thaw( self, log=None, config=None ):
				self.setLog(log)
				if config is None:
						self.log.critical("Configuration option is mandatory")
						exit(-10)
				self.config = config
				
		
		def randomize( self ):
				self.log.info("Starting randomize()")
				if( random.random() <= self.AB_frequency ):
						self.type1 = 'A'
						self.type2 = 'B'
				else:
						self.type1 = 'G'
						self.type2 = 'D'
				self.log.info("Chosen: %s %s", self.type1, self.type2)

				while 1:
						self.V1 = self.config.chooseRandomSegment(self.type1, componentName='V')
						if self.type1 in ('B', 'D'):
								self.D1 = self.config.chooseRandomSegment(self.type1, componentName='D', V=self.V1)
						self.J1 = self.config.chooseRandomSegment(self.type1, componentName='J', V=self.V1, D=self.D1)
						self.C1 = self.config.chooseRandomSegment(self.type1, componentName='C', V=self.V1, D=self.D1, J=self.J1)
								
						sequenceTuple = self.config.recombinate(self.V1, self.D1, self.J1, self.C1)
						if sequenceTuple is not None:
								self.DNA1, self.RNA1 = sequenceTuple
								break
				while 1:
						self.V2 = self.config.chooseRandomSegment(self.type2, componentName='V')
						if self.type2 in ('B', 'D'):
								self.D2 = self.config.chooseRandomSegment(self.type2, componentName='D', V=self.V2)
						self.J2 = self.config.chooseRandomSegment(self.type2, componentName='J', V=self.V2, D=self.D2)
						self.C2 = self.config.chooseRandomSegment(self.type2, componentName='C', V=self.V2, D=self.D2, J=self.J2)
										
						sequenceTuple = self.config.recombinate(self.V2, self.D2, self.J2, self.C2)
						if sequenceTuple is not None:
								self.DNA2, self.RNA2 = sequenceTuple
								self.log.info("randomize() complete Frame shifted")
								return

		# getCDR3Sequences - Return RNA sequences of the CDR3 regions
    # 
    # Arguments: None
    # 
    # Returns:
    # Array to CDR3 nucleotide strings for the two components of this T cell
    # 
    # 
		def getCDR3Sequences( self ):
				return [ self.config.getCDR3Sequence(self.RNA1[3]), self.config.getCDR3Sequence(self.RNA2[3]) ]

						
class tcrRepertoire:

		def __init__( self, config, size, log=None, AB_frequency = 0.9, uniqueCDR3 = False, uniqueChain = False, uniqueTCR = False ):
				if( isinstance(config, tcrConfig) ):
						self.config = config
				else:
						raise ValueError("config object must be a tcrConfig")

				self.setLog(log)

				self.log.info("tcrRepertoire()::__init__ called with size %d, AB ratio %f, unique chains %s, unique CDR3 %s", size, AB_frequency, uniqueChain, uniqueTCR)
			 	self.AB_frequency = AB_frequency
				self.repertoire = [None] * size
				for i in range(0, size):
						self.log.debug("Generating repertoire bucket %d of %d", i + 1, size)
						self.repertoire[i] = tcr(self.AB_frequency, self.config, log=self.log.getChild('tcr'))

						if uniqueCDR3 == True: # Ensure unique CDR3
								unique = False
								while unique == False:
										self.repertoire[i].randomize()
										for j in range(0, i):
												if ( self.config.getCDR3Sequence(self.repertoire[i].RNA1[3]) == self.config.getCDR3Sequence(self.repertoire[j].RNA1[3]) or
														 self.config.getCDR3Sequence(self.repertoire[i].RNA2[3]) == self.config.getCDR3Sequence(self.repertoire[j].RNA2[3]) ) :
														self.log.debug("Duplicate CDR3 at position %d", j)
														break
										else:
												unique = True
												
						elif uniqueChain == True: # Ensure unique chains
								unique = False
								while unique == False:
										self.repertoire[i].randomize()
										for j in range(0, i):
												if ( self.repertoire[i].RNA1 == self.repertoire[j].RNA1 or
														 self.repertoire[i].RNA2 == self.repertoire[j].RNA2 ) :														
														self.log.debug("Duplicate chain at position %d", j)
														break
										else:
												unique = True

						elif uniqueTCR == True: # Ensure unique TCR
								unique = False
								while unique == False:
										self.repertoire[i].randomize()
										for j in range(0, i):
												if ( self.repertoire[i].RNA1 == self.repertoire[j].RNA1 and
														 self.repertoire[i].RNA2 == self.repertoire[j].RNA2 ) :
														self.log.debug("Duplicate CDR3 at position %d", j)
														break
										else:
												unique = True
						else: # No uniqueness constraints
								self.repertoire[i].randomize()
								
						self.log.debug("Finished generating repertoire bucket %d of %d", i + 1, size)
				self.population = [0] * size
				self.population_size = 0
				self.distribution_options = ('flat', 'equal', 'gaussian', 'chisquare')
				self.distribution = None

				return
		
		# setLog - Configure our logging object
		# 
		# Arguments:
		# log - If a logging object, we will use this for our logging
		#       If None, we will configure a new, non-functioning logging object
		#       
		# Returns:
		#  nothing
		# 
		def setLog( self, log ):
				if( isinstance(log, logging.Logger) ):
						self.log = log
				elif log is None:
						self.log = logging.getLogger(__name__)
						self.log.setLevel(99) # A high level, effectively disabling logging
				else:
						raise ValueError("Log object for tcrConfig must be a logging.Logger (or None)")


				
		# rmLog - Remove our logging object
		#         This is usually called prior to serializing this object, as 
		#         logger objects cannot be serialized. see setLog()
		# 
		# Arguments: none
		# Returns: nothing
		# 

		def rmLog( self ):
				self.log = None


		# freeze - Render this self object suitable for pickling (with pickle or
    #          cPickle) Mostly this just involves discarding our self.log 
		#          logger objects, as these contain filehandles that cannot
		#          be serialized.
		# 
		# Arguments:
		# none
		# 
		# Returns:
		# A copy of this object
    #
		def freeze( self ):
				self.rmLog()
				for i in self.repertoire:
						i.freeze()
				self.config = None
				return self

		# thaw - Recover this object after being serialized
		#        Currently, this involves re-establishing the logger objects
		# 
		# Arguments:
		# log - Optional.  A logger object to log to.  If empty or None, will 
		#       disable logging
		# config - Mandatory.  A tcrConfig object to use.
    #
		# Returns: Nothing
		#
		def thaw( self, log=None, config=None ):
				self.setLog(log)
				if config is None:
						self.log.critical("tcrRepertoire.thaw(): Configuration option is mandatory")
						exit(-10)
				for i in self.repertoire:
						i.thaw(self.log.getChild('tcr'), config=config )
				self.config = config
				

		# populate - Populate the repertoire with T cells
		#
		# Arguments:
		# population_size - Integer.  The number of T cells in this population.
		#                   When population_size > repertoire size (as is normal)
		#                   then some TCR CDR3 will be shared between individuals.
		# distribution -    String.  Describes the distribution of individual cells
		#                   across the repertiore.  Current options are:
		#                   flat: Distribute cells across all repertoires equally
		#                   even: Equal odds across all repertoires (may approach flat)
		#                   gaussian: Gaussian distribution with g_cutoff standard
		#                   deviations included in the repertiore
		#                   chisquare: Chi Square distribution.  Takes cs_k and
		#                   cs_cutoff arguments to describe k value and the largest
		#                   x-axis value to include in the distribution, respectively
		#
		# Returns: Nothing
		#		
		def populate( self, population_size, distribution, g_cutoff=3, cs_k=2, cs_cutoff=8 ):
				self.log.info("populate() called...")
				self.log.debug("Arguments: %s", (population_size, distribution, g_cutoff, cs_k, cs_cutoff))
				
				if(distribution not in self.distribution_options):
						raise ValueError("distribution must be one of ", self.distribution_options)
				else:
						self.distribution = distribution

				if(population_size > 0):
						self.population_size = population_size
				else:
						raise ValueError("population size must be a positive integer")
						
				if self.distribution is 'flat':
						for i in range(0, self.population_size):
								random_bin = int(math.floor(random.random() * len(self.repertoire)))
								self.population[random_bin] += 1

				elif self.distribution is 'equal':
						for i in range(0, self.population_size):
								random_bin = int(i % len(self.repertoire))
								self.population[random_bin] += 1
						
				elif self.distribution is 'gaussian':
						if( g_cutoff < 0 ):
								raise ValueError("Argument sd_cutoff for populate must be a positive integer")

						population_generated = 0						
						while( population_generated < self.population_size ):
								self.log.info("Rolling %d individuals", self.population_size - population_generated)
								
								for f in numpy.random.normal(0, 1, self.population_size - population_generated):
										if( abs(f) > g_cutoff ):
												continue
										# Scale the normal values to match our repertoire "buckets"
										f = f + g_cutoff
										bucket_size = g_cutoff * 2.0 / len(self.repertoire)
										self.population[int(f / bucket_size)] += 1
										population_generated += 1

				elif self.distribution is 'chisquare':
						if( cs_k <= 0 or cs_cutoff <= 0 ):
								raise ValueError("Invalid arguments for chi-square distribution.  Must be k > 0 and cutoff >0.  Got %f, %f" % (cs_k, cs_cutoff))
						population_generated = 0
						while population_generated < self.population_size:
								self.log.info("Rolling %d individuals", self.population_size - population_generated)
								for x in numpy.random.chisquare(cs_k, self.population_size - population_generated):
										if( x >= cs_cutoff):
												continue
										# Scale the values to match our repertoire "buckets"
										self.population[int((x/cs_cutoff) * len(self.repertoire))] += 1
										population_generated += 1

				else:
						raise ValueError("Invalid distribution %s" % self.distribution)
				
				self.log.info("populate() complete...")



				
		# simulateRead - Generate reads from a repertoire
		#
		# Arguments:
		# count        - Integer. The total # of reads requested
		# space        - String. The type of sequencing read to generate, 'dna' or
		#                'rna'
		#
		# Optional arguments:
		# paired_end   - Boolean.  True for paired end reads, False for single.
		#                Default is False.
		# distribution - String.  The distribution of the read and inner mate
		#                lengths.  Only option supported currently is 'gaussian'
		# read_length_mean      - Integer.  The mean read length to generate.
		# read_length_sd        - Integer. The standard deviation of the
		#                         distribution of lengths of generated reads.
		#                         Setting this to zero will produce fix-length
		#                         reads of read_length_mean length.  Default
		#                         is 4.
		# read_length_sd_cutoff - Integer.  The number of standard deviations
		#                         from the mean read length to include in generated
		#                         reads.
		# inner_mate_length_mean      - Integer.  The mean length of the inner mate
		#                               portion of paired end reads.  Used only
		#                               when generating paired end reads.  Setting
		#                               this to zero will produced fix-length
		#                               inner mate portions at the mean. Default
		#                               is 100 nucleotides.
		# inner_mate_length_sd        - Integer.  The standard deviation of the
		#                               distribution of the length of the inner
		#                               mate portion of generated reads.  Default
		#                               is 8.
		#
		# inner_mate_length_sd_cutoff - Integer.  The number of standard
		#                               deviations from the mean inner mate length
		#                               to include in generated reads.
		# 
		# Returns:
		#
		# A single 2-tuple (reads, comments), where:
		#
		# reads - If single-end reads, this is an array of strings.  If paired-end
		#         reads, this is an array of 2-tuples of (read1, read2)
		#         with types (string, string).
		#
		# comments - An array of strings.  These are descriptive strings that
		#            provide information regarding the reads in the 'reads'
		#            array, where comments[n] describes reads[n].
		#
		#
		def simulateRead( self, count, space, distribution='gaussian', read_length_mean=25, read_length_sd=4, read_length_sd_cutoff=4, paired_end=False, inner_mate_length_mean=100, inner_mate_length_sd=8, inner_mate_length_sd_cutoff=4):
				self.log.debug("simulateRead() called...")

				if space not in [ 'dna', 'rna' ]:
						self.log.critical("simulateRead() argument 2 must be either 'dna' or 'rna'")
						exit(-10)

				outputReads = []
						
				# Choose an individual to read from (a TCR chain [e.g. alpha or beta] is chosen later)
				readIndividual = None
				for i in range(0, count):
						randIndividual = random.random() * self.population_size
						self.log.debug("Read from individual #%d (of %s)", randIndividual, self.population_size)
						cumulativePopulation = 0
						for j in range(0, len(self.repertoire)):
								cumulativePopulation += self.population[j]
								if randIndividual < cumulativePopulation:
										self.log.debug("Individual is instance of cell %d in repertoire", j)
										readIndividual = j
										break
						outputComment='@STIG:readnum=%d:clone=%d' % (i, j)
								
						# Calculate our read length for this particular read
						readLength = None
						if distribution == 'gaussian':
								if paired_end is False:
										if read_length_sd > 0:
												randLength = 0
												while abs(randLength - read_length_mean) / read_length_sd > read_length_sd_cutoff and randLength <= 0:
														randLength = int(round(numpy.random.normal(read_length_mean, read_length_sd)))
												readLength = (randLength)
										else:
												readLength = read_length_mean
										
								else: # Paired-end reads
										read1Length = 0
										read2Length = 0
										innerMateLength = 0
										if read_length_sd > 0:
												while abs(read1Length - read_length_mean) / read_length_sd > read_length_sd_cutoff and read1Length <= 0:
														read1Length = int(round(numpy.random.normal(read_length_mean, read_length_sd)))
												while abs(read2Length - read_length_mean) / read_length_sd > read_length_sd_cutoff and read2Length <= 0:
														read2Length = int(round(numpy.random.normal(read_length_mean, read_length_sd)))
												while abs(innerMateLength - inner_mate_length_mean) / inner_mate_length_sd > inner_mate_length_sd_cutoff and innerMateLength <= 0:
														innerMateLength = int(round(numpy.random.normal(inner_mate_length_mean, inner_mate_length_sd)))
										else:
												read1Length = read_length_mean
												read2Length = read_length_mean
												innerMateLength = inner_mate_length_mean
										readLength = (read1Length, innerMateLength, read2Length)
														
						else:
								self.log.critical("TODO:Implement non-gaussian distributions?")
								exit(-10)
								#TODO
								# Implement non-gaussian distributions?
										
										
						self.log.debug("Read length for this read will be: %s", readLength)

						# Pick a location within this individual's DNA and generate the read
						totalReadLength = readLength if isinstance(readLength, int) else readLength[0] + readLength[1] + readLength[2]

						# Pick a chain to read from (alpha / beta or gamma / delta)
						receptorCoordinates = None
						if random.random() < 0.5:
								if space == 'dna':
										receptorCoordinates = self.repertoire[readIndividual].DNA1
								elif space == 'rna':
										receptorCoordinates = self.repertoire[readIndividual].RNA1
								outputComment = (outputComment + ":chain=%s" % self.repertoire[readIndividual].type1)
						else:
								if space == 'dna':
										receptorCoordinates = self.repertoire[readIndividual].DNA2
								elif space == 'rna':
										receptorCoordinates = self.repertoire[readIndividual].RNA2
								outputComment = (outputComment + ":chain=%s" % self.repertoire[readIndividual].type2)

						chromosome, sequenceStart, strandStart, sequence, sequenceEnd, strandEnd = receptorCoordinates

						# Determine a location within that chain's sequence and pull the read
						self.log.debug("Choosing between [-100, %d]", len(sequence) + 100 - totalReadLength)
						
						outputSequence = ''
						randomStart = random.choice(range(-100, len(sequence) + 100 - totalReadLength ))
						outputComment = (outputComment + ":relpos=%d" % randomStart)
						
						_5UTRBases = 0
						if   randomStart < 0 and abs(randomStart) > totalReadLength:
								_5UTRBases = totalReadLength
						elif randomStart < 0:
								_5UTRBases = abs(randomStart)

						_3UTRBases = 0
						if   randomStart >= 0 and (len(sequence) - randomStart) < totalReadLength:
								_3UTRBases = totalReadLength - (len(sequence) - randomStart)
						elif randomStart <  0 and (len(sequence) + abs(randomStart)) < totalReadLength:
								_3UTRBases = totalReadLength - (len(sequence) + abs(randomStart))

						self.log.debug("Starting read at position %d, 5p %db 3p %db", randomStart, _5UTRBases, _3UTRBases)

								
						if _5UTRBases > 0:
								outputSequence = self.config.readChromosome(chromosome, sequenceStart - _5UTRBases + 1, sequenceStart, strandStart)

						if _5UTRBases > 0 and _5UTRBases < totalReadLength:
								outputSequence += sequence[0:totalReadLength - _5UTRBases]
						elif _5UTRBases > 0 and _5UTRBases >= totalReadLength:
								outputSequence # Do nothing
						elif _5UTRBases == 0 and _3UTRBases == 0:
								outputSequence += sequence[randomStart:randomStart+totalReadLength]
						elif _3UTRBases > 0:
								outputSequence += sequence[len(sequence) - totalReadLength + _3UTRBases:]
						else:
								self.log.debug("We shouldn't be here")
								raise ValueError("We shouldn't be here")
								
						if _3UTRBases > 0:
								outputSequence += self.config.readChromosome(chromosome, sequenceEnd, sequenceEnd + _3UTRBases - 1, strandEnd)

						if paired_end is False:
								outputReads.append((outputSequence, outputComment))
								if len(outputSequence) != totalReadLength:
										self.log.critical("Read length exception: Expected %d, got %d (read start: %d, sequence length: %d)", totalReadLength, len(outputSequence), randomStart, len(sequence))
						else:
								outputReads.append(((outputSequence[0:readLength[0]], outputSequence[readLength[0] + readLength[1]:]), outputComment))
								if len(outputReads[-1][0]) != read1Length or len(outputReads[-1][0]) != read2Length:
										self.log.critical("Read length exception: Expected (%d:%d), got (%d:%d) (read start: %d, sequence length: %d", read1Length, read2Length, len(outputReads[-1][0]), len(outputReads[-1][1]), randomStart, len(sequence))
								
				return outputReads # end simulateRead()


		
		# Return statistics about this repertoire, suitable for saving to a file
    #
    # Arguments: none
    #
    # Returns:
    # Array with one element for each repertoire clone, formatted as so:
    # [ Clone count,
		#   V allele 1, J allele 1, CDR3 sequence 1, RNA sequence 1, DNA sequence 1,
		#   V allele 2, J allele 2, CDR3 sequence 2, RNA sequence 2, DNA sequence 2 ]
    #
    #
		
		def getStatistics(self, addHeader = False):
				retval = []
				if addHeader == True:
						retval.append(["CLONE,CELL_COUNT,VALLELE_1, JALLELE_1, CDR3_1, RNA_1, DNA_1, VALLELE_2, J_ALLELE_2, CDR3_2, RNA_2, DNA_2"])

				for i in range(0, len(self.repertoire)):
						CDR3_1, CDR3_2 = self.repertoire[i].getCDR3Sequences()

						v1Allele = "%s*%s" % (self.config.receptorSegment[self.repertoire[i].V1[0]]['gene'], self.repertoire[i].V1[1])
						v2Allele = "%s*%s" % (self.config.receptorSegment[self.repertoire[i].V2[0]]['gene'], self.repertoire[i].V2[1])
						j1Allele = "%s*%s" % (self.config.receptorSegment[self.repertoire[i].J1[0]]['gene'], self.repertoire[i].J1[1])
						j2Allele = "%s*%s" % (self.config.receptorSegment[self.repertoire[i].J2[0]]['gene'], self.repertoire[i].J2[1])
				
						stats = [ i, self.population[i],
											v1Allele, j1Allele, CDR3_1, self.repertoire[i].RNA1[3], self.repertoire[i].DNA1[3],
											v2Allele, j2Allele, CDR3_2, self.repertoire[i].RNA2[3], self.repertoire[i].DNA2[3] ]
											
						retval.append(stats)
				return retval
				
