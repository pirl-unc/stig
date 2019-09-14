#! /usr/bin/python3

import stigtools
import unittest
import tempfile
import os
import re
import pprint
import logging
import difflib # NOCOMMIT

config_iterations = 100
config_iterations = 1 # NOCOMMIT

myLog = logging.getLogger('main')
sh = logging.StreamHandler()
#myLog.setLevel(logging.DEBUG) # NOCOMMIT
myLog.setLevel(logging.WARNING)
#sh.setFormatter(logging.Formatter(fmt='%(asctime)s.%(msecs)03d [%(levelname)s] %(name)s %(message)s',
#																	datefmt='%Y%m%d%H%M%S'))
myLog.addHandler(sh)

class TestTcrConfig(unittest.TestCase):

		def setUp(self):
				(self.tempfilehandle, self.tempfilename) = tempfile.mkstemp()
				self.config = stigtools.tcrConfig()
        						
		def tearDown(self):
				os.close(self.tempfilehandle)
				os.remove(self.tempfilename)
				
		def test_config_file_not_exist(self):
				with self.assertRaises(IOError):
						self.config.readTCRConfig('/STIG_invalid_filename_does_not_exist')

		def test_config_file_invalid_format(self):
				os.write(self.tempfilehandle, str.encode("Invalid\ttest\tcomponents"))
				os.lseek(self.tempfilehandle, 0, os.SEEK_SET)
				with self.assertRaises(ValueError):
						self.config.readTCRConfig(self.tempfilename)

class TestTcrConfig_chooseRandomSegment(unittest.TestCase):

		def setUp(self):
				self.config = stigtools.tcrConfig(log = myLog.getChild('tcrConfig'))
				self.config.setWorkingDir('./data')


		def test_alpha_random_1(self):
				# Pick a random TRAV
				vSegNum, vAllele =  self.config.chooseRandomSegment('A', componentName = 'V')
				self.assertTrue(self.config.receptorSegment[vSegNum]['segment_type'] == 'V')
				self.assertTrue(self.config.receptorSegment[vSegNum]['region'] == 'V-REGION')

				# Pick a random TRAJ
				jSegNum, jAllele = self.config.chooseRandomSegment('A', componentName = 'J', V=(vSegNum, '01'))
				self.assertTrue(self.config.receptorSegment[jSegNum]['segment_type'] == 'J')
				self.assertTrue(self.config.receptorSegment[jSegNum]['region'] == 'J-REGION')

		def test_beta_random_1(self):

				# Must pick a random V
				vSegNum, vAllele =  self.config.chooseRandomSegment('B', componentName = 'V')
				self.assertTrue(self.config.receptorSegment[vSegNum]['segment_type'] == 'V')

				# Requests for TRBJ without a defined TRBD fail
				with self.assertRaises(ValueError):
						self.config.chooseRandomSegment('B', componentName = 'J', V=(vSegNum, '01'))

				# Must pick a random TRBD
				dSegNum, jAllele = self.config.chooseRandomSegment('B', componentName = 'D', V=(vSegNum, '01'))
				self.assertTrue(self.config.receptorSegment[dSegNum]['segment_type'] == 'D')

				# Must pick a random TRBJ				
				jSegNum, jAllele = self.config.chooseRandomSegment('B', componentName = 'J', V=(vSegNum, '01'), D=(dSegNum, '01'))
				self.assertTrue(self.config.receptorSegment[jSegNum]['segment_type'] == 'J')

		def test_alpha_random_repeated(self):
				for i in range(config_iterations):
						vSegNum, vAllele = self.config.chooseRandomSegment('A', componentName = 'V')
						jSegNum, jAllele = self.config.chooseRandomSegment('A', componentName = 'J', V=(vSegNum, '01'))
						cSegNum, cAllele = self.config.chooseRandomSegment('A', componentName = 'C', V=(vSegNum, '01'), J=(jSegNum, '01'))

						vSeg = self.config.receptorSegment[vSegNum]
						jSeg = self.config.receptorSegment[jSegNum]
						cSeg = self.config.receptorSegment[cSegNum]

						# Receptor segment types match those requested
						self.assertTrue(vSeg['segment_type'] == 'V')
						self.assertTrue(jSeg['segment_type'] == 'J')
						self.assertTrue(cSeg['segment_type'] == 'C')

						# Receptor segments are on the same chromosome
						self.assertTrue(vSeg['chromosome'] == jSeg['chromosome'] == cSeg['chromosome'])

						# There is no requirement for segments to be on the same strand (see TRBV30, which is reverse strand)

		def test_beta_random_repeated(self):
				for i in range(config_iterations):
						vSegNum, vAllele = self.config.chooseRandomSegment('B', componentName = 'V')
						dSegNum, dAllele = self.config.chooseRandomSegment('B', componentName = 'D', V=(vSegNum, '01'))
						jSegNum, jAllele = self.config.chooseRandomSegment('B', componentName = 'J', V=(vSegNum, '01'), D=(dSegNum, '01'))
						cSegNum, cAllele = self.config.chooseRandomSegment('B', componentName = 'C', V=(vSegNum, '01'), D=(dSegNum, '01'), J=(jSegNum, '01'))

						vSeg = self.config.receptorSegment[vSegNum]
						dSeg = self.config.receptorSegment[dSegNum]
						jSeg = self.config.receptorSegment[jSegNum]
						cSeg = self.config.receptorSegment[cSegNum]

						# Receptor segment types match those requested
						self.assertTrue(vSeg['segment_type'] == 'V')
						self.assertTrue(dSeg['segment_type'] == 'D')
						self.assertTrue(jSeg['segment_type'] == 'J')
						self.assertTrue(cSeg['segment_type'] == 'C')

						# Receptor segments are on the same chromosome
						self.assertTrue(vSeg['chromosome'] == dSeg['chromosome'] == jSeg['chromosome'] == cSeg['chromosome'])

						# There is no requirement for segments to be on the same strand (see TRBV30, which is reverse strand)


		def test_beta_nonrandom(self):
				VDJprobability_orig = 	self.config.VDJprobability
				self.config.VDJprobability = [('TRBV20-1', 1.00),
																			('TRBV20-1', 'TRBD1', 1.00),
																			('TRBV20-1', 'TRBD1', 'TRBJ2-1', 1.00), ]

				vSegNum, vAllele = self.config.chooseRandomSegment('B', componentName = 'V')
				dSegNum, dAllele = self.config.chooseRandomSegment('B', componentName = 'D', V=(vSegNum, '01'))
				jSegNum, jAllele = self.config.chooseRandomSegment('B', componentName = 'J', V=(vSegNum, '01'), D=(dSegNum, '01'))
				cSegNum, cAllele = self.config.chooseRandomSegment('B', componentName = 'C', V=(vSegNum, '01'), D=(dSegNum, '01'), J=(jSegNum, '01'))

				vSeg = self.config.receptorSegment[vSegNum]
				dSeg = self.config.receptorSegment[dSegNum]
				jSeg = self.config.receptorSegment[jSegNum]
				cSeg = self.config.receptorSegment[cSegNum]

				self.assertTrue(vSeg['gene'] == 'TRBV20-1')
				self.assertTrue(dSeg['gene'] == 'TRBD1')
				self.assertTrue(jSeg['gene'] == 'TRBJ2-1')
				self.assertTrue(cSeg['gene'] == 'TRBC2') # TRBC2 is only choice here due to our selection of TRBJ2-1

				self.config.VDJprobability = VDJprobability_orig


class TestTcrConfig_recombinate(unittest.TestCase):
		def setUp(self):
				self.config = stigtools.tcrConfig(log = myLog.getChild('tcrConfig'))
				self.config.setWorkingDir('./data')

		def test_known_beta(self):
				# Set some fixed recombination parameters to make our checking easier
				junctionProbability_orig = self.config.junctionProbability
				self.config.junctionProbability = {
						'D3chewback': [0.2, 0.2, 0.2, 0.2, 0.2],
						'D5chewback':[0.2, 0.2, 0.2, 0.2, 0.2],
						'DJaddition':[0.2, 0.2, 0.2, 0.2, 0.2],
						'Jchewback':[0.2, 0.2, 0.2, 0.2, 0.2],
						'VDaddition':[0.2, 0.2, 0.2, 0.2, 0.2],
						'VJaddition':[0.2, 0.2, 0.2, 0.2, 0.2],
						'Vchewback':[0.2, 0.2, 0.2, 0.2, 0.2],
						}
				
				# TRBV20-1 L-V-GENE-UNIT is contained in HG38 at NC_000007.14:142626727-142627438
				# The *01 allele differs at the lowercase positions here
				# This segment has been edited to remove the V-RS (heptamer, V-spacer and nonamer) tail, as would be done in preparation for inclusion into a functional TCR
				vSegDNA = "ATGCTGCTGCTTCTGCTGCTTCTGGGGCCAGGTATAAGCCTCCTTCTACCTGGGAGCTTGGGTGGGCATGTGCGTGTGTTGGCATGGTCAAGTGGTGGCCAGCAGGGTTGCAATGTGGATTGTTTATGCTCATCGAAGGGAGAGGGAGAGGCCCTGCTCTCTAGAGGTGTAAATGGTAAGGTGAAAGCCGCCGGTCAGAGGAGATGGGGGATTATGGCCCTAGGGAGATGACGGGAAGATTGCACAAAACAAACAGGACTCTCCAGGAGCTGGGAGCACAGGGAGGGAGTGAGGCTCAGCTCTGCCTGGCGTCCTGTCTGACTCGGCTCCCACTGGGCTCTCCTCTCTCTCTGGCTTCTGTCTCAGCAGGCTCCGGGCTTGGTGCTGTCGTCTCTCAACATCCGAGCtGGGTTATCTGTAAGAGTGGAACCTCTGTGAAGATCGAGTGCCGTTCCCTGGACTTTCAGGCCACAACTATGTTTTGGTATCGTCAGTTCCCGAAACAGAGTCTCATGCTGATGGCAACTTCCAATGAGGGCTCCAAGGCCACATACGAGCAAGGCGTCGAGAAGGACAAGTTTCTCATCAACCATGCAAGCCTGACCTTGTCCACTCTGACAGTGACCAGTGCCCATCCTGAAGACAGCAGCTTCTACATCTGCAGTGCTAGAGA"# V-RS: CACAGCGCCAGGAGGGGATCAGACACCGCGGCAAGAACC
				# TRBV20-1 L-PART1+V-EXON from IMGT
				vSegRNA = "atgctgctgcttctgctgcttctggggccaggtataagcctccttctacctgggagcttggcaggctccgggcttggtgctgtcgtctctcaacatccgagctgggttatctgtaagagtggaacctctgtgaagatcgagtgccgttccctggactttcaggccacaactatgttttggtatcgtcagttcccgaaacagagtctcatgctgatggcaacttccaatgagggctccaaggccacatacgagcaaggcgtcgagaaggacaagtttctcatcaaccatgcaagcctgaccttgtccactctgacagtgaccagtgcccatcctgaagacagcagcttctacatctgcagtgctagaga"
				dSeg = "gggactagcggggggg"
				jSeg = "ctcctacaatgagcagttcttcgggccagggacacggctcaccgtgctag"

				# TRBC2 EX1-4 from NCBI NC_000007.14:142801040-142802526
				# n.b. this sequence doesn't exactly match STIG's alleles since the first G here is duplicated between the J and C segments (which IMGT does to preserve full codons in their allele files)
				cSegDNA = "AGGACCTGAAAAACGTGTTCCCACCCAAGGTCGCTGTGTTTGAGCCATCAGAAGCAGAGATCTCCCACACCCAAAAGGCCACACTGGTGTGCCTGGCCACAGGCTTCTACCCCGACCACGTGGAGCTGAGCTGGTGGGTGAATGGGAAGGAGGTGCACAGTGGGGTCAGCACAGACCCGCAGCCCCTCAAGGAGCAGCCCGCCCTCAATGACTCCAGATACTGCCTGAGCAGCCGCCTGAGGGTCTCGGCCACCTTCTGGCAGAACCCCCGCAACCACTTCCGCTGTCAAGTCCAGTTCTACGGGCTCTCGGAGAATGACGAGTGGACCCAGGATAGGGCCAAACCTGTCACCCAGATCGTCAGCGCCGAGGCCTGGGGTAGAGCAGGTGAGTGGGGCCTGGGGAGATGCCTGGAGGAGATTAGGTGAGACCAGCTACCAGGGAAAATGGAAAGATCCAGGTAGCGGACAAGACTAGATCCAGAAGAAAGCCAGAGTGGACAAGGTGGGATGATCAAGGTTCACAGGGTCAGCAAAGCACGGTGTGCACTTCCCCCACCAAGAAGCATAGAGGCTGAATGGAGCACCTCAAGCTCATTCTTCCTTCAGATCCTGACACCTTAGAGCTAAGCTTTCAAGTCTCCCTGAGGACCAGCCATACAGCTCAGCATCTGAGTGGTGTGCATCCCATTCTCTTCTGGGGTCCTGGTTTCCTAAGATCATAGTGACCACTTCGCTGGCACTGGAGCAGCATGAGGGAGACAGAACCAGGGCTATCAAAGGAGGCTGACTTTGTACTATCTGATATGCATGTGTTTGTGGCCTGTGAGTCTGTGATGTAAGGCTCAATGTCCTTACAAAGCAGCATTCTCTCATCCATTTTTCTTCCCCTGTTTTCTTTCAGACTGTGGCTTCACCTCCGGTAAGTGAGTCTCTCCTTTTTCTCTCTATCTTTCGCCGTCTCTGCTCTCGAACCAGGGCATGGAGAATCCACGGACACAGGGGTGTGAGGGAGGCCAGAGCCACCTGTGCACAGGTACCTACATGCTCTGTTCTTGTCAACAGAGTCTTACCAGCAAGGGGTCCTGTCTGCCACCATCCTCTATGAGATCTTGCTAGGGAAGGCCACCTTGTATGCCGTGCTGGTCAGTGCCCTCGTGCTGATGGCCATGGTAAGGAGGAGGGTGGGATAGGGCAGATGATGGGGGCAGGGGATGGAACATCACACATGGGCATAAAGGAATCTCAGAGCCAGAGCACAGCCTAATATATCCTATCACCTCAATGAAACCATAATGAAGCCAGACTGGGGAGAAAATGCAGGGAATATCACAGAATGCATCATGGGAGGATGGAGACAACCAGCGAGCCCTACTCAAATTAGGCCTCAGAGCCCGCCTCCCCTGCCCTACTCCTGCTGTGCCATAGCCCCTGAAACCCTGAAAATGTTCTCTCTTCCACAGGTCAAGAGAAAGGATTCCAGAGGC"
				# TRBC2 EX1-4 from IMGT
				cSegRNA = "aggacctgaaaaacgtgttcccacccgaggtcgctgtgtttgagccatcagaagcagagatctcccacacccaaaaggccacactggtgtgcctggccacaggcttctaccccgaccacgtggagctgagctggtgggtgaatgggaaggaggtgcacagtggggtcagcacagacccgcagcccctcaaggagcagcccgccctcaatgactccagatactgcctgagcagccgcctgagggtctcggccaccttctggcagaacccccgcaaccacttccgctgtcaagtccagttctacgggctctcggagaatgacgagtggacccaggatagggccaaacctgtcacccagatcgtcagcgccgaggcctggggtagagcagactgtggcttcacctccgagtcttaccagcaaggggtcctgtctgccaccatcctctatgagatcttgctagggaaggccaccttgtatgccgtgctggtcagtgccctcgtgctgatggccatggtcaagagaaaggattccagaggc"

				# Chromosome 7: 142796415-142801039, based on IMGT coordinates for end of TRBJ2-1 and TRBC2
				jcSegDNA = "GTAAGAAGGGGGCTCCAGGTGGGAGAGAGGGTGAGCAGCCCAGCCTGCACGACCCCAGAACCCTGTTCTTAGGGGAGTGGACACTGGGCAATCCAGGGCCCTCCTCGAGGGAAGCGGGGTTTGCGCCAGGGTCCCCAGGGCTGTGCGAACACCGGGGAGCTGTTTTTTGGAGAAGGCTCTAGGCTGACCGTACTGGGTAAGGAGGCGGTTGGGGCTCCGGAGAGCTCCGAGAGGGCGGGATGGGCAGAGGTAAGCAGCTGCCCCACTCTGAGAGGGGCTGTGCTGAGAGGCGCTGCTGGGCGTCTGGGCGGAGGACTCCTGGTTCTGGGTGCTGGGAGAGCGATGGGGCTCTCAGCGGTGGGAAGGACCCGAGCTGAGTCTGGGACAGCAGAGCGGGCAGCACCGGTTTTTGTCCTGGGCCTCCAGGCTGTGAGCACAGATACGCAGTATTTTGGCCCAGGCACCCGGCTGACAGTGCTCGGTAAGCGGGGGCTCCCGCTGAAGCCCCGGAACTGGGGAGGGGGCGCCCCGGGACGCCGGGGGCGTCGCAGGGCCAGTTTCTGTGCCGCGTCTCGGGGCTGTGAGCCAAAAACATTCAGTACTTCGGCGCCGGGACCCGGCTCTCAGTGCTGGGTAAGCTGGGGCCGCCGGGGGACCGGGGACGAGACTGCGCTCGGGTTTTTGTGCGGGGCTCGGGGGCCGTGACCAAGAGACCCAGTACTTCGGGCCAGGCACGCGGCTCCTGGTGCTCGGTGAGCGCGGGCTGCTGGGGCGCGGGCGCGGGCGGCTTGGGTCTGGTTTTTGCGGGGAGTCCCCGGGCTGTGCTCTGGGGCCAACGTCCTGACTTTCGGGGCCGGCAGCAGGCTGACCGTGCTGGGTGAGTTTTCGCGGGACCACCCGGGCGGCGGGATTCAGGTGGAAGGCGGCGGCTGCTTCGCGGCACCCGGTCCGGCCCTGTGCTGGGAGACCTGGGCTGGGTCCCCAGGGTGGGCAGGAGCTCGGGGAGCCTTAGAGGTTTGCATGCGGGGGTGCACCTCCGTGCTCCTACGAGCAGTACTTCGGGCCGGGCACCAGGCTCACGGTCACAGGTGAGATTCGGGCGTCTCCCCACCTTCCAGCCCCTCGGTCCCCGGAGTCGGAGGGTGGACCGGAGCTGGAGGAGCTGGGTGTCCGGGGTCAGCTCTGCAAGGTCACCTCCCCGCTCCTGGGGAAAGACTGGGGAAGAGGGAGGGGGTGGGGAGGTGCTCAGAGTCCGGAAAGCTGAGCAGAGGGCGAGGCCACTTTTAATCTTTTTTCTGGGGTGTTTAGAGAGAAGGTGAACGATGGAGGAGAGGATTTGTTAGGACTCTGGGAGAGGCGAGACTGGAGAGGACGAAGGGAAATCCTGGTTTGGGGAATGGGTAGGAGTGGGGGTAACTGCTATTCGTAGGCAAAAAGAGCTGAGCAGGCTGGGAACAGCGCGGGTGGGCAAGGGTCAGCACTGCGGGCAGGCGGGTGGGTGTTAGGGGGCAGAAATCCTGCAGCCGAGGGTGCAGTAGAACACAGAAGAAAAAGCCTGCCAAACAAAAGTGGAACAGAGAAGCCAAAAAGGGAGATGAACATGAGTCAGTGAAGAAAAGAATGAAAGTTTACTGTTTAGCAGTGTGGATCTCTAATCCGACTTAAAACTCCTTGTTCCCGATTCCTATTCCTCCTAAGCCAGAGATCCCTGGGTCCAGGGTGAGGGCACGGCATTCATGCTTACCCACGGGCTGGTCAACAAAGAGGTGCTGACCTGAGAGTAGGGCACATAACCTCAGCCACTGGGGTACACTTACCACCCCCGCCCCCGTGTAGCTCCCTCCCCTATCCTGAAATCTCCCTTAGCACACTAAGTATTCTAGGTTAAACAGCCCAGATGTTCAGGGAGTTCATTCGCCACAAACACACATTAAAATGCAGACAATTTGCCTGTGAGATGAGGAAAATTCTCTGGAAGATTTAGGCCCTGAGAGCTGAAAAGGGACCCTAAACATTACCTGGTGACAACTGCCCTGAGGCCAGAGAAGAGAACTCACAATATTGGTATATTAACCGGTACCATTTGTAGTTAGGCTGTCATTAATCTGGGTGTAATGGGGCTCAGCTACAGAGAAGCGTATCCAGAGGAAATGTGGGGTTCCTGCAGTCAGCTGGGGCACCGAAAAGACCCAGACTTTAGAACCAGATAAAGGCTGAGTTCGAACCTCTGGTTCTTTTGTGATGTGGTGACCTTGGGCAAATTAATGTGTGAACCTCAGTTTCCTCAACTATAAAATGTAATGAACAATACCTACCACTTACTATTGCTGTGAGGAAGAAAAGAGAGTCAACATGTACCGTGTACAGATTATTGATGTAATTCAATGGTTCTTTTCCCCATCCTCCTAGGAGGTCACTGGGGAGACAGGGGGCAGGGTCAGCCCAGTGCCAAAGGATGGGCAGGATCTGAAGTGTGGAAATGGAGTAAGGCTGTGTCTGTGTCAGAGGTGGGTTGGGAAGATGTGAGACAAACATCACAATTTTGCCTAAGGTGAATCCAACCCACAAGTAGAGCACAGGCCAACAGCAGCTCACTAGTACACATACTTACACCAGCAGCTCACTAGTACACACACTTACACCAGACGCTCACTGGTACACACTCACACCAGAAGCTCATTAGAACACACACACCGTCAGCTTGCTAGTACACACACTTACACCAGAAGCTCACTAGTATACACTTACACCAGAAGCTCACTAGTACACACACTTACTCCAGAAGCTCACTAGTACACACACTTACACCCACAGAGACAAGCCCCACACCACACGGACTCACAAATGCAGAAGAAGAGTTGACCCTGCCCTCATGTACAGATAGTATGCTGGTGTGTAGTGAGAGGCGGCTAGTGTTCGCCACGCCTGTTGCATCAATAACTATTCCACAGAACAGCATACATGGTCAAGGAATATTTTTAATATTTACAATGAGGACTGAGTGTTTGCTAAATGAAGAAATCAAAGGAGTACATTCTGGAGGGCTGGTAGACTCCCAGAGCCAGGGTTTTTGAGATGAGAGTGAAAATAAGCAGGCTGAGAGCAGAAAGAAAGAAAGTTCTGCAGATAGAAGTTAGGATATTACCACTTCGGCCCCAGCCCAGCCAGTAAATGTTTAGAAGCATAGTAGTAATTAGCAGGTAGGAGTTGTGGGGAGGAAAGGAAACTGACATGATGGGAAGCAGGACCAGCTACAAAATCTTCAGGACCAGCTCAGAATGAAAATGCCGGGGCTCGTGTTAAATCTTAGGATTTCAACATGGTGACAGCAGAGCACTAAACCAAGCAGGGGGCCCTTCTAGCAGGAGGCCCCATGGGGACCTAGGAAGGAGGTCACATCCATACAAATGTAAAGAAATGTCAGAGTAAATTTCCTCCATTGTCCAGGAGATTCGGAATAGGTTCTCCTAAGACTGATATTTCTTCATTTTAATAGAGTTGCTCAGAAATGAAAAACAATCAATGGGAAGAAAAAGAAAAAAAAAAAAGATAAGATGAGTAGGAGGGCAGGTCCAAAAGAAGTGATTCAGCAAAATGAAAGGGGTCCTCAGGGATTAAAGGGGATGAATTTACCTGTCATCCCTAAGAATCTACAAAGGAGATGCTCAGGACAGAAACTGTATCAACACAACTAGTAGCAAGAAGTTACTCTGATGATATCAGATGTTTATTTGGGAAACTTGCTAGTAGAGAAAGCTACATATAATATTTGGATGCAAAGGGACACAGAAGGTTGAAGAGTCCCTAATTTTGAAATAAGGGAAGATGACTAACTGTCTGAGCTGAGAAAACTCAGGGGTACCTGGAGGCAGAGGAATGGATAAGATGACTTCATGCACCACAAAAAGAAAAAACCTCACATTCTCATGAACGCACTGTAAAACCAAAGGATGTCCTCATGTGAATGCAAAAAATAGGCCATCTGTAAATCCAAAGAAAGCCCCCAGATCTAAAATGTCTCCCTCATCCCAGATTCCCCTTCATTCCTGAGCACCTTAGATTTGGTATAAATAACCTGCTTGGGAGGGGGCTTTTTGAATTCGTACATAATTTAACCTTCACACAGTTTCTGCAAAGTCAGAATGGTGATTATTACCTCACATGCAGAAAAAAGTGATAGGAATTTCTGTCTTAAAAGTCTTGTTGGTGGACAAAGGAAGTTCTAGGATTTGGATCTTGTTTTTTTGGGTTCCAATCCCTTGCTCCAGTTAAAAAACTACCACATAAAATGGTGAGAAGTAGGTAGGCAAGTTTTTATTGATAGAGAGGAAATCAAATAATGGCAATGAGGAGACATCACCTGGAATGTTAGGCAGTGCCTAACTGGGGGATGGACAGACAATGGGCAGTGCCAACCCATAGGGTGGATACAAAAGACAGGCAAGGAAGGGGTAGAACCATCAAAGAGGAATAGGCTGGTGACCCCAAAGCAAGGAGGACCTAGTAACATAATTGTGCTTCATTATGGTCCTTTCCCGGCCTTCTCTCTCACACATACACAGAGCCCCTACCAGGACCAGACAGCTCTTAGAGCAACCCTAGCCCCATTACCTCTTCCCTTTCCAG"
				
				# Locate the index values of some known gene segments
				vSegNum = dSegNum = jSegNum = cSegNum = None
				for i in range(0, len(self.config.receptorSegment)):
						if re.match('^[VDJ]-REGION|EX1', self.config.receptorSegment[i]['region']):
								if self.config.receptorSegment[i]['gene'] == 'TRBV20-1':
										vSegNum = i
								elif self.config.receptorSegment[i]['gene'] == 'TRBD2':
										dSegNum = i
								elif self.config.receptorSegment[i]['gene'] == 'TRBJ2-1':
										jSegNum = i
								elif self.config.receptorSegment[i]['gene'] == 'TRBC2':
										cSegNum = i

				# Make sure we found all of them
				self.assertTrue(vSegNum is not None)
				self.assertTrue(dSegNum is not None)
				self.assertTrue(jSegNum is not None)
				self.assertTrue(cSegNum is not None)

				# Recombinate may return None if our recombination fails (early stops, frameshift, etc), so we loop here:
				for i in range(0, 500):
						myLog.setLevel(logging.DEBUG) # NOCOMMIT
						x = self.config.recombinate((vSegNum, '01'), (dSegNum, '01'), (jSegNum, '01'), (cSegNum, '01'))
						if x is not None:
								DNA = x[0][3]
								RNA = x[1][3]
								dnaRegex = ("^" +
														vSegDNA[:-4] + # Segment
														"((" + ")|(".join([vSegDNA[-4:], vSegDNA[-4:-1], vSegDNA[-4:-2], vSegDNA[-4:-3]]) + "))?" + # Chewback
														"[CTAG]{,4}" + # insertions
														"((" + ")|(".join([dSeg[:4], dSeg[1:4], dSeg[2:4], dSeg[3:4]]) + "))?" +
														dSeg[4:-4] + 
														"((" + ")|(".join([dSeg[-4:], dSeg[-4:-1], dSeg[-4:-2], dSeg[-4:-3]]) + "))?" +
														"[CTAG]{,4}" +
														"((" + ")|(".join([jSeg[:4], jSeg[1:4], jSeg[2:4], jSeg[3:4]]) + "))?" +
														jSeg[4:] +
														jcSegDNA +
														cSegDNA +
														"$" ).upper()
								dnaRegex2 = ("^ATG" +
#														 vSegDNA[:-4] +
														 
														 ".*$").upper()

								dnaRegex3 = (vSegDNA[:-4] + # Segment
														"((" + ")|(".join([vSegDNA[-4:], vSegDNA[-4:-1], vSegDNA[-4:-2], vSegDNA[-4:-3]]) + "))?" + # Chewback
														"[CTAG]{,4}" + # insertions
														"((" + ")|(".join([dSeg[:4], dSeg[1:4], dSeg[2:4], dSeg[3:4]]) + "))?" +
														dSeg[4:-4] + 
														"((" + ")|(".join([dSeg[-4:], dSeg[-4:-1], dSeg[-4:-2], dSeg[-4:-3]]) + "))?" ).upper()

								dnaRegex4 = (vSegDNA[:-4] + # Segment
														"((" + ")|(".join([vSegDNA[-4:], vSegDNA[-4:-1], vSegDNA[-4:-2], vSegDNA[-4:-3]]) + "))?" + # Chewback
														"[CTAG]{,4}" + # insertions
														"((" + ")|(".join([dSeg[:4], dSeg[1:4], dSeg[2:4], dSeg[3:4]]) + "))?" +
														dSeg[4:-4] + 
														"((" + ")|(".join([dSeg[-4:], dSeg[-4:-1], dSeg[-4:-2], dSeg[-4:-3]]) + "))?" +
														"[CTAG]{,4}" +
														"((" + ")|(".join([jSeg[:4], jSeg[1:4], jSeg[2:4], jSeg[3:4]]) + "))?" +
														jSeg[4:]).upper()

								rnaRegex = ("^" +
														vSegRNA[:-4] + # Segment
														"((" + ")|(".join([vSegRNA[-4:], vSegRNA[-4:-1], vSegRNA[-4:-2], vSegRNA[-4:-3]]) + "))?" + # Chewback
														"[CTAG]{,4}" + # insertions
														"((" + ")|(".join([dSeg[:4], dSeg[1:4], dSeg[2:4], dSeg[3:4]]) + "))?" +
														dSeg[4:-4] + 
														"((" + ")|(".join([dSeg[-4:], dSeg[-4:-1], dSeg[-4:-2], dSeg[-4:-3]]) + "))?" +
														"[CTAG]{,4}" +
														"((" + ")|(".join([jSeg[:4], jSeg[1:4], jSeg[2:4], jSeg[3:4]]) + "))?" +
														jSeg[4:] +
														cSegRNA +
														"$").upper()
								rnaRegex2 = ("^" +
														vSegRNA[:-4] + # Segment
														"((" + ")|(".join([vSegRNA[-4:], vSegRNA[-4:-1], vSegRNA[-4:-2], vSegRNA[-4:-3]]) + "))?" + # Chewback
														".*").upper()
								rnaRegex3 = ("^" +
														vSegRNA[:-4] + # Segment
														"((" + ")|(".join([vSegRNA[-4:], vSegRNA[-4:-1], vSegRNA[-4:-2], vSegRNA[-4:-3]]) + "))?" + # Chewback
														"[CTAG]{,4}" + # insertions
														"((" + ")|(".join([dSeg[:4], dSeg[1:4], dSeg[2:4], dSeg[3:4]]) + "))?" +
														dSeg[4:-4] + 
														".*$").upper()
								rnaRegex4 = ("^" +
														vSegRNA[:-4] + # Segment
														"((" + ")|(".join([vSegRNA[-4:], vSegRNA[-4:-1], vSegRNA[-4:-2], vSegRNA[-4:-3]]) + "))?" + # Chewback
														"[CTAG]{,4}" + # insertions
														"((" + ")|(".join([dSeg[:4], dSeg[1:4], dSeg[2:4], dSeg[3:4]]) + "))?" +
														dSeg[4:-4] + 
														"((" + ")|(".join([dSeg[-4:], dSeg[-4:-1], dSeg[-4:-2], dSeg[-4:-3]]) + "))?" +
														"[CTAG]{,4}" +
														"((" + ")|(".join([jSeg[:4], jSeg[1:4], jSeg[2:4], jSeg[3:4]]) + "))?" +
														jSeg[4:] +
														".*$").upper()


								self.assertTrue(re.match(dnaRegex2, DNA))
								self.assertTrue(re.match(dnaRegex3, DNA))
								self.assertTrue(re.match(dnaRegex4, DNA))
								self.assertTrue(re.match(dnaRegex, DNA))
								self.assertTrue(re.match(rnaRegex, RNA))

								self.config.junctionProbability = junctionProbability_orig


				# If our loop runs out, fail the test with a message for user/developer
				self.assertTrue("Recombinate failed too many times, perhaps retry unit tests?" == '')


class TestTcr(unittest.TestCase):
		def setUp(self):
				self.config = stigtools.tcrConfig()
				#self.tcr = stigtools.tcr(1.0, self.config)


				
if __name__ == '__main__':
		unittest.main()
