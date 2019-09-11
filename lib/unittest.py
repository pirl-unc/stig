#! /usr/bin/python3

import stigtools
import unittest
import tempfile
import os
import pprint
import logging

myLog = logging.getLogger('main')
sh = logging.StreamHandler()
myLog.setLevel(99)
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
				vSegNum, vAllele =  self.config.chooseRandomSegment('A', componentName = 'V')
				self.assertTrue(self.config.receptorSegment[vSegNum]['segment_type'] == 'V')
				jSegNum, jAllele = self.config.chooseRandomSegment('A', componentName = 'J', V=(vSegNum, '01'))
				self.assertTrue(self.config.receptorSegment[jSegNum]['segment_type'] == 'J')

		def test_beta_random_1(self):
				vSegNum, vAllele =  self.config.chooseRandomSegment('B', componentName = 'V')

				# Must pick a random V
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
				for i in range(100):
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
				for i in range(100):
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
				self.config.VDJprobability = [('TRBV20-1', 1.00),
																			('TRBV20-1', 'TRBD1', 1.00),
																			('TRBV20-1', 'TRBD1', 'TRBJ2-1', 1.00),
																			]

				vSegNum, vAllele = self.config.chooseRandomSegment('B', componentName = 'V')
				dSegNum, dAllele = self.config.chooseRandomSegment('B', componentName = 'D', V=(vSegNum, '01'))
				jSegNum, jAllele = self.config.chooseRandomSegment('B', componentName = 'J', V=(vSegNum, '01'), D=(dSegNum, '01'))
				cSegNum, cAllele = self.config.chooseRandomSegment('B', componentName = 'C', V=(vSegNum, '01'), D=(dSegNum, '01'), J=(jSegNum, '01'))
				# WE'RE WORKING HERE!


class TestTcr(unittest.TestCase):
		def setUp(self):
				self.config = stigtools.tcrConfig()
				#self.tcr = stigtools.tcr(1.0, self.config)


				
if __name__ == '__main__':
		unittest.main()
