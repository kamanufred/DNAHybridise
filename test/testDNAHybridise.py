#!/usr/bin/env python

import sys
import os
import glob
import unittest
sys.path.append('..')
import DNAHybridise

class TestGetGenomeFunctions(unittest.TestCase):
	
	def setUp(self):
		self.data_path = os.path.abspath("genomes/")
		self.ref = os.path.abspath("genomes/sample1.fasta")
		self.query = os.path.abspath("genomes/sample2.fasta")
		self.ref_size = 10740
		self.query_size = 10740
		self.nucmer_c = 100
		self.wrapper_object = DNAHybridise.NUCmerWrapper(self.ref,self.query,self.ref_size,self.query_size,self.nucmer_c)
		self.parsedata_object = DNAHybridise.ParseData(self.data_path)

	def test_which(self):
		self.pythonpath = 'python'
		self.assertEqual(DNAHybridise.Which(self.pythonpath), 1)
		
	def test_GetFastaSize(self):
		self.sample_fasta = os.path.abspath("genomes/sample1.fasta")
		self.assertEqual(DNAHybridise.GetFastaSize(self.sample_fasta), 10740)
	
	def test_GetUniq(self):
		sample_list = [1,1,3,4,4]
		self.assertEqual(self.wrapper_object.GetUniq(sample_list), [1,3,4])
	
	def test_RunNUCmer(self):
		self.assertTrue(self.wrapper_object.RunNUCmer())
	
	def test_ParseNUCmer(self):
		self.assertEqual(self.wrapper_object.ParseNUCmer(), 0.029992821697200467)
		all_temp = glob.glob('temp_*')		
		for i in all_temp:
			os.system("rm -rf %s" % (i))
			
	def test_GetData(self):
		self.assertEqual(len(self.parsedata_object.GetData()),3)
		
	def test_CreatePairs(self):
		self.assertEqual(len(self.parsedata_object.CreatePairs()),3)

	
if __name__ == '__main__':
    unittest.main()

    