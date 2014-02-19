#!/usr/bin/env python

# Copyright 2014 Frederick Kinyua Kamanu, Ph.D [frederick dot kamanu at gmail.com]
# All Rights Reserved.

#*****************************************************************************************
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>.
#*****************************************************************************************

import sys
import os
import random
import string
import re
import shutil
from optparse import OptionParser

#Main class setup
#-----------------------------------------------------------------------------------------
class GetGenome(object):
	"""
	A class for downloading fasta formated dna sequences given the accession number
	"""
	def __init__(self,acc_id,output_dir,email):
		"""
		initialize the acc id and email
		"""
		self.acc_id = acc_id
		self.output_dir = output_dir
		self.email = email
	
	def DownloadData(self):
		"""
		Download and save sequence by provided ID
		"""
		from Bio import Entrez
		Entrez.email = self.email
		handle = Entrez.efetch(db='nucleotide',id=self.acc_id,rettype='fasta',retmode='text')
		local_file=open(self.output_dir + '/' + self.acc_id, 'w')
		local_file.write(handle.read())
		handle.close()
		local_file.close()
		

#Parsing commandline options
#-----------------------------------------------------------------------------------------
def commandline_options():

    parser = OptionParser(usage="usage: %prog [options]")
    parser.add_option("-i", "--ids",
                      type="string",
                      dest="accession_list",
                      help="A list of accesion Ids to download")
    parser.add_option("-o", "--output",
                      type="string",
                      dest="output_directory",
                      help="Output directory for downloaded sequences")
    parser.add_option("-e", "--email",
                      type="string",
                      dest="email_address",
                      help="An email address that is needed by the NCBI")
                      
    (options, args) = parser.parse_args()
    options_args_parser = [options,args,parser]

    return options_args_parser

#Main
#-----------------------------------------------------------------------------------------
def main():

    options, args, parser = commandline_options()
    id_list = options.accession_list
    output_dir = options.output_directory
    email_address = options.email_address
    if id_list is None or email_address is None or output_dir is None:
        sys.stderr.write("\nError: A mandatory option is missing!\n")
        parser.print_help()
        sys.exit(-1)

    else:
    	try:
    		os.mkdir(output_dir)
    	except OSError:
    		shutil.rmtree(output_dir)
    		os.mkdir(output_dir)
    		
    	ids = re.split('[\s+,]',id_list)
    	for id in ids:
    		sys.stderr.write("downloading %s \n" % (id))
	    	GetGenome(id,output_dir,email_address).DownloadData()


#Eexecute main
#-----------------------------------------------------------------------------------------
if __name__ == "__main__":
    main() 
