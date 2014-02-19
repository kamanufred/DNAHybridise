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
import subprocess
import re
import itertools
from optparse import OptionParser

from mpi4py import MPI

#DNAHybridise main class
#-----------------------------------------------------------------------------------------
class NUCmerWrapper(object):
    """
    A NUCmer wrapper class for computing pairwise genome similarity
    between a query and the reference genome. Returns a similarity score
    """

    def __init__(self,ref,query,ref_size,query_size,nucmer_c):
        """
        Intialize the user-defined parameters for NUCmer analysis
        """
        self.ref = ref
        self.query = query
        self.ref_size = ref_size
        self.query_size = query_size
        self.nucmer_c = nucmer_c

    def GetUniq(self,seq):
        """
        generate a unique list
        """
        seen = set()
        seen_add = seen.add
        return [ i for i in seq if i not in seen and not seen_add(i)]

    def RunNUCmer(self):
        """
        Run a NUCmer pairwise comparison.
        """
        #Create a random directory name
        lst = [random.choice(string.ascii_letters + string.digits) for n in xrange(10)]
        random_string = "".join(lst)
        random_dir_name = "temp_" + random_string
        os.mkdir(random_dir_name)
        clean_ref = self.ref.split("/")[-1]
        clean_query =  self.query.split("/")[-1]


        #execute NUCmer
        os.system("cp %s %s %s" % (self.ref,self.query,random_dir_name))
        os.system("cd %s; nucmer -maxmatch -c %d %s %s 2>nucmer.log; \
                show-coords -rclT out.delta >nucmer.coords; rm out.delta" % \
                (random_dir_name,self.nucmer_c,clean_ref,clean_query ))
        return random_dir_name

    def ParseNUCmer(self):
        """
        Parse the NUCmer coords file and calculate the similarity score
        """
        nucmer_dir = self.RunNUCmer()
        nucmer_file = nucmer_dir + '/nucmer.coords'
        all_files = nucmer_dir + '/*'
        f = open(nucmer_file, "r")
        contig_and_size = dict()
        contig_start_end = []

        for line in f:
            split_line = (line.strip().split('\t'))
            if len(split_line) == 13:
                start_pos = int(split_line[0])
                end_pos = int(split_line[1])
                contig_id = split_line[11]
                contig_size = int(split_line[7])
                similarity = float(split_line[6])

                contig_and_size[contig_id] = contig_size
                contig_start_end.append((contig_id,start_pos,end_pos))

        total_match_size = 0
        for contig in contig_and_size:
            match_range = []
            for coords in contig_start_end:
                if coords[0] == contig:
                    synteny = range(coords[1],coords[2]+1)
                    for i in synteny:
                        match_range.append(i)

            uniq_match_range = self.GetUniq(match_range)
            match_size = len(uniq_match_range)
            total_match_size += match_size

        match_proportion1 = float(total_match_size) / float(self.ref_size)
        match_proportion2 = float(total_match_size) / float(self.query_size)
        
        match_p_mult = (match_proportion1 * match_proportion2)
        os.system("rm -rf %s" % (all_files)) #clean up
        os.system("rm -rf %s" % (nucmer_dir)) #clean up
        
        return match_p_mult

#Parse object class
#-----------------------------------------------------------------------------------------
class ParseData(object):
    """
    Extract and organise genomic data into pairs for similarity
    analyis
    """

    def __init__(self,data_path):

        """
        Initialize the genomic data path
        """
        self.data_path = data_path

    def GetData(self):
        """
        Get genomic sequences locally
        """
        species = []
        formated_path = self.data_path.rstrip('/')
        genomes = os.listdir(formated_path)

        if len(genomes) == 0:
            sys.exit("\nError!: The provided directory is empty\n")
        else:
            species = [formated_path + '/' + genome for genome in genomes]
            return species

    def CreatePairs(self):
        """
        Create genome pairs for a pair-wise comparison
        """
        genomes = self.GetData()
        pair_list = [",".join(map(str,comb)) for comb in itertools.combinations(genomes, 2)]
        return pair_list


#Dependency check
#-----------------------------------------------------------------------------------------
def Which(program):
    """
    Check if a program has been installed in a *NIX system and is in the path
    """
    status = 0
    try:
        # pipe output to /dev/null for silence
        null = open("/dev/null", "w")
        subprocess.Popen(program, stdout=null, stderr=null)
        null.close()
        status = 1

    except OSError:
        status = 0

    return status


#Fasta file size
#-----------------------------------------------------------------------------------------
def GetFastaSize(in_file):
    """
    Determine genome size and check for empty or improperly formated
    fasta files
    """
    size_count = 0
    f = open(in_file, "r")
    start_characters = dict()
    empty_line = re.compile('^$')
    for line in f:
        if not re.match(empty_line,line):
            stripped_line = line.strip()
            start_character = stripped_line[0]
            start_characters[start_character] = ''

            if not stripped_line.startswith('>'):
                line_size = len(stripped_line)
                size_count += line_size
    if '>' not in start_characters or size_count == 0:
        sys.exit("\nError!: Empty or improperly formated fasta file encountered\n")
    else:
        return size_count


#Write output matrix
#-----------------------------------------------------------------------------------------
def CreateMatrix(scores,output_file):
    """
    Read pair-wise scores and tabulate them into an upper triangular matrix for n genomes
    """
    #get uniq names from scores input
    uniq_names = sorted(list(set([i[0] for i in scores] + [i[1] for i in scores])))

    data_d = dict()
    count = 1
    for line in scores:
        ref,query,score = line
        data_d[count] = {}
        data_d[count]['ref'] = ref
        data_d[count]['query'] = query
        data_d[count]['score'] = score
        count += 1
    
    #create a matrix structure
    f_out = open(output_file,'w')
    header = ",".join(uniq_names)
    f_out.write("%s\n" % (header))
    for i in uniq_names:
        pairs = []
        for j in uniq_names:
            score = ''
            for key in data_d:
                if data_d[key]['ref'] == i and data_d[key]['query'] == j:
                    score = str(data_d[key]['score'])
            if i == j:
                pairs.append('1')
            else:
                pairs.append(score)
        pair_string = ','.join(pairs)                
        f_out.write("%s,%s\n" % (i,pair_string))
        

#Setting commandline options
#-----------------------------------------------------------------------------------------
def commandline_options():
    """
    Commandline arguments input
    """
    parser = OptionParser(usage="usage: %prog [options]")
    parser.add_option("-g", "--genomes",
                      type="string",
                      dest="genome_path",
                      help="A path to the list of genomes for analysis")
    parser.add_option("-o", "--output",
                      action="store",
                      type="string",
                      dest="similarity_file",
                      default="matrix.txt",
                      help="The similarity matrix file (default matrix.txt)")
    parser.add_option("-c","--mincluster",
                      type="int",
                      dest="nucmer_c",
                      default=100,
                      help="Minimum cluster size for nucmer (default 100)")

    (options, args) = parser.parse_args()
    options_args_parser = [options,args,parser]

    return options_args_parser


#Processing genome pairs
#-----------------------------------------------------------------------------------------
def ProcessGenomePairs(ref_query_nucmer_c):
    """"
    Analyse genomes in pairwise fashion and return the similarity score
    between them
    """
    ref,query,nucmer_c = ref_query_nucmer_c
    ref_size = GetFastaSize(ref)
    query_size = GetFastaSize(query)

    wrapper_object = NUCmerWrapper(ref,query,ref_size,query_size,nucmer_c)
    raw_score = wrapper_object.ParseNUCmer()

    clean_ref = ref.split('/')[-1]
    clean_query = query.split('/')[-1]
    return (clean_ref,clean_query,raw_score)


#Main
#-----------------------------------------------------------------------------------------
def main():

    options, args, parser = commandline_options()
    genome_file_path = options.genome_path
    if genome_file_path is None:
        sys.stderr.write("\nError: A mandatory option is missing!\n")
        parser.print_help()
        sys.exit(-1)
    nucmer_c = options.nucmer_c
    matrix_file = options.similarity_file
   
    mummer_status = Which('mummer')
    if mummer_status == 0:
        sys.exit("\nError!: MUMmer is not installed\n")
    else:
        comm = MPI.COMM_WORLD
        size = comm.Get_size()
        rank = comm.Get_rank()
        
        pair_list = ParseData(genome_file_path).CreatePairs()
        ref_query_nucmer_c = []
        for i in pair_list:
            ref,query = i.split(',')
            ref_query_nucmer_c.append((ref,query,nucmer_c))
        
        if rank == 0:
            sys.stderr.write("\nStarting analysis. Number of processes: %d\n" % (size))
            data = ref_query_nucmer_c
            pieces = [[] for _ in range(size)]
            for i, piece in enumerate(data):
                pieces[i % size].append(piece)
        else:
            data = None
            pieces = None
        data = comm.scatter(pieces, root=0)
        rank_results = []
        for i in data:
            rank_results.append(ProcessGenomePairs(i))
        
        comm.Barrier()
        all_rank_results = comm.gather(rank_results,root=0)
        if rank == 0:
            final_data = [i for i in itertools.chain.from_iterable(all_rank_results)]
            CreateMatrix(final_data,matrix_file)
            sys.stderr.write("\nAnalysis complete, similarity matrix written to %s\n" % (matrix_file))
            

#Execute main
#-----------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()

