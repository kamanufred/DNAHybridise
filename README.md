# DNAHybridise
 A parallel computational framework that simulates DNA-DNA hybridisation in calculating distance matrices, and inferring microbial phylogeny using whole-genome DNA sequences.   
   

#Requirements
- Unix / Mac OS 
- [Python](http://www.python.org/) (tested with Python 2.6/2.7)   
 Python Modules ([mpi4py](https://code.google.com/p/mpi4py/), [itertools](http://docs.python.org/2/library/itertools.html), [optparse](http://docs.python.org/2/library/optparse.html), [subprocess](http://docs.python.org/2/library/subprocess.html))
- [R](http://www.r-project.org/)   
R Packages ([ape](http://cran.r-project.org/web/packages/ape/), [getopt](http://cran.r-project.org/web/packages/getopt/index.html))
- [MUMmer](http://mummer.sourceforge.net/)  
 
#Installation
\- Install and configure a Message Passing Interface (MPI) Library e.g [Open MPI](http://www.open-mpi.org/)      
\- Checkout the source: `git clone https://github.com/kamanufred/DNAHybridise.git`   
\- This avails the 3 scripts, `GetGenome.py`, `DNAHybridise.py` and `GenomeCluster.R`  which can be ran independently.  
     
#Usage

###1. Calculating Genomic Distances
The `DNAHybridise.py` is a parallel program that computes distance matrices when given a set of genomes as input.  

To compute the similarity matrix for a set of bacterial genomes contained in a directory named 'genomes' using 10 cores, run the command;  
   
`mpirun -n 10 python DNAHybridise -g /path/to/genomes -o matrix.txt -c 100`

 Run `python DNAHybridise.py -h` to see all the options.    

    Usage: DNAHybridise.py [options]

    Options:
     -h, --help            show this help message and exit
     -g GENOME_PATH, --genomes=GENOME_PATH
                        A path to the list of genomes for analysis
     -o SIMILARITY_FILE, --output=SIMILARITY_FILE
                        The similarity matrix file (default matrix.txt)
     -c NUCMER_C, --mincluster=NUCMER_C
                        Minimum cluster size for nucmer (default 100)


###2. Inferring Phylogeny
The `GenomeCluster.R` script infers phylogeny using the similarity matrix which is output in step 1 above.

To infer phylogeny using the manhattan distance and a bootstrap of 500 replicates, run the command;  

   `Rscript GenomeCluster.R -i matrix.txt -d manhattan -b 500 -o results -m pdf `

 Run `Rscript GenomeCluster.R -h` to see all the options.  
 
     Usage: GenomeCluster.R [options]
	  --help,h		Print usage options
	  --input,i		Give the input file name (e.g matrix.txt)
	  --dist_method,d		Give the distance method (default manhattan)
				Other distance methods are: [maximum, manhattan, canberra, binary, minkowski]
	  --bootstrap,b		Give the bootstrap value (default 1000)
	  --output,o		Give output image file name
	  --image_format,m	Give output image file format (default PDF)

  

###3. Automatic Genomic Sequence Download
Genome sequences can be  downloaded automatically from GenBank/RefSeq database using the `GetGenome.py` script.   

__Example usage__: `python GetGenome.py -i a,b,c,d -o genomes -e myemail@email.com`
where `a`,`b`,`c` and `d` represent sample RefSeq/GenBank Accession IDs
    
   Run ` python GetGenome.py -h` to see all the options    

    Usage: GetGenome.py [options]

    Options:
     -h, --help            show this help message and exit
     -i ACCESSION_LIST, --ids=ACCESSION_LIST
                           A list of accesion Ids to download
     -o OUTPUT_DIRECTORY, --output=OUTPUT_DIRECTORY
                           Output directory for downloaded sequences
     -e EMAIL_ADDRESS, --email=EMAIL_ADDRESS
                           An email address that is needed by the NCBI   


#Running Tests    
From the test directory, run `python testDNAHybridise.py`   

#Citation
A manuscript detailing this work has been submitted for publication.   

#Contact
Contact `frederick(dot)kamanu(at)gmail.com` for feedback regarding the software.  

#License    
This software is provided under the GNU General Public License. See the included file.