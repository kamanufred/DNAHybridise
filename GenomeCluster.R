#!/usr/bin/env Rscript

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

library(ape)
library(getopt)


#Parse commandline options
#-----------------------------------------------------------------------------------------
get_options <- function()
{
	spec = matrix(c(
		'help' , 'h', 0, "logical",
		'input' , 'i', 1, "character",
		'output' , 'o', 1, "character",
		'dist_method' , 'd', 1, "character",
		'bootstrap', 'b', 1, "integer",
		'image_format', 'm', 1, "character"
		), byrow=TRUE, ncol=4)
}	

#Create square matrix from an upper triangle matrix
#-----------------------------------------------------------------------------------------
create_square_m <- function(datafile)
{
    mat <- data.matrix(read.csv(datafile))
    mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
    return (mat)
}

#Draw heatmap
#-----------------------------------------------------------------------------------------
draw_heatmap <- function(inputfile,outfile,image_format)
{
	data <- create_square_m(inputfile)
	genome_matrix <- 1-cor(data)
	
	if(image_format == 'pdf')
	{
		pdf(file=paste(outfile,'_heatmap.pdf',sep=""))
	}
	else if (image_format == 'png')
	{
		png(file=paste(outfile,'_heatmap.png',sep=""))
	}
	else
	{
		writeLines("\n Wrong image format. Please specify pdf or png\n")
		q(status=1)
	}
	
	heatmap(genome_matrix,cexRow=0.7,cexCol=0.7,margins=c(13,13))

}

#Infer the nj tree	
#-----------------------------------------------------------------------------------------
infer_tree <- function(inputfile,dist_method,bootstrap,outfile,image_format)
{
	data <- create_square_m(inputfile)
	my_dist <- dist(data, method = dist_method)
	tree <- ladderize(nj(my_dist))
	if(image_format == 'pdf')
	{
		pdf(file=paste(outfile,'_tree.pdf',sep=""))
	}
	else if (image_format == 'png')
	{
		png(file=paste(outfile,'_tree.png',sep=""))
	}
	else
	{
		writeLines("\n Wrong image format. Please specify pdf or png\n")
		q(status=1)
	}
	plot(tree,cex=0.55,use.edge.length=FALSE)
	bp <- boot.phylo(tree, data, FUN=function(xx) nj(dist(xx, method = dist_method)), B=bootstrap,rooted=F)
	nodelabels(bp,adj = c(1.2, -0.2), frame = "n", cex = 0.4)
	tree$node.label<-bp
	write.tree(tree,paste(outfile,'_tree.nwk',sep=""),append = FALSE,digits = 2, tree.names = FALSE)

}


#Main
#-----------------------------------------------------------------------------------------
main <- function()
{
	spec <- get_options()
	options <-  getopt(spec)
	
	usage_information <- function()
	{
		writeLines("Usage: GenomeCluster.R [options]")
		writeLines("\t--help,h\t\tPrint usage options")
		writeLines("\t--input,i\t\tGive the input file name (e.g matrix.txt)")
		writeLines("\t--dist_method,d\t\tGive the distance method (default manhattan)")
		writeLines("\t\t\t\tOther distance methods are: [maximum, manhattan, canberra, binary, minkowski]")
		writeLines("\t--bootstrap,b\t\tGive the bootstrap value (default 1000)")
		writeLines("\t--output,o\t\tGive output image file name")
		writeLines("\t--image_format,m\tGive output image file format (default PDF)")
	}
	
	if (!is.null(options$help)) 
	{
		usage_information()
		q(status=1)
	}
	else if (is.null(options$input) && is.null(options$output))
	{
		writeLines("\nError!: Missing a mandatory option!\n")
		usage_information()
		q(status=1)
	}

	#Preset some defaults
	if ( is.null(options$dist_method) ) { options$dist_method = 'manhattan' }
	if ( is.null(options$bootstrap ) ) { options$bootstrap = 500 }
	if ( is.null(options$image_format) ) { options$image_format = 'pdf' }
	
	draw_heatmap(options$input,options$output,options$image_format)
	
	infer_tree(options$input,options$dist_method,options$bootstrap,
		   options$output,options$image_format)
		   
}

#Execute main
#-----------------------------------------------------------------------------------------
main()