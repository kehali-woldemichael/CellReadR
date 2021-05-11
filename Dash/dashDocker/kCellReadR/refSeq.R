#!/usr/bin/env Rscript

# Reading in paramaters for fetching references sequences 
args <- commandArgs(trailingOnly=TRUE)

# Generate reference sequeces if necessary 
path_functions <- 'kCellReadR/functions_biomaRt.R'
source(path_functions)
generate_refSeq(args[1], args[2], 'wikigene_name')
