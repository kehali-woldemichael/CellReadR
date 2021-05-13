library(biomaRt)
library(Biostrings)

basePath <- '/home/user1/Dropbox/Research/Neurobiology_PhD/Rotations/Huang/Projects/CellReadR/Code/Output/biomaRt/'

generate_refSeq <- function(geneName, species, typeEnsembl) {
    # Just renaming to ensure that right species found 
    if(species == 'Mouse') {
        searchName <- "Mouse genes"
    } else {
        searchName <- species
    }

    totalPath <- paste(basePath, species, '/Reverse_', geneName, '.fasta', sep='')

    # Only if does not already exist ... otherwise just extends exists ... resulting in repeats in fasta file
    if(file.exists(totalPath) == FALSE){ 
        # Get references sequences 
        referenceSeq <- fetchSequences(geneName, searchName, typeEnsembl)
        # Save reference sequences 
        saveFasta(geneName, species, referenceSeq@seq, referenceSeq@seqCDS, referenceSeq@seqCDNA)
    }
}

# Class for storing sequences ... 
# Created mostly because R does not allow functions to have multiple outputs 
setClass(Class="ReferenceSeq", 
         representation(
             seq='data.frame', 
             seqCDS='data.frame',
             seqCDNA='data.frame'
         )
)

# Just in case need to check ensembl dataset 
getDataset <- function(geneName) {
    ensembl <- useMart("ensembl")
    datasets <- listDatasets(ensembl)

    return(datasets) 
}

# Fetching sequences from Ensembl
fetchSequences <- function(geneName, species, typeEnsembl){
    ensembl <- useMart("ensembl")
    datasets <- listDatasets(ensembl)
    
    nameDataset <- datasets[grep(species, datasets$description), ]$dataset
    ensembl <- useDataset(nameDataset, mart = ensembl)
    
    # Getting sequences of exons for gene 
    seq <- getSequence(id = geneName, 
                      type = typeEnsembl, 
                      seqType = "gene_exon", 
                      mart = ensembl)
    
    # Getting sequences of exons for gene 
    seqCDS <- getSequence(id = geneName, 
                      type = typeEnsembl, 
                      seqType = "coding", 
                      mart = ensembl)
    
   # Getting sequences of exons for gene 
    seqCDNA <- getSequence(id = geneName, 
                      type = typeEnsembl, 
                      seqType = "cdna", 
                      mart = ensembl) 
    
    if(nrow(seq) == 0) {
        stop('Sequence not fetched')
    }
    else {
        # Returning ReferenceSeq object 
        return(new("ReferenceSeq", seq=seq, seqCDS=seqCDS, seqCDNA=seqCDNA))
    }
}

# Generating array of sequences tbat are reverse complement of reference array of sequences 
convert_reverseComplement <- function(arraySeq) {
    temp_seq <- arraySeq 
    for(i in 1:nrow(arraySeq)) {
        temp_seq[i,1] <- toString(reverseComplement(DNAString(arraySeq[i,1])))
    }
    return(temp_seq)
}


# Generating array of complement sequences to array of reference sequence 
convert_complement <- function(arraySeq) {
    temp_seq <- arraySeq 
    for(i in 1:nrow(arraySeq)) {
        temp_seq[i,1] <- toString(complement(DNAString(arraySeq[i,1])))
    }
    return(temp_seq)
}

# Converting to RNA ...
convert_RNA <- function(arraySeq) {
    temp_seq <- arraySeq 
    for(i in 1:nrow(arraySeq)) {
        temp_seq[i,1] <- toString(RNAString(DNAString(arraySeq[i,1])))
    }
    return(temp_seq)
}

# Function for saving Fasta files 
saveFasta <- function(geneName, species, seq, seqCDS, seqCDNA) {
    #saveBaseName <- paste(getwd(), '/Output/biomaRt/', species, '/', sep = '' )
    saveBaseName <- paste(basePath, species, '/', sep='')
    
    # Create Ouput directory if it does not exist in current working directory 
    if(dir.exists("Output") == FALSE) {
        dir.create("Output")
    }
    # Create biomaRt directory if it does not exist 
    if(dir.exists("Output/biomaRt") == FALSE) {
        dir.create("Output/biomaRt")
    }
    # Create directory for species if it does not exist
    if(dir.exists(saveBaseName) == FALSE) {
        dir.create(saveBaseName)
    }
    
    # Saving exon sequences if greater than 200 bp
    passedSeq <- seq[nchar(seq$gene_exon) > 200, ]
    fileName <- paste(saveBaseName, 'Template_', geneName, '.fasta', sep = '')
    exportFASTA(passedSeq, file = fileName)
    
    # Generating reverse complement FASTA file  
    reverseSeq <- convert_reverseComplement(passedSeq)
    reverse_fileName <- paste(saveBaseName , 'Reverse_', geneName, '.fasta', sep = '')
    exportFASTA(reverseSeq, file = reverse_fileName)
    
    # Generating complement FASTA file  
    complementSeq <- convert_complement(passedSeq)
    complement_fileName <- paste(saveBaseName, 'Complement_', geneName, '.fasta', sep = '')
    exportFASTA(complementSeq, file = complement_fileName)
    
    # Generating CDS FASTA file  
    cds_fileName <- paste(saveBaseName, 'CDS_', geneName, '.fasta', sep = '')
    exportFASTA(seqCDS, file = cds_fileName)
    
    # Generating RNA CDS RNA FASTA file  
    cds_RNA_seq <- convert_RNA(seqCDS)
    cds_RNA_fileName <- paste(saveBaseName, 'CDS_RNA_', geneName, '.fasta', sep = '')
    exportFASTA(cds_RNA_seq, file = cds_fileName)

    # Generating cDNA FASTA file  
    cDNA_fileName <- paste(saveBaseName, 'cDNA_', geneName, '.fasta', sep = '')
    exportFASTA(seqCDNA, file = cDNA_fileName)

    # Generating RNA cDNA FASTA file  
    cDNA_RNA_seq <- convert_RNA(seqCDNA)
    cDNA_RNA_fileName <- paste(saveBaseName, 'cDNA_RNA_', geneName, '.fasta', sep = '')
    exportFASTA(cDNA_RNA_seq, file = cDNA_RNA_fileName)

}
