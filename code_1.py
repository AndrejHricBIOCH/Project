{
 "cells": [],
 "metadata": {},
 "nbformat": 4,
 "nbformat_minor": 4
}

def read_fasta_seq():
    """
    Import Biopython modules and read fasta files
    """
    #import Biopython modules
    import Bio
    from Bio import SeqIO
    from Bio import Phylo
    from Bio.Phylo.TreeConstruction import DistanceMatrix
    from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

    # read the database and mystery sequence
    database = {record.id : str(record.seq) for record in SeqIO.parse("dog_breeds.fa", "fasta")}
    unknown = {record.id : str(record.seq) for record in SeqIO.parse("mystery.fa", "fasta")}
    

    def hamming_distance(seq1, seq2):
        """
        Determine hamming distance between two sequences i.e., number of mutations between two sequences
        """
        return sum([1 for x, y in zip(seq1, seq2) if x!=y])

def closest_sequence():
    """
    Comparing mystery sequence to dog_breed database and finiding the closest sequence
    """
    # take the mystery sequence from unknown dict
    mystery = unknown['gb|KM061522.1|']
    pairs = {}
    # iterate over all database sequences and compare with mystery sequence
    for Id, seq in database.items():
        pairs[(Id, seq)] = hamming_distance(seq, mystery)
        
    # find out the closest breed and print the result
    closest_breed = sorted(pairs.items(), key=lambda x: x[1])[0]
    Id, seq, mut = closest_breed[0][0], closest_breed[0][1], closest_breed[1]

    # print truncated sequence of the closest breed sequence () +sequence id
    return(f"The closest breed (mutation={mut}) to unknown sequence in the database is:\n{seq[:90]}...\t{Id}")

def seq_phylogeny():
    """
    Computation of a phylogenetic tree, via construction of a distance matrix and using phylo module
    """
    # extract ids and sequences from dict record
    ids = list(database.keys()) + list(unknown.keys())
    sequences = list(database.values()) + list(unknown.values())

    # construct distance matrix to draw phylogenetic tree with Bio.Phylo module
    distM = [[0]*i for i in range(1, len(sequences)+1)]
    for i, x in enumerate(sequences):
        for j, y in zip(range(i), sequences):
            if i != j:
                distM[i][j] = hamming_distance(x, y)
    
    # construct and draw the phylogenetic tree with upgma method
    dm = DistanceMatrix(ids, distM)
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(dm)
    return (Phylo.draw_ascii(tree))