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



