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
    

    

