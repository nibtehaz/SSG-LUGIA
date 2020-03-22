"""
    Codes to handle file operation
"""

from Bio import SeqIO

def loadGenome(file_path):
    """
    Reads a fasta file and loads the genome sequence
    
    Arguments:
        file_path {string} -- path to the fasta file
    
    Returns:
        string -- the genome sequence
    """

    genome = SeqIO.read(file_path,'fasta')
    
    return str(genome.seq)

