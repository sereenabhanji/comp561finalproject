"""
This class creates a Mutation object, which holds the rsid, variant sequence,
and information about the variant's location on the CFTR gene
"""


class Mutation:
    def __init__(self, rsid, start, end, sequence):
        self.rsid = rsid
        self.start = start
        self.end = end
        self.sequence = sequence


