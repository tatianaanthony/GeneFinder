# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Tatiana Anthony

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('T')
    'A'
    >>> get_complement('G')
    'C'
    """
    # TODO: implement this
    """Use if/else if/else to check which nucleotide it is.
    If it's not any of the 4, throw an eror
    """

    """

    """

    if nucleotide == 'T':
        return 'A'
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide=='G':
        return 'C'
    elif nucleotide == 'A':
        return 'T' 
    else:
        print(nucleotide + " isn't a nucleotide letter!")
        return 'N/A'


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    rc_dna='' #Initiate reverse complement dna
    for element in dna: #go through each element in the dna from front to back
        rc_dna = get_complement(element) + rc_dna #Add element to the front of the list (building from back to front)

    return rc_dna


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("CATGAATGTAGATAGATGTGCCC")
    'CATGAATGTAGA'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    # TODO: implement this
 
    i=0
    while i*3+3<=len(dna):
        codon = dna[3*i:3*i+3] #Identifies an in-frame codon
        if aa_table[codon] == '|': #Checks if the codon is a stop codon
            return dna[:3*i] #returns everything before the stop codon
        i=i+1 # go to the next codon
    return dna


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    i=0
    all_ORFs=[]
    while i*3+3<=len(dna):
        codon = dna[3*i:3*i+3] #Identifies an in-frame codon
        if aa_table[codon] == 'M': #Checks if the codon is a stop codon
            new_dna = dna[3*i+3:] #takes the rest of the DNA minus the start codon
            new_ORF = codon + rest_of_ORF(new_dna) #Adds the start codon to the ORF found
            all_ORFs.append(new_ORF) #Add ORF to list of ORFs

            num_codons = int(len(new_ORF)/3) #Figure out how many codons the ORF is
            # print(num_codons)
            # print(all_ORFs)
            i= i+num_codons #skip to the end of the ORF to start searching again
        else:
            i=i+1 #Go to the next codon
    return all_ORFs



def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    all_ORFs = []
    i=0
    while i<3: #Each start position
        test_dna = dna[i:] #Get rid of first element to move frame over 1
        all_ORFs = all_ORFs + find_all_ORFs_oneframe(test_dna) #Add ORFs from this frame to all_ORFs
        i+=1

    return all_ORFs


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    # TODO: implement this
    all_ORFs=find_all_ORFs(dna) #find ORFs of one strand
    all_ORFs = all_ORFs + find_all_ORFs(get_reverse_complement(dna)) #Find and add ORFs of reverse complement strand
    return all_ORFs


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    all_ORFs=find_all_ORFs_both_strands(dna)
    # print("__________________________________")
    # print(all_ORFs)
    # print("__________________________________")
    current_longest = 0
    for i in range(len(all_ORFs)):
        if len(all_ORFs[i])>current_longest:
            current_longest = i
    return all_ORFs[current_longest]



def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    # TODO: implement this
    pass


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    # TODO: implement this
    pass


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # TODO: implement this
    pass

if __name__ == "__main__":
    import doctest
    doctest.run_docstring_examples(longest_ORF, globals(),verbose=True)
