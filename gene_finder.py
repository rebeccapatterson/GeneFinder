# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE:
GENE FINDER CODE AS OF 02-07-16.  HAS ALL FUNCTIONALITY TO MY KNOWLEDGE.

@author: REBECCA PATTERSON

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
    """
    if nucleotide == 'A' :
       return 'T'
    elif nucleotide == 'T' :
       return 'A'
    elif nucleotide == 'C' :
        return 'G'
    elif nucleotide == 'G' :
        return 'C'
    else :
        return 'not a DNA nucleotide'


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
    reverse_complement= ''
    i= len(dna)-1
    while i >=0 :
        reverse_complement+=(get_complement(dna[i]))
        i= i-1
    return reverse_complement	
        

def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """ 
    i=0
    while i<len(dna) :
    	if dna[i:i+3]== 'TAA' or dna[i:i+3]== 'TAG' or dna[i:i+3]== 'TGA':
    		return dna[:i]
    	else :
    	    i+=3
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
    all_frame_orfs= []
    i=0
    while i<len(dna):
        if dna[i:i+3]== 'ATG' :
            all_frame_orfs.append(rest_of_ORF(dna[i:]))
            i= i+len(rest_of_ORF(dna[i:]))
        else:
            i+=3
    return all_frame_orfs    


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
    all_orfs= []
    i=0
    while i<3 :
        all_orfs.extend(find_all_ORFs_oneframe(dna[i:]))
        i+=1
    return all_orfs


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    all_both_orfs=[]
    reverse= get_reverse_complement(dna)
    all_both_orfs.extend(find_all_ORFs(dna))
    all_both_orfs.extend(find_all_ORFs(reverse))
    return all_both_orfs


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    longest=''
    i=0
    orfs= find_all_ORFs_both_strands(dna)
    while i<len(orfs) :
        if orfs[i]> len(longest) :
            longest= orfs[i]
        i+=1
    return longest        


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    longest_non= ''
    i=0
    while i<num_trials :
        shuffled_orf= longest_ORF(shuffle_string(dna))
        if len(shuffled_orf)>len(longest_non) :
            longest_non= shuffled_orf
        i+=1
    return len(longest_non)        


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
    amino= ''
    i=0
    while i+3<len(dna) :
        a_acid= aa_table[dna[i:i+3]]
        amino+= a_acid
        i+=3
    return amino    


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    print 'finding threshold'
    threshold= longest_ORF_noncoding(dna,1500)
    print threshold
    print 'finding all orfs'
    all_orfs= find_all_ORFs_both_strands(dna)
    i=0
    amino_sequences=[]
    print 'starting while loop'
    while i<len(all_orfs) :
        if len(all_orfs[i])>threshold :
            amino= coding_strand_to_AA(all_orfs[i])
            amino_sequences.append(amino)
        i+=1
    print 'found all amino sequences'
    print amino_sequences    
    

#    import doctest
#    doctest.run_docstring_examples(get_complement,             globals(), verbose=True)
#    doctest.run_docstring_examples(get_reverse_complement,     globals(), verbose=True)
#    doctest.run_docstring_examples(rest_of_ORF,                globals(), verbose=True)
#    doctest.run_docstring_examples(find_all_ORFs_oneframe,     globals(), verbose=True)
#    doctest.run_docstring_examples(find_all_ORFs,              globals(), verbose=True)
#    doctest.run_docstring_examples(find_all_ORFs_both_strands, globals(), verbose=True)
#    doctest.run_docstring_examples(longest_ORF,                globals(), verbose=True)
#    doctest.run_docstring_examples(coding_strand_to_AA,        globals(), verbose=True)



gene_finder(load_seq("./data/X73525.fa"))