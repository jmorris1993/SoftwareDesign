# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 09:10:26 2014

@author: julian
"""


""""
HW 3 PART 1
"""




"""
Found this instead of the amino_acid dictionary from
http://zientzilaria.herokuapp.com/blog/2007/05/24/translating-dna-into-proteins/
"""
aminoacids = { 'ATA':'I',
    'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T',
    'ACT':'T', 'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S',
    'AGT':'S', 'AGA':'R', 'AGG':'R', 'CTA':'L', 'CTC':'L', 'CTG':'L',
    'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 'CAC':'H',
    'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R',
    'CGT':'R', 'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A',
    'GCC':'A', 'GCG':'A', 'GCT':'A', 'GAC':'D', 'GAT':'D', 'GAA':'E',
    'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 'TCA':'S',
    'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L',
    'TTG':'L', 'TAC':'Y', 'TAT':'Y', 'TAA':'\_', 'TAG':'\_', 'TGC':'C',
    'TGT':'C', 'TGA':'_', 'TGG':'W', }
              

"""the join function can do what collapse does"""

def collapse(L):
    """ Converts a list of strings to a string by concatenating all elements of the list """
    output = ""
    for s in L:
        output = output + s
    return output


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents a protein coding region).
        
        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment
    """

    protein = ''  #makes an empty string called protein

    for n in range(0,len(dna),3):                #loops from 0 to the length of DNA sequence in steps of 3
        if aminoacids.has_key(dna[n:n+3]) == True:  #checks to make sure dna is in triplets
            protein += aminoacids[dna[n:n+3]]       #adds the corresponding protein to existing list
    return protein
    
#print coding_strand_to_AA('AGCCGT')
    


def coding_strand_to_AA_unit_tests():
    """ Unit tests for the coding_strand_to_AA function """
    
    print "input:"+str("ATGCGA")
    print "expected output:" + str('MR') 
    print "actual output:" + str(coding_strand_to_AA("ATGCGA"))
    

   
#coding_strand_to_AA_unit_tests()


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence
    
        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    """
    """http://zientzilaria.herokuapp.com/blog/2008/03/11/
    fasta-module-generating-reverse-complement-of-dna-sequences/
    I got the idea for how to do this from this website, but I modified it
    so that it all fits into one function instead of two.
    """

    dna = list(dna)    #turns dna into a list
    dna.reverse()      #reverses the list
    reverse=''.join(dna)      #joins the list into a string
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}  #dictionary of complements
    complemented = [complement[base] for base in reverse]  
    #I used list comprehension because people online said it's the most efficient way.
    #It replaces a letter in the reversed list with the assigned value 
    #in the dictionary for every letter.
    return ''.join(complemented)   #joins the list back into a string

#print get_reverse_complement("ATCGGG")
    
def get_reverse_complement_unit_tests(dna):
    """ Unit tests for the get_complement function """
        
    print "input:"+str(dna)
    print "expected output:" + str('AAAGCGGGCAT') 
    print "actual output:" + str(get_reverse_complement(dna))
    
#get_reverse_complement_unit_tests('ATGCCCGCTTT')


"""
HW 3 PART 2
"""


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start codon and returns
        the sequence up to but not including the first in frame stop codon.  If there
        is no in frame stop codon, returns the whole string.
        
        dna: a DNA sequence
        returns: the open reading frame represented as a string
    """
  
    result = ""

    for i in range(0,len(dna)):
        if i%3==0 and dna[i:i+3]!="TAA" and dna[i:i+3]!="TAG" and dna[i:i+3]!="TGA":
            result=result+str(dna[i:i+3])
        elif i%3==0 and dna[i:i+3]=="TAA":
            break
        elif i%3==0 and dna[i:i+3]=="TAG":
            break
        elif i%3==0 and dna[i:i+3]=="TGA":
            break
    
    return result

#print rest_of_ORF("ATGAGATAGG")

def rest_of_ORF_unit_tests(dna):
    """ Unit tests for the rest_of_ORF function """
    print "input:"+str(dna)
    print "expected output:" + "ATGAGA"
    print "actual output:" + str(rest_of_ORF(dna))
    
#rest_of_ORF_unit_tests("ATGAGATAGG")
        
           
def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence and returns
        them as a list.  This function should only find ORFs that are in the default
        frame of the sequence (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.
        
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    """

    openings=[]
    for x in range(0,len(dna)-2,3):
        if dna[x:x+3]=="ATG":
            openings.append(rest_of_ORF(dna[x:]))
    
    return openings

#print find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")

     
def find_all_ORFs_oneframe_unit_tests(dna):
    """ Unit tests for the find_all_ORFs_oneframe function """

    print "input:"+str(dna)
    print "expected output:" + str(['ATGCATGAATGTAGA', 'ATGTGCCC'])
    print "actual output:" + str(find_all_ORFs_oneframe(dna))

#find_all_ORFs_oneframe_unit_tests("ATGCATGAATGTAGATAGATGTGCCC")

def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in all 3
        possible frames and returns them as a list.  By non-nested we mean that if an
        ORF occurs entirely within another ORF and they are both in the same frame,
        it should not be included in the returned list of ORFs.
        
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    """
     
    openings=[]
    for x in range(0,len(dna)):
        if dna[x:x+3]=="ATG":
            openings.append(rest_of_ORF(dna[x:]))
    return openings

#print find_all_ORFs("ATGCATGAATGTAG")


def find_all_ORFs_unit_tests():
    """ Unit tests for the find_all_ORFs function """
    print "input:"+str("ATGCATGAATGTAG")
    print "expected output:" + str(['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG'])
    print "actual output:" + str(find_all_ORFs("ATGCATGAATGTAG"))
    
#find_all_ORFs_unit_tests()


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.
        
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    """
    print "ORF for original strand:"+str(find_all_ORFs(dna))
    print "ORF for reverse complement strand:"+str(find_all_ORFs(get_reverse_complement(dna)))

#find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")

def find_all_ORFs_both_strands_unit_tests():
    """ Unit tests for the find_all_ORFs_both_strands function """

    print "input:"+str("ATGCGAATGTAGCATCAAA")
    print "expected output:" 
    print str("ORF for original strand:['ATGCGAATG', 'ATG']")
    print str("ORF for reverse complement strand:['ATGCTACATTCGCAT']")
    print "actual output:" 
    print str(find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA"))

#find_all_ORFs_both_strands_unit_tests()


"""
HW 3 PART 3
"""



def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string"""    
        
    """found idea for code from http://stackoverflow.com/questions/873327/
    pythons-most-efficient-way-to-choose-longest-string-in-list"""
        
    bothORF= find_all_ORFs(dna)+find_all_ORFs(get_reverse_complement(dna))
    print max(bothORF, key=len)
    
#longest_ORF("ATGCGAATGTAGCATCAAA")


def longest_ORF_unit_tests():
    """ Unit tests for the longest_ORF function """

    print "input:"+str("ATGCGAATGTAGCATCAAA")
    print
    print "expected output:" 
    print str("ATGCTACATTCGCAT") 
    print
    print "actual output:" 
    print str(longest_ORF("ATGCGAATGTAGCATCAAA"))
    
#longest_ORF_unit_tests()


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence
        
        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    total=[]
    while num_trials>0:
        a=list(dna)              #converts dna to list
        import random
        random.shuffle(a)
        stringed=''.join(a)      #used the built-in function join instead of the collapse function
        b=find_all_ORFs(stringed)
        
        if len(b)>0:             #to prevent error
            c=max(b, key=len)    #finds the string with maximum length
            total.append(c)
        num_trials=num_trials-1
    longestORF=max(total, key=len)        #finds the string with maximum length from all the trials
    return longestORF

    
    
#print longest_ORF_noncoding("ATGCGAATGTAGCATCAAA",10)


def gene_finder(dna, threshold):
    """returns ORF for original strand"""
    originalORF=find_all_ORFs(dna)
    originalORF=[x for x in originalORF if len(x)>=threshold]
    
    """returns ORF for reverse complement strand"""
    reversecompORF=find_all_ORFs(get_reverse_complement(dna))
    reversecompORF=[x for x in reversecompORF if len(x)>=threshold]
    

    proteinoriginal=[]
    proteinreversecomp=[]
    
    for x in originalORF:
        a=coding_strand_to_AA(x)
        proteinoriginal.append(a)
    for x in reversecompORF:
        b=coding_strand_to_AA(x)
        proteinreversecomp.append(b)
    print "protein sequence for OriginalORF:"    
    print proteinoriginal
    print "protein sequence for reversecompORF:"
    print proteinreversecomp
    
    
    
#gene_finder("ATGCGAATGTAATAGGGGATGTTAGTGACTGCATGCGATGCAAA",5)
    

"""
Next, use your longest_ORF_noncoding on the Salmonella DNA sequence and 
compute a conservative threshold for distinguishing between genes and 
non-genes by running longest_ORF_noncoding for 1500 trials.  Make a note of 
this number as it will be used with the gene_finder function.

Next, use your gene_finder function with the original Salmonella DNA sequence 
and the threshold computed above to get a list of candidate genes.

Finally, take the amino acid sequences produced by your gene finder and 
search for them using protein-BLAST.  What types of genes appear to be in 
this DNA sequence?  Record your findings in a file called salmonella.txt.
"""

"""
import load
dna = load.load_seq("./data/X73525.fa")
"""

import load
dna=load.load_seq("./X73525.fa")
counter=-1
current=dna[counter]
while current=="C"or current=="A"or current=="G" or current=="T":
    counter-=1
    current=dna[counter]
dna=dna[counter+1:]



print len(longest_ORF_noncoding(dna, 1500))
"""it gave back 639"""



#print gene_finder(dna,639)

    
    
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    