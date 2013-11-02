#!usr/bin/python
# Filename = translate.py

f = file('nucleotide.fa')
#sequence = f.readline()

tranalate_code = {    'AAA' : 'K', 'AAC' : 'N', 'AAG' : 'K', 'AAT' : 'N',
                      'ACA' : 'T', 'ACC' : 'T', 'ACG' : 'T', 'ACT' : 'T',
                      'AGA' : 'R', 'AGC' : 'S', 'AGG' : 'R', 'AGT' : 'S',
                      'ATA' : 'I', 'ATC' : 'I', 'ATG' : 'M', 'ATT' : 'I',
                      'CAA' : 'Q', 'CAC' : 'H', 'CAG' : 'Q', 'CAT' : 'H',
                      'CCA' : 'P', 'CCC' : 'P', 'CCG' : 'P', 'CCT' : 'P',
                      'CGA' : 'R', 'CGC' : 'R', 'CGG' : 'R', 'CGT' : 'R',
                      'CTA' : 'L', 'CTC' : 'L', 'CTG' : 'L', 'CTT' : 'L',
                      'GAA' : 'E', 'GAC' : 'D', 'GAG' : 'E', 'GAT' : 'D',
                      'GCA' : 'A', 'GCC' : 'A', 'GCG' : 'A', 'GCT' : 'A',
                      'GGA' : 'G', 'GGC' : 'G', 'GGG' : 'G', 'GGT' : 'G', 
                      'GTA' : 'V', 'GTC' : 'V', 'GTG' : 'V', 'GTT' : 'V', 
                      'TAA' : '*', 'TAC' : 'Y', 'TAG' : '*', 'TAT' : 'Y', 
                      'TCA' : 'S', 'TCC' : 'S', 'TCG' : 'S', 'TCT' : 'S', 
                      'TGA' : '*', 'TGC' : 'C', 'TGG' : 'W', 'TGT' : 'C', 
                      'TTA' : 'L', 'TTC' : 'F', 'TTG' : 'L', 'TTT' : 'F'
}

def translate():
    """Translate a nucleotide sequence into protein sequence"""

    global sequence
    n = len(sequence) / 3
    sequence_protein = ''
    for i in range(0, n):
        code = sequence[i * 3 : i * 3 + 3]
        sequence_protein = sequence_protein + tranalate_code[code]
    return sequence_protein

def choose_orf(c = 1):
    """This function choose the ORF for translate function

    c is the code for ORF choosing, c can be (+1, +2, +3, -1, -2, -3, 6)
    c = +1 for translate in the forward sequence from 1st nucleotide
    c = +2 for translate in the forward sequence from 2nd nucleotide
    c = +3 for translate in the forward sequence from 3rd nucleotide
    c = -1 for translate in the backward sequence from 1st nucleotide
    c = -2 for translate in the backward sequence from 2nd nucleotide
    c = -3 for translate in the backward sequence from 3rd nucleotide
    c = 6 for translate in all the 6 conditions above"""

    global sequence
    sequence = f.readline()
    print c
    if c == 1 or c == 2 or c ==3:
        sequence = sequence[c - 1 :]
        return translate()
print choose_orf(2)
print choose_orf(2)
print choose_orf(3)
print choose_orf()