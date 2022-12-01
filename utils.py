def fasta_reader(s):
    with open(s, 'r') as fájl:
        inpt = fájl.read().split('\n')
        dikt = {}
        prev = ''
        for s in inpt:
            if len(s) == 0:
                break
            if s[0] == '>':
                s = s[1:]
                dikt[s] = ''
                prev = s
            elif len(s) != 0:
                dikt[prev] += s
        return dikt

def gc_content(s):
    return len([i for i in s if i == 'G' or i == 'C']) / len(s) * 100

def dna_to_prot(valami):
    codon_dict = {'UUU': 'F', 'CUU': 'L', 'AUU': 'I', 'GUU': 'V', 'UUC': 'F', 'CUC': 'L', 'AUC': 'I', 'GUC': 'V', 'UUA': 'L', 'CUA': 'L', 'AUA': 'I', 'GUA': 'V', 'UUG': 'L', 'CUG': 'L', 'AUG': 'M', 'GUG': 'V', 'UCU': 'S', 'CCU': 'P', 'ACU': 'T', 'GCU': 'A', 'UCC': 'S', 'CCC': 'P', 'ACC': 'T', 'GCC': 'A', 'UCA': 'S', 'CCA': 'P', 'ACA': 'T', 'GCA': 'A', 'UCG': 'S', 'CCG': 'P', 'ACG': 'T', 'GCG': 'A', 'UAU': 'Y', 'CAU': 'H', 'AAU': 'N', 'GAU': 'D', 'UAC': 'Y', 'CAC': 'H', 'AAC': 'N', 'GAC': 'D', 'UAA': 'Stop', 'CAA': 'Q', 'AAA': 'K', 'GAA': 'E', 'UAG': 'Stop', 'CAG': 'Q', 'AAG': 'K', 'GAG': 'E', 'UGU': 'C', 'CGU': 'R', 'AGU': 'S', 'GGU': 'G', 'UGC': 'C', 'CGC': 'R', 'AGC': 'S', 'GGC': 'G', 'UGA': 'Stop', 'CGA': 'R', 'AGA': 'R', 'GGA': 'G', 'UGG': 'W', 'CGG': 'R', 'AGG': 'R', 'GGG': 'G'}
    prot_string = ''
    for i in range(0,len(valami)-3,3):
        act = valami[i:i+3]
        if act in codon_dict.keys():
            prot_string += codon_dict[valami[i:i+3]]
        else:
            prot_string += " "
    return prot_string

def rna_to_dna(s):
    return s.replace("U", "T")

def dna_to_rna(s):
    return s.replace("T", "U")

def rev_comp(dna):
    rev = []
    for i in range(len(dna)-1, -1,-1):
        if dna[i] == 'A':
            rev.append('T')
        elif dna[i] == 'T':
            rev.append('A')
        elif dna[i] == 'G':
            rev.append('C')
        elif dna[i] == 'C':
            rev.append('G')
    return ''.join(rev)

def hamming(s,t):
    d = 0
    for i in range(len(s)):
        if s[i] != t[i]:
            d += 1
    return d

def rand_from_gc(s,x):
    val = 1
    for i in range(len(s)):
        if s[i] == "G" or s[i] == "C":
            val *= (x / 2)
        else:
            val *= (1 - x) / 2
    return val

def subs(s,t):
    indices = []
    for i in range(len(s)-len(t)):
        if s[i:i+len(t)] == t:
            indices.append(i)
    return indices

def overlap(s,t):
    return max([i for i in range(len(s)) if t[:i] == s[-i:]])

def pcov(l):
    sorrend = [l[0]]
    l.remove(l[0])
    while len(l) != 0:
        for i in l:
            if sorrend[-1][1:] == i[:-1]:
                sorrend.append(i)
                l.remove(i)

    k = len(sorrend[0])
    n = len(sorrend)

    v = sorrend[0]
    for i in range(n - k):
        v += sorrend[i + 1][-1]
    return v