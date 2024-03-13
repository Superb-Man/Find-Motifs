import random
import operator
import numpy as np

def getRandomMotifs(dna, k, t):
    motifs = []
    for i in range(t):
        j = random.randrange(len(dna[0])-k+1)
        motifs.append(dna[i][j:j+k])
    return motifs

def Score(motifs):
    score = 0
    k = len(motifs[0])
    t = len(motifs)
    for i in range(k):
        count = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        for motif in motifs:
            count['A'] += 1 if motif[i] == 'A' else 0
            count['C'] += 1 if motif[i] == 'C' else 0
            count['G'] += 1 if motif[i] == 'G' else 0
            count['T'] += 1 if motif[i] == 'T' else 0
        score += t - max(count.values())
    return score

def Entropy(motifs):
    entropy = 0.0
    k = len(motifs[0])
    t = len(motifs)
    for i in range(k):
        count = {'A': 1, 'C': 1, 'G': 1, 'T': 1}
        for motif in motifs:
            count[motif[i]] += 1
        l = t + 4 # Laplace's rule of succession
        for key in count:
            count[key] /= l
        entropy += -sum([p * np.log2(p) for p in count.values()])

    return entropy


def profileMatrix(motifs):
    k = len(motifs[0]) # pattern length
    t = len(motifs) #total dna sequences

    profile = {'A': [], 'C': [], 'G': [], 'T': []}
    for i in range(k):
        for nuc in profile.keys():
            profile[nuc].append(1)
    for motif in motifs:
        for j in range(k):
            profile[motif[j]][j] += 1
    for nuc in profile.keys():
        for i in range(k):
            profile[nuc][i] /= t + 4
    return profile

def probableKmer(text, k, profile):
    max_prob = -1
    kmer = text[0:k]
    l = len(text)
    for i in range(l - k + 1):
        pattern = text[i:i+k]
        prob = 1
        for j in range(k):
            prob *= profile[pattern[j]][j]
        if prob > max_prob:
            max_prob = prob
            kmer = pattern
    return kmer


# number of DNA reads
# randomly select k-mers Motifs = (Motif_1, , Motif_t) in each string from Dna

def gibbsSampler(dna, k, N):
    """
    Gibbs sampler for motif discovery
    Args:
    dna: list of strings
    k: int, length of motif
    N: int, number of iterations
    Returns:
    best_motifs: list of strings, best motifs found
    """
    t = len(dna)
    motifs = getRandomMotifs(dna, k, t)
    bestMotifs = motifs
    for j in range(N):
        i = random.randrange(t)
        motifs.pop(i)
        profile = profileMatrix(motifs)
        motifs.insert(i, probableKmer(dna[i], k, profile))
        if Entropy(motifs) < Entropy(bestMotifs):
            bestMotifs = motifs
    return bestMotifs

def MotifsFind(Profile, Dna):
    k = len(Profile['A'])
    t = len(Dna)
    motifs = []
    for i in range(t):
        motif = probableKmer(Dna[i], k, Profile)
        motifs.append(motif)
    return motifs

def RandomizedMotifSearch(dna, k, N):
    t = len(dna)
    # randomly select k-mers Motifs = (Motif1, ..., Motift) in each string from Dna
    Motifs = getRandomMotifs(dna, k, t)
    BestMotifs = Motifs
    
    for j in range(N):
        Profile = profileMatrix(Motifs)
        Motifs = MotifsFind(Profile, dna)
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
        else:
            return BestMotifs

def multipleSeedSearch(dna, numSeeds, k, N,func):
    bestMotifs = []
    for i in range(numSeeds):
        motifs = func(dna, k, N)
        if len(bestMotifs) == 0 or Entropy(motifs) < Entropy(bestMotifs):
            bestMotifs = motifs
        # print(bestMotifs,Entropy(bestMotifs))
    return bestMotifs


if __name__ == '__main__':
   
    dna = []
    # read the input file
    with open('hm03.txt', 'r') as f:
        dna = f.read().splitlines()
    
    
    # motif= gibbsSampler(dna, 8, 1000)
    motif = multipleSeedSearch(dna, 100, 8, 1000,gibbsSampler)
    

    print(motif, Entropy(motif))
    print(Score(motif))
    motif = multipleSeedSearch(dna, 100, 8, 1000,RandomizedMotifSearch)
    print(motif, Entropy(motif))
    print(Score(motif))

    # RSAT = ["CGGGCCCC","TGGGCCGC","CGGGCCGC","CCGGCCGC","CGGGCCGC","CGGGCCCC",
    #         "CGGGCCGC","CGGGCCGC","TGGGCCCC","CGGGCCCC"]
    # # print ENTROPY AND SCORE OF THE RSAT
    # print("RSAT ENTROPY AND SCORE:")
    # print(Entropy(RSAT),Score(RSAT))

    # output = ['AAAATAAA', 'AAAATAAA', 'AAAATAAA', 'AAAAAAAA', 'AAAAGCAA', 'AAAATAAA', 'AAAATCAA', 'AAAATAAA', 'AAAATAAA', 'AAAATCAA'] 
    # # print ENTROPY AND SCORE OF THE OUTPUT
    # print(Entropy(output),Score(output))

    # print(len(dna[0]))
# ['AAAATAAA', 'AAAATAAA', 'AAAATAAA', 'AAAAAAAA', 'AAAAGCAA', 'AAAATAAA', 'AAAATCAA', 'AAAATAAA', 'AAAATAAA', 'AAAATCAA'] 9.540852816243516
    
		
	
