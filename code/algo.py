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
            count[motif[i]] += 1
        score += t - max(count.values())
    return score


def getHighestEntropy(t):
    # total neucleotide = 4 {A,C,G,T}
    # make almost equal distribution of neucleotide--->entropy tends to be maximized
    # assume 25% of each neucleotide will appear at least
    # t >=4 for now
    entropy = 0.0
    count = [int(t/4),int(t/4),int(t/4),int(t/4)]
    # count = [5,1,1,1]
    rem = t%4
    for i in range(rem):
        count[i]+=1
    for i in range(4):
        count[i]/= t 
    entropy = -sum([p * np.log2(p) for p in count])

    return entropy


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


#using gain Information

def gain(entropy, maxEntropy, k):
    return ((maxEntropy - entropy)/k)


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
    highestEntropy = getHighestEntropy(t)
    motifs = getRandomMotifs(dna, k, t)
    bestMotifs = motifs
    gainInfo = gain(Entropy(motifs),highestEntropy*k,k)
    for j in range(N-1):
        i = random.randrange(t)
        motifs.pop(i)
        profile = profileMatrix(motifs)
        motifs.insert(i, probableKmer(dna[i], k, profile))
        entropy = Entropy(motifs)
        if gain(entropy,highestEntropy*k,k) > gainInfo :
            bestMotifs = motifs
            gainInfo = gain(entropy,highestEntropy*k,k)
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
    highestEntropy = getHighestEntropy(t)
    # randomly select k-mers Motifs = (Motif1, ..., Motift) in each string from Dna
    Motifs = getRandomMotifs(dna, k, t)
    BestMotifs = Motifs
    gainInfo = gain(Entropy(Motifs),highestEntropy*k,k)
    
    while True:
        Profile = profileMatrix(Motifs)
        Motifs = MotifsFind(Profile, dna)
        entropy = Entropy(Motifs)
        if gain(entropy,highestEntropy*k,k) > gainInfo:
            BestMotifs = Motifs
            gainInfo = gain(entropy,highestEntropy*k,k)
        else:
            return BestMotifs

def multipleSeedSearch(dna, numSeeds, k, N,func):
    bestMotifs = []
    highestEntropy = getHighestEntropy(len(dna))
    gainInfo = 0
    for i in range(numSeeds):
        motifs = func(dna, k, N)
        entropy = Entropy(motifs)
        if len(bestMotifs) == 0 or gain(entropy,highestEntropy*k,k) > gainInfo:
            bestMotifs = motifs
            gainInfo = gain(entropy,highestEntropy*k,k)
        # print(bestMotifs,Entropy(bestMotifs))
    return bestMotifs



def gibbsSamplerByDifScore(dna, k, N):
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
    best_motifs = motifs
    for j in range(N):
        i = random.randrange(t)
        motifs.pop(i)
        profile = profileMatrix(motifs)
        motifs.insert(i, probableKmer(dna[i], k, profile))
        if Score(motifs)/len(motifs[0]) < Score(best_motifs)/len(best_motifs[0]):
            best_motifs = motifs
    return best_motifs


def RandomizedMotifSearchDifScore(dna, k, N):
    t = len(dna)
    highestEntropy = getHighestEntropy(t)
    # randomly select k-mers Motifs = (Motif1, ..., Motift) in each string from Dna
    Motifs = getRandomMotifs(dna, k, t)
    BestMotifs = Motifs
    while True:
        profile = profileMatrix(Motifs)
        Motifs = MotifsFind(profile, dna)
        if Score(Motifs)/len(Motifs[0]) < Score(BestMotifs)/len(BestMotifs[0]):
            BestMotifs = Motifs
        else:
            return BestMotifs


def multipleSeedSearchByDifScore(dna, numSeeds, k, N,func):
    bestMotifs = []
    for i in range(numSeeds):
        motifs = func(dna, k, N)
        if len(bestMotifs) == 0 or Score(motifs)/len(motifs[0]) < Score(bestMotifs)/len(bestMotifs[0]):
            bestMotifs = motifs
    return bestMotifs
		
	
