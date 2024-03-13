from algo import *

def main():
    dna = []
    # read the input file
    with open('hm03.txt', 'r') as f:
        dna = f.read().splitlines()
    
    for i in range(8,25):
        print(i)
        motif = multipleSeedSearch(dna, 100, i, 1000,gibbsSampler)
        print(motif, Entropy(motif),Score(motif))
        print("\n")