from algo import algo

def main():
    dna = []
    # read the input file
    with open('data/hm03.txt', 'r') as f:
        dna = f.read().splitlines()
    
    for i in range(8,25):
        print(i)
        motif = algo.multipleSeedSearch(dna, 100, i, 1000,algo.gibbsSampler)
        print(motif, algo.Entropy(motif),algo.Score(motif))
        print("\n")