import  numpy as np
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
    print(count)
    entropy = -sum([p * np.log2(p) for p in count])

    return entropy

print(getHighestEntropy(10)) 
