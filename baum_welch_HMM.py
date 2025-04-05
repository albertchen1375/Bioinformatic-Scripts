# python3
import sys
from math import *
import numpy as np



class BaumWelch:

    
    def __init__(self):

        def parse_file(bounds_fn):
            try:
                with open(bounds_fn, 'r') as rFile:
                    in_bound = False
                    concatenated_reads = ""  # Initialize an empty string to store all reads
                    current_read = ""
                    for line in rFile:
                        line = line.strip()
                        if line.startswith('>'):
                            if current_read:
                                concatenated_reads += current_read  # Concatenate current read
                                current_read = ""
                            in_bound = True
                            continue
                        if in_bound:
                            current_read += line
                    # Append the last read after the loop
                    if current_read:
                        concatenated_reads += current_read  # Concatenate the last read
                    return concatenated_reads
            except IOError:
                print("Could not read file: ", bounds_fn)
                return None
            
        

        emission = parse_file(r'input.fasta')
        initial_transition = np.array([[0.1, 0.9], [0.2, 0.8]])
        initial_emission = np.array([[0.01, 0.10, 0.10, 0.79] , [0.80, 0.12, 0.04, 0.04]])

        States = ['G', 'N']
        Emissionstates = ['n', 'x', 'y', 'z']

        iterations  = 2


        emmision1 = "yzzxzzzyxyxzyxzzyyzzxxzyzyyyyyyxzyxxyzzzyzxyxxxxyxzzzzyzxyxxzyyyyxyzyxzzyzyxyxzzxyxxxzxyxyyxzxxyzxzz"
        States1 = ['A', 'B', 'C', 'D']
        Emissionstates1= ['x', 'y', 'z']

        initial_transition1 = np.array([ [0.056, 0.443,0.334, 0.167], [0.826, 0.052,0.011, 0.111], [0.022, 0.163, 0.811, 0.004], [0.141, 0.696, 0.055, 0.108]])
        initial_emission1 = np.array([[0.34, 0.084, 0.576] , [0.196, 0.592, 0.212], [0.344, 0.609, 0.047], [0.32, 0.487, 0.193]])
    
    
        transition, emission = self.BaumWelchLearning(emission, initial_transition, initial_emission, Emissionstates, States, iterations)
        self.saveTransitionAndEmission(Emissionstates, States, transition, emission)

    

    def softDecode(self, xList, transition, emission):
        n = len(xList)
        l = transition.shape[0]
        forward = [[0 for _ in range(l)] for __ in range(n)]
        backward = [[0 for _ in range(l)] for __ in range(n)]


        epsilon = 1e-10
        for k in range(l):
            forward[0][k] = emission[k, xList[0]]/l
        for i in range(1, n):
            for k in range(l):
                forward[i][k] = sum([forward[i-1][kpre]*transition[kpre, k]*emission[k, xList[i]] for kpre in range(l)])
        fsink = sum(forward[n-1]) + epsilon 

        for k in range(l):
            backward[n-1][k] = 1
        for i in range(n-2, -1, -1):
            for k in range(l):
                backward[i][k] = sum([backward[i+1][kpre]*transition[k, kpre]*emission[kpre, xList[i+1]] for kpre in range(l)])
        
        Pr = np.zeros((l, n), dtype = float)
        for i in range(n):
            for k in range(l):
                Pr[k, i] = forward[i][k]*backward[i][k]/fsink
        
        Pr2 = np.zeros((l, l, n-1), dtype = float)
        for k1 in range(l):
            for k2 in range(l):
                for i in range(n-1):
                    Pr2[k1, k2, i] = forward[i][k1]*transition[k1, k2]*emission[k2, xList[i+1]]*\
                    backward[i+1][k2]/fsink

        return Pr, Pr2
    
    def estimateParameters(self, xList, Pr, Pr2, nAlphabet):
        n = len(xList)
        l = Pr2.shape[0]
        transition = np.zeros((l, l), dtype = float)
        emission = np.zeros((l, nAlphabet), dtype = float)
        for k1 in range(l):
            for k2 in range(l):
                transition[k1, k2] = sum(Pr2[k1, k2, :])
        for k in range(l):
            for i in range(n):
                emission[k, xList[i]] += Pr[k, i]

        for i in range(l):
            sum1 = sum(transition[i,:])
            if 0 == sum1:
                transition[i,:] += 1/l
            else:
                transition[i,:] /= sum1
            sum2 = sum(emission[i,:])
            if 0 == sum2:
                emission[i,:] += 1/nAlphabet
            else:
                emission[i,:] /= sum2
        return transition, emission

    def BaumWelchLearning(self, x, transition, emission, alphabet, states, iterNo):
        x2ind = {alphabet[i]:i for i in range(len(alphabet))}
        xList = [x2ind[x[i]] for i in range(len(x))]
        for _ in range(iterNo):
            Pr, Pr2 = self.softDecode(xList, transition, emission)
            transition, emission = self.estimateParameters(xList, Pr, Pr2, len(alphabet))
        return transition, emission

    def saveTransitionAndEmission(self, alphabet, states, fullTransition, emission):
        f = open('result.txt', 'w')
        print(' '.join([' '] + states))
        f.write('\t'+'\t'.join(states) + '\n')
        for i in range(fullTransition.shape[0]):
            print(' '.join([states[i]] + ['{:.3f}'.format(a) for a in fullTransition[i, :]]))
            f.write('\t'.join([states[i]] + ['{:.3f}'.format(a) for a in fullTransition[i, :]]) + '\n')
        print('--------')
        f.write('--------'+'\n')
        print(' '.join([' '] + alphabet))
        f.write('\t'+'\t'.join(alphabet)+'\n')
        for i in range(emission.shape[0]):
            print(' '.join([states[i]] + ['{:.3f}'.format(a) for a in emission[i, :]]))
            f.write('\t'.join([states[i]] + ['{:.3f}'.format(a) for a in emission[i, :]])+'\n')
        f.close()

if __name__ == '__main__':
    BaumWelch()