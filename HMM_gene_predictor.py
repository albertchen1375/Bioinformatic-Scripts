import sys
import random
import numpy as np
import mpmath


np.set_printoptions(suppress=True, threshold=np.inf)

#turns our reads fasta file into a list of reads
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


def ForwardBackward(symbols, transitionmatrix, emissionmatrix):
        
        emission_length = len(symbols)
        num_states = transitionmatrix.shape[0]

    #create forward and backward storage
        forward = [[0] * num_states for _ in range(emission_length)]
        backward = [[0] * num_states for _ in range(emission_length)]   


    #forward
        for y in range(num_states):
            forward[0][y] = emissionmatrix[y, symbols[0]]/num_states

        for x in range(1, emission_length):
            for y in range(num_states):
                forward[x][y] = sum([forward[x-1][state]*transitionmatrix[state, y]*emissionmatrix[y, symbols[x]] for state in range(num_states)])
     
            #normalization
            sum_forward_i = sum(forward[x])
            for y in range(num_states):
                forward[x][y] /= sum_forward_i

    #backward
        for y in range(num_states):
            backward[emission_length-1][y] = 1
        for x in range(emission_length-2, -1, -1):
            for y in range(num_states):
                backward[x][y] = sum([backward[x+1][state]*transitionmatrix[y, state]*emissionmatrix[state, symbols[x+1]] for state in range(num_states)])
        
            #normalization
            sum_backward_i = sum(backward[x])
            for y in range(num_states):
                backward[x][y] /= sum_backward_i


        #calculation of state probabilities and transition probabilities
        state_prob = np.zeros((num_states, emission_length))
        for x in range(emission_length):

            for y in range(num_states):
                
                    state_prob[y, x] = forward[x][y]*backward[x][y] 
                    
        
        transition_prob = np.zeros((num_states, num_states, emission_length-1))

        for state1 in range(num_states):
            for state2 in range(num_states):
                for x in range(emission_length-1):
                        transition_prob[state1, state2, x] = forward[x][state1]*transitionmatrix[state1, state2]*emissionmatrix[state2, symbols[x+1]]
                        backward[x+1][state2] 
        
        return state_prob, transition_prob
    
def estimateParameters(symbols, state_prob, joint_transition_prob, num_emission_states):
        
        #create storage and find how many emission states and emitted symbols
        size = joint_transition_prob.shape[0]

        transition = np.zeros((size, size))
        emission = np.zeros((size, num_emission_states))

        n = len(symbols)

        #update parameters
        for v in range(size):
            for j in range(size):
                transition[v, j] = sum(joint_transition_prob[v, j, :])

        for m in range(size):
            for v in range(n):
                emission[m, symbols[v]] += state_prob[m, v]


        #normalization of transition and emission while accounting for 0
        for v in range(size):

            small_number = 1e-10

            emission_sums = sum(emission[v,:])
            if emission_sums == 0:
                emission[v,:] += small_number
            else:
                emission[v,:] /= emission_sums


            state_totals = sum(transition[v,:])
            if state_totals == 0:

                transition[v,:] += small_number

            else:
                transition[v,:] /= state_totals

        return transition, emission

def BaumWelchLearning( x, transition, emission, alphabet,  iterNo):
        
        char_map = {}
        for i in range(len(alphabet)):
            char_map[alphabet[i]] = i

        x_numbers = []
        for char in x:
            x_numbers.append(char_map[char])

        z= len(alphabet)
        for _ in range(iterNo):
            
            post_probab, joint_probab = ForwardBackward(x_numbers, transition, emission)
            transition, emission = estimateParameters(x_numbers, post_probab, joint_probab, z)
            
        return transition, emission
    

def main():

    emission = parse_file(r'input.fasta')
    

    initial_transition = np.array([[0.5, 0.5], [0.5, 0.5]])
    initial_emission = np.array([[0.22, 0.37, 0.19, 0.22] , [0.98, 0.005,0.014, 0.001]])

    States = ['G', 'N']
    Emissionstates = ['n', 'x', 'y', 'z']

    iterations  = 150

    
    transitionmatrix, emissionmatrix = BaumWelchLearning(emission, initial_transition, initial_emission, Emissionstates, iterations)
  
    print(transitionmatrix)
    print(emissionmatrix)
    
    char_map = {Emissionstates[i]:i for i in range(len(Emissionstates))}
    symbol_in_nums = [char_map[emission[i]] for i in range(len(emission))]
    posterior = ForwardBackward(symbol_in_nums, initial_transition, initial_emission)

    posterior = posterior[0]
  
    prob_gene = []
    for i in range(len(posterior[0]) - 1):

        prob_gene.append(posterior[0][i] / (posterior[0][i] + posterior[1][i] ))

    top_positions = np.argsort(prob_gene)[::-1]

    file = open('predictions.csv', 'w')
    for i in range(0, 50000):
        file.write( str(top_positions[i] + 1) + '\n') 
    file.close()

if __name__ == "__main__":
    main()
    