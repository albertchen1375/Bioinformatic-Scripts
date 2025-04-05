#!/usr/bin/env python
import sys
import random
import numpy as np

#turns our reads fasta file into a list of reads
def parse_bounds_file(bounds_fn):
    try:
        with open(bounds_fn, 'r') as rFile:
            in_bound = False
            all_reads = []
            current_read = ""
            for line in rFile:
                line = line.strip()
                if line.startswith('>'):
                    if current_read:
                        all_reads.append(current_read)
                        current_read = ""
                    in_bound = True
                    continue
                if in_bound:
                    current_read += line
            # Append the last read after the loop
            if current_read:
                all_reads.append(current_read)
            return all_reads
    except IOError:
        print("Could not read file: ", bounds_fn)
        return None
    

symbol_to_number_dict = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

def GetScore(motif, PWM):
    score = 1
    for i in range(len(motif)):
        rowtomult = symbol_to_number_dict[motif[i]]
        score *= PWM[rowtomult][i]
    return score


def initialize_motif(sequences, motif_length):
    motif_instances = []
    for sequence in sequences:
        start_index = random.randint(0, len(sequence) - motif_length)
        motif_instance = sequence[start_index:start_index + motif_length]
        motif_instances.append(motif_instance)
    return motif_instances

def dividebysum(matrix):
    matrixdiv = matrix/(matrix.sum(axis=0, keepdims=True))
    return matrixdiv



#1 add up all instances and add AAA CCC,....
#2 build initial matrix
#3. gibbs , start of iterations
    #4. choose a sequence to delete
    #5 recompute the pwm
    #6. look at sequence two look at all possible motifs and find probabilities
    #7. from probabilities, get the motif to add and recompute the matrix


def Gibbs_Sampling(sequences, motif_length):

    motif_instances = initialize_motif(sequences, motif_length)
    
    #create the PWM with our motif_instances
    motif_matrix = np.zeros((4, motif_length))
    motif_array = np.array([list(instance) for instance in motif_instances])
    motif_matrix[0, :] = np.sum(motif_array == 'A', axis=0)
    motif_matrix[1, :] = np.sum(motif_array == 'C', axis=0)
    motif_matrix[2, :] = np.sum(motif_array == 'G', axis=0)
    motif_matrix[3, :] = np.sum(motif_array == 'T', axis=0)
    
    
    motif_matrix += 1  # Add pseudocounts


    iterationsper = [0] * len(sequences)
    
    for i in range(len(sequences) * 300):
        if i == 10000 or i == 100000 or i == 500000 or i == 1000000 or i == 1200000:
            print(str(i) + " done")
        
        #select sequence and its associated sequence to remove
        sequence_index_to_remove = random.randint(0, len(sequences) - 1)

        iterationsper[sequence_index_to_remove] += 1
        
        motif_to_remove = motif_instances[sequence_index_to_remove]

        for i in range(len(motif_to_remove)):
            motif_matrix[symbol_to_number_dict[motif_to_remove[i]]][i] -=1


        #compute likelihood
        motif_to_add = ""
        probabilities = [0] * (len(sequences[sequence_index_to_remove]) - motif_length)

        for i in range(len(sequences[sequence_index_to_remove]) - motif_length):
            subsequence = sequences[sequence_index_to_remove][i : i + motif_length]
            subscore = GetScore(subsequence, dividebysum(motif_matrix))
            probabilities[i] = subscore


        #choose the motif to add based on probabilities
        total_score = sum(probabilities)
        probabilities = [score / total_score for score in probabilities]
        motif_to_add_index = np.random.choice(np.arange(len(sequences[sequence_index_to_remove]) -  motif_length),p = probabilities)
        motif_to_add = sequences[sequence_index_to_remove][motif_to_add_index : motif_to_add_index + motif_length]

        motif_instances[sequence_index_to_remove] = motif_to_add
        
        #update motif matrix 
        for i in range(len(motif_to_add)):
            motif_matrix[symbol_to_number_dict[motif_to_add[i]]][i] += 1

        
    motif_matrix /= (motif_matrix.sum(axis=0, keepdims=True))
    
    return motif_matrix


def reverse_complement(pattern: str) -> str:
    """Calculate the reverse complement of a DNA pattern."""
    str = ""
    comp_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    for i in reversed(pattern) : 
         str += comp_dict[i]
            
    return str
    
def main():
   
   boundscentered = parse_bounds_file(r'boundcentered.fasta')
   boundsrandomoffset = parse_bounds_file(r'boundrandomoffset.fasta')

   motif_length = 10
   PWM = Gibbs_Sampling(boundscentered, motif_length)
   
  
   motif_instances = []
   for sequence in boundscentered:
       max_score = 0
       index_with_maxscore = 0
       for i in range(16, 187):
           possible_motif = sequence[i:i+motif_length]
           tempscore = GetScore(possible_motif, PWM)
           if tempscore > max_score:
               max_score = tempscore
               index_with_maxscore = i
       motif_instances.append(sequence[index_with_maxscore: index_with_maxscore + motif_length])
       
   motif_matrix = np.zeros((4, motif_length))

   

# Convert to NumPy array
   motif_array = np.array([list(instance) for instance in motif_instances])

   motif_matrix[0, :] = np.sum(motif_array == 'A', axis=0)
   motif_matrix[1, :] = np.sum(motif_array == 'C', axis=0)
   motif_matrix[2, :] = np.sum(motif_array == 'G', axis=0)
   motif_matrix[3, :] = np.sum(motif_array == 'T', axis=0)

   motif_matrix /= motif_matrix.sum(axis=0, keepdims=True)

   print("done creating pwm")

   
   
   #step 2: going through each boundrandomoffset and locating peaks
   index_sequence = []
   for sequence in boundsrandomoffset:
       max_score = 0
       index_with_maxscore = 0
       for i in range(16, 187):
           possible_motif = sequence[i:i+motif_length]
           tempscore = GetScore(possible_motif, motif_matrix)
           if tempscore > max_score:
               max_score = tempscore
               index_with_maxscore = int(i + (motif_length/2))


       reversecomplement = reverse_complement(sequence)       
       for i in range(16, 187):
           possible_motif = reversecomplement[i:i+motif_length]
           tempscore = GetScore(possible_motif, PWM)
           if tempscore > max_score:
               max_score = tempscore
               index_with_maxscore = int(201 - i - (motif_length/2))

       index_sequence.append(index_with_maxscore)

   file = open('predictions.csv', 'w')
   for index_list, peakindex in enumerate (index_sequence):
        file.write('seq' + str(index_list+1) + '\t' + str(peakindex) + '\n') 
   file.close()

if __name__ == "__main__":
    main()
