#!/usr/bin/env python
import sys
from collections import Counter

#turns our reads fasta file into a list of reads
def parse_reads_file(reads_fn):
    try:
        with open(reads_fn, 'r') as rFile:
            first_line = True
        
            count = 0
            all_reads = []
            for line in rFile:
                if first_line:
                    first_line = False
                    continue
                #get rid of all lines that are not the actual read
                if "read" in line:
                    continue
                ends = line.strip().split(',')
                all_reads.append(ends[0])
        return all_reads
    except IOError:
        print("Could not read file: ", reads_fn)
        return None

#turns our reference genome into one big string
def parse_ref_file(ref_fn):
    
    try:
        with open(ref_fn, 'r') as gFile:
            first_line = True
            ref_genome = ''
            for line in gFile:
                if first_line:
                    first_line = False
                    continue
                ref_genome += line.strip()
        return ref_genome
    except IOError:
        print("Could not read file: ", ref_fn)
        return None

#creates minimizers for reference
def referenceminimizers(reference:str):
    minimizer_hash = {}
    windowsize = 25
    minsize = 10
    #for each window
    for i in range(len(reference) - windowsize):
        minimizer = reference[i : i + minsize]
        index = i
        window = reference[i : i + windowsize]
    
    #calculate the minimizer for the window
        for j in range(windowsize - minsize):
            if(window[j:j+minsize] < minimizer): 
                minimizer = window[j:j+minsize]
                index = j + i

        minimizer_key = minimizer

        if minimizer_key not in minimizer_hash:
            minimizer_hash[minimizer_key] = []

        if index not in minimizer_hash[minimizer_key]:
            minimizer_hash[minimizer_key].append(index)

    return minimizer_hash

    
genomes = []
for i in range (100):
    genomes.append(parse_ref_file(r'project1c_genome_' + str(i) + '.fasta'))

readspath = r'project1c_reads.fasta'
    
reads = parse_reads_file(readspath)

#create list of minimizer hashes
genomeminimizers= []
for i in range(len(genomes)):
    genomeminimizers.append(referenceminimizers(genomes[i]))

#create a list containing index
readgenomes= [-1] * len(reads)

for index, read in enumerate(reads):
        windowsize = 25
        minsize = 10
        max_d = len(read)/10 - 1
        if max_d < 1 or len(read) < 25:
            continue  

        #for each read create a list containing how many minimizers in for coverage
        amt_minimizers_genomelist = [0] * len(genomes)

        #calculate the minimizers for each window in the read
        readminimizerdict = {}
        for i in range(len(read) - windowsize):
            minimizer = read[i : i + minsize]
            read_index = i
            window = read[i : i + windowsize]
            for j in range(windowsize - minsize):
                if(window[j:j+minsize] < minimizer): 
                    minimizer = window[j:j+minsize]
                    read_index = j + i

            if minimizer not in readminimizerdict:
                readminimizerdict[minimizer] = read_index

        #check if the minimizer is in each genome minimizer list and increment in amt_minimizer_list if present      
            for z, genomemin in enumerate(genomeminimizers):
                if minimizer in genomemin:
                    amt_minimizers_genomelist[z] += 1

        #whichever genome has most minimizers in the read is the genome, genome_number_ contains indices of max
        max_value = max(amt_minimizers_genomelist)
        Genome_Number_ = [i for i, num in enumerate(amt_minimizers_genomelist) if num == max_value]

        #if there is one that has max, that is the genome
        if len(Genome_Number_ ) == 1:
            readgenomes[index] = Genome_Number_[0]

        #case where multiple genomes have multiple minimizers    
        else:

            #create new dictionary of genomes that have same minimizers
            equalgenomes = {}
            for genomeindex in Genome_Number_:
                equalgenomes[genomeindex] = genomeminimizers[genomeindex]

            for readminimizer in readminimizerdict:
                readindex = readminimizerdict[readminimizer]

                for genomenumbers in equalgenomes:
                        zebra = equalgenomes[genomenumbers]

                        if readminimizer not in genomeminimizers[genomenumbers]:
                            continue

                        for ref_min_index in genomeminimizers[genomenumbers][readminimizer]:
                        
                            reference_position = ref_min_index - read_index
                            mismatch_number = 0

                        #compare read with reference starting at the referece_position for substitutions
                            for k in range(len(read)):
                                if (reference_position + k > len(zebra) - 1):
                                    break
                    
                                if zebra[reference_position + k] != read[k]:
                                    mismatch_number += 1

                                #check for indels
                                if mismatch_number == 1 and k > 2:

                                    #check for insertion
                                    reference_w_ins = zebra[reference_position + k : reference_position + len(read)]
                                    readinsertion = read[k + 1:]
                                    if (len(reference_w_ins) == len(readinsertion)):
                                        for nucleotide in range(len(readinsertion)):
                                        
                                            if reference_w_ins[nucleotide] != readinsertion[nucleotide]:
                                                break
                                            if nucleotide == len(readinsertion) - 1:
                                                readgenomes[index] = genomenumbers
                                            
                                    #check for deletion
                                    reference_w_deletion = zebra[reference_position + k + 1 : reference_position + len(read)]
                                    readdeletion = read[k:]
                                    if (len(reference_w_deletion) == len(readdeletion)):
                                        for nucleotide in range(len(reference_w_deletion)):

                                            if reference_w_deletion[nucleotide] != readdeletion[nucleotide]:
                                                break
                                            if nucleotide == len(reference_w_deletion) - 1:
                                                readgenomes[index]= genomenumbers

                                    if mismatch_number > max_d:
                                        break 

                            readgenomes[index] = genomenumbers
                            break
                        break
                break

            if readgenomes[index] == -1:
                readgenomes[index] = Genome_Number_[0]


file = open('predictions.txt', 'w')

for p in range(len(readgenomes)):
    if not readgenomes[p]:
        genomeprediction = 0
    else:
        genomeprediction = readgenomes[p]
    file.write('>read_' + str(p) + ' Genome_Number_' + str(genomeprediction) + '\n') 
file.close()