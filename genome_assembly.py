#!/usr/bin/env python
import sys
from typing import List, Dict, Iterable
import random
from copy import deepcopy
from collections import Counter
from collections import defaultdict
import networkx as nx

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

#this breaks down each read into kmers and removes any kmers with frequency less than k 
def generate_spectrum(reads: List[str], k: int):
    
    all_kmers = {}  # Store indices of reads containing each k-mer
    spectrum = []
    spectrum_index = {}
    
    for index, read in enumerate(reads):
        spectrum_index[read] = index
        for i in range(len(read) - k + 1):
            kmer = read[i:i+k]
            if kmer not in all_kmers:
                all_kmers[kmer] = []
            all_kmers[kmer].append((index, i))

    for kmer in all_kmers:
        if len(all_kmers[kmer]) > 3:
            spectrum.append(kmer)
      
        
     
    return spectrum, spectrum_index

def de_bruijn_kmers(k_mers: List[str]) -> Dict[str, List[str]]:
    """Forms the de Bruijn graph of a collection of k-mers."""
    adj_list= {}
    for kmer in k_mers:
        if kmer[:-1] not in adj_list:
            adj_list[kmer[:-1]] = []
        adj_list[kmer[:-1]].append(kmer[1:])
                                  
    return adj_list


def get_node_degrees(adjacency_list):
    node_degrees = {}
    
    # Increment the out-degree for each incoming node by the number of outgoing nodes
    for incoming_node, outgoing_nodes in adjacency_list.items():
        if incoming_node not in node_degrees:
            node_degrees[incoming_node] = [len(outgoing_nodes), 0]  # Update out-degree
        else:
            node_degrees[incoming_node][0] = len(outgoing_nodes)  # Update out-degree
        
        # Increment the in-degree for each outgoing node by 1
        for outgoing_node in outgoing_nodes:
            if outgoing_node not in node_degrees:
                node_degrees[outgoing_node] = [0, 1]  # Update in-degree
            else:
                node_degrees[outgoing_node][1] += 1  # Update in-degree

    return node_degrees


def findstartnodes(Node_dict_degrees):
    start_nodes = []
    for node, degrees in Node_dict_degrees.items():
        if degrees[1] < degrees[0]:
            start_nodes.append(node)
    return start_nodes



def find_path(adjacency_list, start_node):
    # Initialize a set to keep track of visited nodes
    visited = set()
    path = []
    path.append(start_node)

    currnode = start_node
    while currnode in adjacency_list:
        path.append(adjacency_list[currnode][0])
        currnode = adjacency_list[currnode][0]


    return path


def genome_path(path: List[str]) -> str:
    """Forms the genome path formed by a collection of patterns."""
    output = ""
    output+=(path[0])
    for i in range(1, len(path)):
        output+=(path[i][len(path[i]) - 1])
                              
    return output

def hamming_distance(p: str, q: str) -> int:
    """Calculate the Hamming distance between two strings."""
    hamming_distance = 0
    minimum = min(len(q), len(p))
    for i in range(len(p)):
        if p[i] != q[i]:
            hamming_distance += 1
            
    return hamming_distance

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

 
def readmap(reference:str, reads):

    aligned_reads = {}

    ref_minimizerdict = referenceminimizers(reference)

    for indexread, read in enumerate (reads):

        windowsize = 25
        minsize = 10
        max_d = len(read)/10 - 1
        if max_d < 1 or len(read) < 25:
            continue
    
        #for each window
        for i in range(len(read) - windowsize):

            minimizer = read[i : i + minsize]
            read_index = i
            window = read[i : i + windowsize]

        #calculate the minimizer for the window and position in the read
            for j in range(windowsize - minsize):
                if(window[j:j+minsize] < minimizer): 
                    minimizer = window[j:j+minsize]
                    read_index = j + i

            if minimizer in ref_minimizerdict:
                for ref_min_index in ref_minimizerdict[minimizer]:
                    reference_position = ref_min_index - read_index
                    if 0 <= reference_position <= len(reference) - len(read):
                        if hamming_distance(read, reference[reference_position:reference_position + len(read)]) < 4:
                            aligned_reads[indexread] = reference_position
                    

    finalreadlist = []
    sorted_reads = sorted(aligned_reads.items(), key=lambda x: x[1])
    finalreadlist = [key for key, _ in sorted_reads]
    return finalreadlist               
                    
                        
def main():
    readsfile = r'project2b_reads.fasta'

    reads = parse_reads_file(readsfile)
    readhash = {}
    for index, read in enumerate(reads):
        readhash[hash(read)] = index

    spectrum, spectrum_index = generate_spectrum(reads, 25)
    db= de_bruijn_kmers(spectrum)
    nodedegrees = get_node_degrees(db)
    startnodes = findstartnodes(nodedegrees)

    paths = []
    for startnode in startnodes:
        path_kmers = find_path(db, startnode)
        genome = genome_path(path_kmers)
        paths.append(genome)

    path = ""
    paths = [string for string in paths if len(string) >= 50]
    max_length_string = max(paths, key=len)
    path = max_length_string

    paths.remove(max_length_string)
    
    for pathsingle in paths:
        path+=pathsingle

    aligned_reads = readmap(path, reads)
    

    file = open('predictions.txt', 'w')
    for read in aligned_reads:
        file.write('>read_' + str(read) + '\n') 
    file.close()

if __name__ == "__main__":
    main()
