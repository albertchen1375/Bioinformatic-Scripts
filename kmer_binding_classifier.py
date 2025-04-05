#!/usr/bin/env python
import sys
import random
import numpy as np

from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, classification_report

import itertools



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
    
def generate_kmer_table(sequence, k):
    characters = ['A', 'T', 'C', 'G']
    sequences = [''.join(seq) for seq in itertools.product(characters, repeat=k)]
    sequence_dict = {sequence: 0 for sequence in sequences}
    for i in range(len(sequence) - k + 1):
        sequence_dict[sequence[i: i+k]] +=1

    return sequence_dict

def feature_extraction(table):
    feature_vector = [count for count in table.values()]
    
    return feature_vector

    
def main():

# 1. Read in Data
    
    bounds = parse_bounds_file(r'bound.fasta')
    notbound = parse_bounds_file(r'notbound.fasta')
    test = parse_bounds_file(r'test.fasta')


# 2. Extract data for training
    features = []

    kmer_length = 6

    for i in range(len(bounds)):
        freqtable = generate_kmer_table(bounds[i], kmer_length)
        features.append(feature_extraction(freqtable))

    for i in range(len(notbound)):
        freqtable = generate_kmer_table(notbound[i], kmer_length)
        features.append(feature_extraction(freqtable))

    
    print("done1")

# 3. Create Labels
    
    bounds_labels = [1] * len(bounds)  # Bound sequences are labeled as 1
    notbound_labels = [0] * len(notbound)  # Not bound sequences are labeled as 0
    labels = bounds_labels + notbound_labels

    X_train, X_test, y_train, y_test = train_test_split(features, labels, test_size=0.2, random_state=38)

# 4. Model Training
    logistic_regression = LogisticRegression()
    logistic_regression.fit(X_train, y_train)

    features = None

    print("done2")
# 5. Model Evaluation
    y_prob = logistic_regression.predict_proba(X_test)
    y_pred = (y_prob[:, 1] > 0.5).astype(int)
    accuracy = accuracy_score(y_test, y_pred)
    report = classification_report(y_test, y_pred)
    print("Accuracy:", accuracy)
    print("Classification Report:\n", report)


# 6. Apply to Test Data
    
    probabilities = []
    for i in range(len(test)):
        freqtable = generate_kmer_table(test[i], kmer_length)
        features = feature_extraction(freqtable)
        prob = logistic_regression.predict_proba(np.array([features]))
        probabilities.append(prob[0, 1])
    

    print("donepredicting")
    top_indices = np.argsort(probabilities)[::-1][:6000]
    print(probabilities[48896])
    print(probabilities[32338])

    file = open('predictions.csv', 'w')
    for index in top_indices:
        file.write('seq' +  str(index + 1) + '\n') 
    file.close()

    

if __name__ == "__main__":
    main()