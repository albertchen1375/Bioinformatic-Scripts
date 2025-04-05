# 🧬 Sequence Analysis Toolkit

A comprehensive toolkit for analyzing DNA sequences using computational and probabilistic methods. This collection includes scripts for genome assembly, motif discovery, transcription factor binding site classification, Hidden Markov Model (HMM) training and inference, and minimizer-based sequence alignment.

Each script is built from scratch using core Python and NumPy to provide an educational, interpretable implementation of key algorithms in bioinformatics.

---

## 📜 Script Descriptions
### genome_assembly.py
Purpose:
Implements a simplified genome assembly pipeline using De Bruijn graphs, followed by minimizer-based read mapping.
This script simulates the core mechanics of genome assemblers by breaking reads into k-mers, building an overlap graph, and traversing it to reconstruct genome contigs. The reconstructed sequence is then used as a reference for aligning reads via fast, heuristic-based minimizers

Input: A collection of short sequencing reads

Output: A list of read indices that align to the reconstructed genome

### minimapaligner.py
Purpose:
A minimal, custom implementation of read mapping using minimizers, a strategy commonly used in fast mappers like minimap2.
This script identifies representative k-mers (minimizers) in windows across each read and quickly compares them to those from the reference, significantly reducing the computational cost of alignment.

Input: A reference sequence and a list of reads

Output: Predicted alignment positions of reads within the reference

### kmer_binding_classifier.py
Purpose:
Predicts whether DNA sequences are transcription factor (TF) binding sites using a logistic regression model trained on k-mer frequency features.
This approach represents each sequence as a vector of k-mer counts and trains a classifier to distinguish between bound and unbound sequences. Ideal for tasks like ChIP-seq motif prediction or enhancer detection.

Input: Labeled bound and unbound DNA sequences, plus unlabeled test sequences

Output: Top-ranked test sequences most likely to be bound (e.g., by a transcription factor)

### HMM_gene_predictor.py
Purpose:
Identifies gene-coding regions in a DNA sequence using a Hidden Markov Model (HMM) and the Baum-Welch algorithm.
After training the HMM to differentiate between coding and non-coding regions, the script uses posterior decoding to calculate the probability of each base being in a gene. The top scoring positions are reported.

Input: A symbolic sequence (e.g., using characters like x, y, z, n)

Output: Positions with the highest posterior probability of being in a gene-coding state

### baum_welch_HMM.py
Purpose:
A class-based, modular implementation of the Baum-Welch algorithm for training Hidden Markov Models.
It iteratively refines the transition and emission probabilities to maximize the likelihood of observed sequences. Designed to be reused, imported, or extended for other probabilistic modeling tasks.

Input: A sequence and initial transition/emission probabilities

Output: Trained transition and emission matrices saved to a file

### gibbs_sampling.py
Purpose:
Performs motif discovery using Gibbs Sampling, a probabilistic technique ideal for identifying conserved patterns in unaligned DNA sequences.
This script randomly initializes motifs and iteratively refines a Position Weight Matrix (PWM) that represents the motif consensus, accounting for both the forward and reverse complement strands.

Input: Sequences with centrally located motifs at unknown positions

Output: Predicted binding site positions for each input sequence

## ✍️ Use Cases
Task	Recommended Script
1. Genome assembly from short reads	genome_assembly.py
2. Fast read alignment	minimapaligner.py
3. Predict TF binding from k-mers	kmer_binding_classifier.py
4. Gene prediction using HMMs	HMM_gene_predictor.py
5. Train HMM from sequence data	baum_welch_HMM.py
6. Motif discovery in DNA	gibbs_sampling.py


