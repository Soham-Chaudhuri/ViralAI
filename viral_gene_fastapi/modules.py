from collections import Counter
import numpy as np
import pandas as pd
import dill
from sklearn.feature_extraction.text import TfidfVectorizer
with open('./assets/bacteria&virus.pkl', 'rb') as f:
    model_list = dill.load(f)

model = model_list["model"]
tfidf = model_list["tfidf"]
scaler = model_list["scaler"]

def gc_content(sequence):
    g_count = sequence.count('G')
    c_count = sequence.count('C')
    total_bases = len(sequence)
    gc_percentage = (g_count + c_count) / total_bases * 100
    return gc_percentage

def at_content(sequence):
    a_count = sequence.count('A')
    t_count = sequence.count('T')
    total_bases = len(sequence)
    at_percentage = (a_count + t_count) / total_bases * 100
    return at_percentage

def kmer_frequencies(sequence, k=3):
    kmers = [sequence[i:i+k] for i in range(len(sequence) - k + 1)]
    kmer_counts = Counter(kmers)
    return kmer_counts

def molecular_weight(sequence):
    nucleotide_weights = {
        'A': 331.2,
        'T': 322.2,
        'C': 307.2,
        'G': 347.2
    }
    weight = sum(nucleotide_weights[nuc] for nuc in sequence if nuc in nucleotide_weights)
    return weight

def hydrophobicity(sequence):
    hydrophobic_nucleotides = {'A', 'T'}
    hydrophilic_nucleotides = {'C', 'G'}
    hydrophobic_score = sum(1 for nuc in sequence if nuc in hydrophobic_nucleotides)
    hydrophilic_score = sum(1 for nuc in sequence if nuc in hydrophilic_nucleotides)
    return hydrophobic_score / len(sequence), hydrophilic_score / len(sequence)

def net_charge(sequence):
    basic_charge = {'A': 1, 'T': 1}
    acidic_charge = {'C': -1, 'G': -1}
    charge = sum(basic_charge.get(nuc, 0) for nuc in sequence) + \
            sum(acidic_charge.get(nuc, 0) for nuc in sequence)
    return charge

def dinucleotide_frequencies(seq):
    dinucs = [seq[i:i+2] for i in range(len(seq)-1)]
    freq = Counter(dinucs)
    return {dinuc: freq[dinuc] / len(dinucs) for dinuc in freq}

def sequence_entropy(seq):
    freq = Counter(seq)
    probs = [freq[base] / len(seq) for base in freq]
    return -sum(p * np.log2(p) for p in probs)
def Kmers_contribution(seq, size=6):
    return [seq[x:x+size].lower() for x in range(len(seq) - size + 1)]

def feature_extraction_pipeline(seq):
    df = pd.DataFrame(columns=['Sequence', 'GC', 'AT', 'Molecular Wt', 'Hydrophobic Score', 'Hydrophilic Score', 'Sequence Entropy'])
    seq = seq.replace('\n', '')
    words = Kmers_contribution(seq)
    joined_words = ' '.join(words)
    gc = gc_content(seq)
    at = at_content(seq)
    mw = molecular_weight(seq)
    hydrophobic_score, hydrophilic_score = hydrophobicity(seq)
    s_entropy = sequence_entropy(seq)
    df.loc[len(df)] = [joined_words, gc, at, mw, hydrophobic_score, hydrophilic_score, s_entropy]
    sequence_features = tfidf.transform(df['Sequence']).toarray()
    sequence_feature_names = tfidf.get_feature_names_out()
    numeric_features = df.drop(columns=['Sequence']).values
    combined_features = np.hstack((sequence_features, numeric_features))
    x_dense = pd.DataFrame(combined_features, columns=sequence_feature_names.tolist() + df.drop(columns=['Sequence']).columns.tolist())
    return x_dense
