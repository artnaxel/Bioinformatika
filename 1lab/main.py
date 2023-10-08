import os
import numpy as np
import collections
from Bio import SeqIO
from itertools import product

MIN_CODON_LENGTH = 100

START_CODON = 'ATG'

STOP_CODONS = {'TAA', 'TAG', 'TGA'}

AMINO_ACID_DICTIONARY = {
    'ttt': 'Phe', 'tct': 'Ser', 'tat': 'Tyr', 'tgt': 'Cys',
    'ttc': 'Phe', 'tcc': 'Ser', 'tac': 'Tyr', 'tgc': 'Cys',
    'tta': 'Leu', 'tca': 'Ser', 'taa': '*', 'tga': '*',
    'ttg': 'Leu', 'tcg': 'Ser', 'tag': '*', 'tgg': 'Trp',
    'ctt': 'Leu', 'cct': 'Pro', 'cat': 'His', 'cgt': 'Arg',
    'ctc': 'Leu', 'ccc': 'Pro', 'cac': 'His', 'cgc': 'Arg',
    'cta': 'Leu', 'cca': 'Pro', 'caa': 'Gln', 'cga': 'Arg',
    'ctg': 'Leu', 'ccg': 'Pro', 'cag': 'Gln', 'cgg': 'Arg',
    'att': 'Ile', 'act': 'Thr', 'aat': 'Asn', 'agt': 'Ser',
    'atc': 'Ile', 'acc': 'Thr', 'aac': 'Asn', 'agc': 'Ser',
    'ata': 'Ile', 'aca': 'Thr', 'aaa': 'Lys', 'aga': 'Arg',
    'atg': 'Met', 'acg': 'Thr', 'aag': 'Lys', 'agg': 'Arg',
    'gtt': 'Val', 'gct': 'Ala', 'gat': 'Asp', 'ggt': 'Gly',
    'gtc': 'Val', 'gcc': 'Ala', 'gac': 'Asp', 'ggc': 'Gly',
    'gta': 'Val', 'gca': 'Ala', 'gaa': 'Glu', 'gga': 'Gly',
    'gtg': 'Val', 'gcg': 'Ala', 'gag': 'Glu', 'ggg': 'Gly'
}

SEQUENCES = [
    'Lactococcus',
    'Escherichia',
    'Streptococcus',
    'Cellulophaga',
    'Coronavirus',
    'Adenovirus',
    'Variolavirus',
    'Herpesvirus',
] 

def filter_short_codons(codons, min_length):
    return [codon for codon in codons if len(codon) >= min_length]

def find_codons_in_frame(frame, start_codon, stop_codons, min_length):
    codons = []
    i = 0
    while i < len(frame) - 1:
        if frame[i] == start_codon:
            for j in range(i, len(frame) - 1):
                if frame[j] in stop_codons:
                    codons.append(''.join(map(str, frame[i:j+1])))
                    i = j
                    break
        i += 1
    return filter_short_codons(codons, min_length)

def get_all_frames(sequence):
    return [
        [sequence[i:i + 3] for i in range(start, len(sequence), 3)]
        for start in range(3)
    ] + [
        [sequence.reverse_complement()[i:i + 3] for i in range(start, len(sequence), 3)]
        for start in range(3)
    ]

def find_all_codons(sequence, start_codon, stop_codons, min_length):
    frames = get_all_frames(sequence)
    codons = [find_codons_in_frame(frame, start_codon, stop_codons, min_length) for frame in frames]
    return [item for sublist in codons for item in sublist]

def find_codon_frequency(sequence, amino_acid_dict):
    codon_data = []
    sequence_string = ''.join(sequence).lower()
    sequence_triplets = [sequence_string[i:i + 3] for i in range(0, len(sequence_string), 3)]
    codon_count = collections.Counter(amino_acid_dict.get(c, 'X') for c in sequence_triplets if c in amino_acid_dict)
    codon_number = len(sequence_triplets)
    for codon in codon_count:
        codon_data.append((round(codon_count[codon] / codon_number, 4)))
    return codon_data

def find_dicodon_frequency(sequence, amino_acid_dict):
    sequence_string = ''.join(sequence).lower()
    
    if len(sequence_string) % 6 != 0:
        sequence_string = sequence_string[:-(len(sequence_string) % 6)]

    sequence_dicodons = [sequence_string[i:i + 6] for i in range(0, len(sequence_string), 6)]
    
    dicodon_aa_pairs = [
        (amino_acid_dict.get(sequence_dicodons[i][:3], 'X'), 
         amino_acid_dict.get(sequence_dicodons[i][3:], 'X')) 
        for i in range(len(sequence_dicodons))
        if sequence_dicodons[i][:3] in amino_acid_dict and sequence_dicodons[i][3:] in amino_acid_dict
    ]

    dicodon_count = collections.Counter(dicodon_aa_pairs)
    all_possible_dicodons = [
        (aa1, aa2) for aa1 in amino_acid_dict.values() for aa2 in amino_acid_dict.values()
    ]
    dicodon_frequency_dict = {dicodon: 0 for dicodon in all_possible_dicodons}
    
    dicodon_number = len(dicodon_aa_pairs)
    for dicodon in dicodon_count:
        frequency = round(dicodon_count[dicodon] / dicodon_number, 4)
        dicodon_frequency_dict[dicodon] = frequency
    dicodon_data = np.array(list(dicodon_frequency_dict.values()))
    return dicodon_data

n = 64
m = 4096

def process_files(folder, files, amino_acid_dict, start_codon, stop_codons, min_length):
    codon_frequencies = []
    dicodon_frequencies = []
    codon_distance_matrix = []
    dicodon_distance_matrix = [] 
    
    for file in files:
        filepath = os.path.join(folder, file)
        if not os.path.exists(filepath):
            print(f"File {filepath} does not exist.")
            continue
        
        for record in SeqIO.parse(filepath, "fasta"):
            codons = find_all_codons(record.seq, start_codon, stop_codons, min_length)
            codon_frequency = find_codon_frequency(codons, amino_acid_dict)
            dicodon_frequency = find_dicodon_frequency(codons, amino_acid_dict)
            codon_frequencies.append(np.array(codon_frequency))
            dicodon_frequencies.append(np.array(dicodon_frequency))
    
    for i in range(len(files)):
        codon_distances = []
        dicodon_distances = []
        for j in range(len(files)):
            codon_distances.append(np.sqrt(np.sum(np.power((codon_frequencies[i] - codon_frequencies[j]), 2)) / n))
            dicodon_distances.append(np.sqrt(np.sum(np.power((dicodon_frequencies[i] - dicodon_frequencies[j]), 2)) / m))
        
        codon_distance_matrix.append(codon_distances)
        dicodon_distance_matrix.append(dicodon_distances)
    

    with open("codon_distance_matrix.txt", "w") as file:
        file.write(str(len(SEQUENCES)) + '\n')
        for i, distances in enumerate(codon_distance_matrix):
            file.write(SEQUENCES[i] + ' ' + ' '.join(map(str, distances)) + '\n')
    
    with open("dicodon_distance_matrix.txt", "w") as file:
        file.write(str(len(SEQUENCES)) + '\n')
        for i, distances in enumerate(dicodon_distance_matrix):
            file.write(SEQUENCES[i] + ' ' + ' '.join(map(str, distances)) + '\n')
            
if __name__ == "__main__":
    folder = "data"
    files = [
                "bacterial1.fasta", "bacterial2.fasta", "bacterial3.fasta", "bacterial4.fasta",
                "mamalian1.fasta", "mamalian2.fasta", "mamalian3.fasta", "mamalian4.fasta"
            ]
    process_files(folder, files, AMINO_ACID_DICTIONARY, START_CODON, STOP_CODONS, MIN_CODON_LENGTH)