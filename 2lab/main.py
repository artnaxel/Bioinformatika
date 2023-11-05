import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import numpy as np
from Bio.Blast import NCBIWWW, NCBIXML

ENCODINGS_DICTIONARY = {
    "Sanger Phred+33": (33, 0, 40),
    "Solexa Solexa+64": (64, -5, 40),
    "Illumina 1.3+ Phred+64": (64, 0, 40),
    "Illumina 1.5+ Phred+64": (64, 2, 41),
    "Illumina 1.8+ Phred+33": (33, 0, 41)
}

def determine_encoding(fastq_file):
    """Determine the encoding of a FASTQ file."""
    with open(fastq_file, 'r') as f:
        quality_scores = [line.strip() for idx, line in enumerate(f) if idx % 4 == 3]

    all_scores = ''.join(quality_scores)
    min_val = ord(min(all_scores))
    max_val = ord(max(all_scores))
    diff_sums_dict = {}
    for encoding, (offset, min_score, max_score) in ENCODINGS_DICTIONARY.items():
        diff_from_expected_min = abs(min_val - offset - min_score)
        diff_from_expected_max = abs(max_val - offset - max_score)
        diff_sum = diff_from_expected_min + diff_from_expected_max
        diff_sums_dict[encoding] = diff_sum

    return min(diff_sums_dict, key=diff_sums_dict.get)

def cg_percentage_in_sequence(sequence):
    """Calculate the CG percentage in a sequence."""
    cg_count = sequence.count('C') + sequence.count('G')
    return (cg_count / len(sequence)) * 100

def blast_search(sequence):
    """Perform a BLAST search for a given sequence."""
    result_handle = NCBIWWW.qblast("blastn", "nt", sequence,
                    expect=10,
                    word_size=11,
                    nucl_reward=2,
                    hitlist_size=1,
                    nucl_penalty=-3,
                    gapcosts="5 2",
                    entrez_query="Bacteria [Organism]")
    blast_record = NCBIXML.read(result_handle)
    first_alignment = blast_record.alignments[0]
    first_match = first_alignment.title.split('|')[4].split(' ')[1:]
    organism = ' '.join(first_match)
    return organism

def write_results_to_file(results, filename="results1.txt"):
    """Write the analysis results to a file."""
    with open(filename, 'w') as f:
        f.write("Read ID\tMicroorganism\n")
        for result in results:
            f.write(f"{result['read_id']}\t{result['microorganism']}\n")

def analyze_cg_distribution(fastq_file):
    """Analyze the CG distribution in a FASTQ file and perform BLAST searches on sequences from CG peaks."""
    sequences_and_percentages = []
    with open(fastq_file, 'r') as f:
        lines = f.readlines()
        for i in range(0, len(lines), 4):
            read_id = lines[i].strip()[1:]
            sequence = lines[i+1].strip()
            percentage = cg_percentage_in_sequence(sequence)
            sequences_and_percentages.append((read_id, sequence, percentage))

    percentages = [percentage for _, _, percentage in sequences_and_percentages]
    hist_data, bin_edges = np.histogram(percentages, bins=20)
    peaks, _ = find_peaks(hist_data, prominence=100)
    peak_bin_ranges = [(bin_edges[i], bin_edges[i+1]) for i in peaks]
    results = []

    for bin_range in peak_bin_ranges:
        start, end = bin_range
        sequences_in_bin = [(read_id, seq) for read_id, seq, perc in sequences_and_percentages if start <= perc <= end]
        for read_id, seq in sequences_in_bin[:5]:
            organism = blast_search(seq)
            print(f"Read ID: {read_id}\nSequence: {seq}\nOrganism: {organism}\n")
            results.append({'read_id': read_id, 'sequence': seq, 'microorganism': organism})

    write_results_to_file(results)
    plot_cg_distribution(percentages, hist_data, bin_edges, peaks)

def plot_cg_distribution(percentages, hist_data, bin_edges, peaks):
    """Plot the distribution of CG percentages in the sequences."""
    plt.hist(percentages, bins=20, edgecolor='black', alpha=0.7, color='skyblue')
    for peak in peaks:
        plt.plot(bin_edges[peak], hist_data[peak], "ro")
        plt.annotate(f"{bin_edges[peak]:.2f}%", (bin_edges[peak], hist_data[peak]), textcoords="offset points", xytext=(0,10), ha='center')
    plt.title('Distribution of C/G Nucleotides in Reads', fontsize=16)
    plt.xlabel('Proportion of C/G Nucleotides (%)', fontsize=14)
    plt.ylabel('Number of Reads', fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.show()
    print(f"Number of peaks detected: {len(peaks)}")

def main():
    fastq_file = "reads_for_analysis.fastq"
    print(determine_encoding(fastq_file))
    analyze_cg_distribution(fastq_file)

if __name__ == "__main__":
    main()
