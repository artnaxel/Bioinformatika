from Bio import SeqIO, SearchIO
from Bio.Blast import NCBIWWW
from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO
from Bio.Align import substitution_matrices
from os import path

HOMO_SAPIENS_ALBUMIN_FILE = 'homo_sapiens_albumin.fasta'
HUMAN_ID = "sp|P02768.2|"

def blast_search_similar_sequences(homo_sapiens_albumin):
    if not path.exists("blast.xml"):
        result = NCBIWWW.qblast(
                program="blastp",
                database="swissprot",
                sequence=homo_sapiens_albumin.seq,
                entrez_query='mammals[Organism]')
        with open(".\\blast.xml", "w") as out:
            out.write(result.read())
    blast_result = SearchIO.read("blast.xml", "blast-xml")
    return blast_result


def read_homo_sapiens_albumin_file():
    for record in SeqIO.parse(HOMO_SAPIENS_ALBUMIN_FILE, "fasta"):
        homo_sapiens_albumin = record
    return homo_sapiens_albumin

def load_align(file_path):
    return AlignIO.read(file_path, "fasta")

def write_blast_results_to_fasta(blast_results, sequence_len):
    with open("parsed_blast_sequences.fasta", "w") as out_file:
        for hit in blast_results.hits:
            hsp = hit.hsps[0]
            align_length = hsp.query_end - hsp.query_start
            query_cover = (align_length / sequence_len) * 100
            if "Albumin" in hit.description and query_cover >= 80:
                record = hit.hsps[0].hit
                SeqIO.write(record, out_file, "fasta")

def find_sequences(human_seq, alignment):
    window_size = 15
    most_unique_score = float('inf')
    most_similar_score = -1
    most_unique_sequence = None
    most_similar_sequence = None
    blosum62 = substitution_matrices.load("BLOSUM62")
    for i in range(len(human_seq) - window_size + 1):
        human_window = human_seq[i:i + window_size].replace('-', '')
        if len(human_window) < window_size:
            human_window = human_seq[i:i+window_size + window_size - len(human_window)].replace('-', '')
        current_score = 0

        for record in alignment:
            if record.id == HUMAN_ID:
                continue

            mammal_window = record.seq[i:i+window_size].replace('-', '')
            if len(mammal_window) < window_size:
                mammal_window = record.seq[i:i+window_size + window_size - len(mammal_window)].replace('-', '')
            score = sum(blosum62.get((human_aa, mammal_aa), 0) for human_aa, mammal_aa in zip(human_window, mammal_window))
            current_score += score

        if current_score < most_unique_score:
            most_unique_score = current_score
            most_unique_sequence = human_window

        if current_score > most_similar_score:
            most_similar_score = current_score
            most_similar_sequence = human_window

    return most_unique_sequence, most_unique_score, most_similar_sequence, most_similar_score


def mafft_alignment():
    mafft_exe = "mafft-win\mafft.bat"
    in_file = "parsed_blast_sequences.fasta"
    mafft_cline = MafftCommandline(mafft_exe, input=in_file)
    stdout, stderr = mafft_cline()
    if not path.exists("mafft_L-INS-i.fasta"):
        with open("mafft_L-INS-i.fasta", "w") as output:
            output.write(stdout)

homo_sapiens_albumin = read_homo_sapiens_albumin_file()
blast_result = blast_search_similar_sequences(homo_sapiens_albumin)

write_blast_results_to_fasta(blast_result, len(homo_sapiens_albumin.seq))
mafft_alignment()

alignemnt = load_align("mafft_L-INS-i.fasta")
most_unique_sequence, most_unique_score, most_similar_sequence, most_similar_score = find_sequences(alignemnt[0].seq, alignemnt)

print('Most unique sequence: ', most_unique_sequence, most_unique_score, '\nMost similar sequence: ', most_similar_sequence, most_similar_score)

