from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.SubsMat.MatrixInfo import blosum62
from sys import argv
from pathlib import Path
from multiprocessing import Pool
from dataclasses import dataclass
from typing import Type
import time


@dataclass
class AlignmentJobArgs:
    protein_file: Type[Path]
    dna_string: str

def get_proteins_filepaths():
    protein_cache_dir = Path('protein_cache')
    if not protein_cache_dir.is_dir():
        print('protein_cache dir does not exist. Please run `make fetch_proteins`')
    return [fasta_file for fasta_file in protein_cache_dir.iterdir() if fasta_file.is_file() and fasta_file.suffix == '.fasta']


def get_sequence_file(args):
    if len(args) == 0:
        return None

    if '/' in args[0]:
        sequence_file = Path(args[0])
    else:
        sequence_file = Path("protein_cache/" + args[0] + ".fasta")

    if sequence_file.is_file():
        return sequence_file

def get_protein_sequence(filepath):
    return SeqIO.read(filepath, "fasta")

def stat_sequence(sequence_filepath):
    for seq_record in SeqIO.parse(sequence_filepath, "fasta"):
        print(seq_record.id)
        print(seq_record.description)
        print(repr(seq_record.seq))
        print(len(seq_record))

def align(seq_b, seq_a):
    aligner = PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 1
    aligner.mismatch_score = -2
    aligner.gap_score = -3
    return aligner.align(seq_a, seq_b)
    

def score(seq_a, seq_b):
    aligner = PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 1
    aligner.mismatch_score = -2
    aligner.gap_score = -3

    score = aligner.score(seq_a, seq_b)
    return score


def multiprocessing_find_best_protein_match(dna_string: str):
    jobs = [AlignmentJobArgs(protein_file, dna_string) for protein_file in get_proteins_filepaths()]
    with Pool() as p:
        results = p.map(find_best_protein_match, jobs)
    print(results)
    return max(results, key=lambda x: x.score)

def multiprocessing_find_best_protein_match_by_score(dna_string: str):
    jobs = [AlignmentJobArgs(protein_file, dna_string) for protein_file in get_proteins_filepaths()]
    with Pool() as p:
        scores = p.map(find_protein_match_scores, jobs)
    return sorted(scores, key=lambda x: x[0], reverse=True)

def find_protein_match_scores(alignment_data: Type[AlignmentJobArgs]):
    start = time.process_time()
    candidate_sequence = Seq(alignment_data.dna_string)

    print("Aligning sequence to ", alignment_data.protein_file.name)
    protein_seq = get_protein_sequence(alignment_data.protein_file)
    score_result = score(candidate_sequence, protein_seq)
    print("alignment of %s complete." % protein_seq.id)
    end = time.process_time()
    print("Scoring %s took %s" % (protein_seq.id , str(end - start)))
    return (score_result, protein_seq.id)

def find_best_protein_match(alignment_data: Type[AlignmentJobArgs]):
    start = time.process_time()

    candidate_sequence = Seq(alignment_data.dna_string)
    
    print("Aligning sequence to ", alignment_data.protein_file.name)
    protein_seq = get_protein_sequence(alignment_data.protein_file)
    alignments = align(candidate_sequence, protein_seq)
    print("alignment of %s complete." % protein_seq.id)
    alignment = max(alignments, key=lambda alignment: alignment.score)
    print("Best alignment found: ", alignment)
    end = time.process_time()
    print("Completed alignment with " + protein_seq.id + " in " + str(end - start) + "seconds")
    return alignment

def single_threaded_find_best_protein_match(dna_string: str):
    candidate_sequence = Seq(dna_string)

    best_protein, best_alignment = None, None
    best_alignment_score = -10000

    for protein in get_proteins_filepaths():
        protein_seq = get_protein_sequence(protein)
        alignments = align(candidate_sequence, protein_seq)
        for alignment in alignments:
            if alignment.score > best_alignment_score:
                best_protein = protein_seq.id
                best_alignment = alignment
                best_alignment_score = alignment.score
    
    return (best_protein, best_alignment, best_alignment_score)


def main():
    print(get_proteins_filepaths())
    sequence_file = get_sequence_file(argv[1:])
    if len(argv[1:]) == 0:
        print("Sequence file name not specified.")
        return
    elif sequence_file is None:
        print("Could not locate file " + sequence_file)
        return
    # Load and print out details of the sequence
    stat_sequence(sequence_file)

    dna_string = input('Paste the dna string: ')
    dna_seq = Seq(dna_string)
    protein_record = get_protein_sequence(sequence_file)

    idx = protein_record.seq.find(dna_seq)
    print("Found at index ", idx)
    print(protein_record.seq[idx:idx+len(dna_seq)])

if __name__ == '__main__':
    main()
