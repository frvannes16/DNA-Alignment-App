import unittest
import random
from Bio import SeqIO
from alignment import get_proteins_filepaths, align, score

random.seed('ginkgo')


class TestSequenceAlignment(unittest.TestCase):
    def test_alignment_on_same_protein(self):
        protein_paths = get_proteins_filepaths()
        for f in protein_paths:
            print('aligning: ' + f.name)
            self.check_alignment_of_protein(f)

    def test_alignment_on_different_protein(self):
        protein_paths = get_proteins_filepaths()
        file_distance = random.randint(len(protein_paths))
        for i in range(len(protein_paths)):
            target_protein = protein_paths[i]
            alternate_protein = protein_paths[i-file_distance]
            target_sequence = [x for x in SeqIO.parse(target_protein, "fasta")][0].seq
            alternate_subsequence = self.get_subsequence_of_protein_file(alternate_protein, 500)
            match_score, mismatch_score = score(target_sequence, alternate_subsequence)
            self.assertTrue(match_score < mismatch_score)


    def check_alignment_of_protein(self, protein_path):
        sequence = [x for x in SeqIO.parse(protein_path, "fasta")][0]
        subsequence = self.get_random_sequence_substring(sequence, 500)
        align(sequence, subsequence)
        match_score, mismatch_score = score(sequence, subsequence)
        self.assertTrue(match_score > mismatch_score)

    def get_subsequence_of_protein_file(self, path, size):
        sequence = [x for x in SeqIO.parse(path, "fasta")][0]
        return self.get_random_sequence_substring(sequence, 500)

    def get_random_sequence_substring(self, seq, size):
        seq_len = len(seq)
        start_idx = random.randint(0, seq_len - size)
        return seq[start_idx:start_idx + size]


if __name__ == '__main__':
    unittest.main()
