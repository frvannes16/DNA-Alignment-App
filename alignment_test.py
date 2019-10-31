import unittest
import random
from pathlib import Path
from Bio import SeqIO
from alignment import get_proteins_filepaths, multiprocessing_find_best_protein_match_by_score, multiprocessing_find_best_protein_match

random.seed('ginkgo')


class TestSequenceAlignment(unittest.TestCase):

    # This method tests the find_best_protein_match method.
    def test_search_NC_000852(self):
        protein_file = self.get_sequence_file("NC_000852")
        random_subsequence = str(self.get_subsequence_of_protein_file(protein_file, 50).seq)
        result = multiprocessing_find_best_protein_match(random_subsequence)
        print(result)
        best_protein = result
        self.assertTrue("NC_000852" in best_protein.target)

    # This method tests the find_best_protein_match method.
    def test_search_NC_007346(self):
        protein_file = self.get_sequence_file("NC_007346")
        random_subsequence = str(self.get_subsequence_of_protein_file(protein_file, 50).seq)
        result = multiprocessing_find_best_protein_match(random_subsequence)
        print(result)
        best_protein = result
        self.assertTrue("NC_007346" in best_protein.target)

    # This method tests the find_best_protein_match method.
    def test_search_NC_023640(self):
        protein_file = self.get_sequence_file("NC_023640")
        random_subsequence = str(self.get_subsequence_of_protein_file(protein_file, 50).seq)
        result = multiprocessing_find_best_protein_match(random_subsequence)
        print(result)
        best_protein = result
        self.assertTrue("NC_023640" in best_protein.target)
        

    @staticmethod
    def get_subsequence_of_protein_file(path, size):
        sequence = [x for x in SeqIO.parse(path, "fasta")][0]
        return TestSequenceAlignment.get_random_sequence_substring(sequence, 500)

    @staticmethod
    def get_random_sequence_substring(seq, size):
        seq_len = len(seq)
        start_idx = random.randint(0, seq_len - size)
        return seq[start_idx:start_idx + size]

    @staticmethod
    def get_sequence_file(protein_name):
        sequence_file = Path("protein_cache/" + protein_name + ".fasta")

        if sequence_file.is_file():
            return sequence_file


if __name__ == '__main__':
    unittest.main()
