from Bio import SeqRecord, SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAlignment
from pathlib import Path
from typing import Type, List
from dataclasses import dataclass
from multiprocessing import Pool
from enum import Enum


class ProteinStore:
    protein_paths = [fasta_file for fasta_file in Path('protein_cache').iterdir(
    ) if fasta_file.is_file() and fasta_file.suffix == '.fasta']

    @classmethod
    def list_all_proteins(cls) -> List[Type[Path]]:
        return cls.protein_paths

    @classmethod
    def get_protein_record(cls, protein_path: type(Path)) -> type(SeqRecord):
        if protein_path in cls.protein_paths:
            return SeqIO.read(protein_path, "fasta")
        return None


class SearchType(Enum):
    STRING_MATCH = 1
    ALIGNMENT = 2


@dataclass
class SearchResult:
    match_found: bool
    protein_rec: type(SeqRecord)
    start_idx: int
    end_idx: int
    search_method: Type[SearchType]
    alignment: type(PairwiseAlignment)


class DnaSearchTool:
    @classmethod
    def search(cls, query: str) -> Type[SearchResult]:
        query_seq = Seq(query)
        # Create search tasks, one for each protein file.
        search_tasks: List[Type[SearchTask]] = [SearchTask(
            protein_filepath, query_seq) for protein_filepath in ProteinStore.list_all_proteins()]
        # Issue the tasks to a pool
        with Pool() as thread_pool:
            results: List[Type[SearchResult]] = thread_pool.map(
                _search_map_func, search_tasks)

        return cls._best_match(results)

    @classmethod
    def _best_match(cls, search_results: List[Type[SearchResult]]) -> SearchResult:
        max_score: int = 0
        best_result: SearchResult = None
        for result in search_results:
            if result.match_found:
                if result.search_method is SearchType.STRING_MATCH:
                    best_result = result
                    break
                elif result.alignment.score > max_score:
                    max_score = result.alignment.score
                    best_result = result
        return best_result


@dataclass
class SearchTask:
    protein_file: Type[Path]
    query_seq: Type[Seq]

    def executeSearch(self) -> Type[SearchResult]:
        # Load the protein sequence.
        target_rec = ProteinStore.get_protein_record(self.protein_file)
        if target_rec is not None:
            # Perform a direct string search for potentially fastest results.
            start_idx: int = target_rec.seq.find(self.query_seq)
            if start_idx > -1:
                # Successful string match!
                return SearchResult(match_found=True,
                                    protein_rec=target_rec,
                                    start_idx=start_idx,
                                    end_idx=start_idx+len(self.query_seq),
                                    search_method=SearchType.STRING_MATCH,
                                    alignment=None)
            # Perform alignment as backup.
            matching_alignment = self.align(target_rec)
            if matching_alignment is not None:
                # Successful alignment!
                return SearchResult(match_found=True,
                                    protein_rec=target_rec,
                                    start_idx=0,
                                    end_idx=0,
                                    search_method=SearchType.ALIGNMENT,
                                    alignment=matching_alignment)
        # No results found
        return SearchResult(match_found=False, search_method=SearchType.ALIGNMENT, protein_rec=None, start_idx=0, end_idx=0, alignment=None)

    def align(self, target_rec: type(SeqRecord)) -> type(PairwiseAlignment):
        print('Pairwise alignment is not yet supported.')
        return None


def _search_map_func(task: Type[SearchTask]) -> Type[SearchResult]:
    return task.executeSearch()


def main():
    input_dna: str = input(
        'Enter your DNA string to find a matching protein: ')
    best_match = DnaSearchTool.search(input_dna)
    if best_match is None:
        print("Match could not be found :(")
    else:
        print("Match found")
        print(best_match)


if __name__ == '__main__':
    main()
