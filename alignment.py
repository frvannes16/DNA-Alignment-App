from Bio import SeqIO
from sys import argv
from pathlib import Path


def get_proteins_filepaths():
    protein_cache_dir = Path('protein_cache')
    if not protein_cache_dir.is_dir():
        print('protein_cache dir does not exist. Please run `make fetch_proteins`')
    return [fasta_file for fasta_file in protein_cache_dir.iterdir() if fasta_file.is_file() and fasta_file.suffix == '.fasta']


def get_sequence(args):
    if len(args) == 0:
        return None

    if '/' in args[0]:
        sequence_file = Path(args[0])
    else:
        sequence_file = Path("protein_cache/" + args[0] + ".fasta")

    if sequence_file.is_file():
        return sequence_file

def stat_sequence(sequence_filepath):
    for seq_record in SeqIO.parse(sequence_filepath, "fasta"):
        print(seq_record.id)
        print(seq_record.description)
        print(repr(seq_record.seq))
        print(len(seq_record))

def main():
    print(get_proteins_filepaths())
    sequence_file = get_sequence(argv[1:])
    if len(argv[1:]) == 0:
        print("Sequence file name not specified.")
        return
    elif sequence_file is None:
        print("Could not locate file " + sequence_file)
        return
    # Load and print out details of the sequence
    stat_sequence(sequence_file)

if __name__ == '__main__':
    main()
