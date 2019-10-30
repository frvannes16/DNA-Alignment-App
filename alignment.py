from Bio import SeqIO
from sys import argv
from pathlib import Path


def get_sequence(args):
    if len(args) == 0:
        return None

    if '/' in args[0]:
        sequence_file = Path(args[0])
    else:
        sequence_file = Path("protein_cache/" + args[0] + ".fasta")
    
    if sequence_file.is_file():
        return sequence_file


def main():
    sequence_file = get_sequence(argv[1:])
    if sequence_file is None:
        print("Could not locate file " + sequence_file)

    for seq_record in SeqIO.parse(sequence_file, "fasta"):
        print(seq_record.id)
        print(seq_record.description)
        print(repr(seq_record.seq))
        print(len(seq_record))

if __name__ == '__main__':
    main()
