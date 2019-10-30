fetch_proteins:
	mkdir -p protein_cache
	cat proteins | xargs -I {} wget -O protein_cache/{}.fasta 'https://www.ncbi.nlm.nih.gov/search/api/sequence/{}/?report=fasta'

test:
	python -m unittest alignment_test