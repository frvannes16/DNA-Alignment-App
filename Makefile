fetch_proteins:
	mkdir -p protein_cache
	cat proteins | xargs -I {} wget -O protein_cache/{}.fasta 'https://www.ncbi.nlm.nih.gov/search/api/sequence/{}/?report=fasta'

clean_proteins:
	rm -rf protein_cache

clean: clean_proteins

all: fetch_proteins

	