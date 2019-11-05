fetch_proteins:
	mkdir -p dna_alignment_manager/api/data/protein_cache
	cat proteins | xargs -I {} wget -O dna_alignment_manager/api/data/protein_cache/{}.fasta 'https://www.ncbi.nlm.nih.gov/search/api/sequence/{}/?report=fasta'

clean_proteins:
	rm -rf protein_cache

clean: clean_proteins

docker-run: 
	docker run --rm -d --name dna_docker -p 8000:8000/tcp -e DNA_ALIGNER_SECRET=$DNA_ALIGNER_SECRET dna-alignment-app:latest

docker-stop:
	docker stop dna_docker

docker-build:
	docker build --rm -f "Dockerfile" -t dna-alignment-app:latest .

start-celery-worker:
	cd dna_alignment_manager && celery -A dna_alignment_manager worker -l info

start-redis:
	docker start redis-instance

stop-redis:
	docker stop redis-instance


all: fetch_proteins

	