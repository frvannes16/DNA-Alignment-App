FROM python:3-alpine

COPY dna_alignment_manager /opt/app/dna_alignment_manager

WORKDIR /opt/app/dna_alignment_manager

RUN apk add --update --no-cache yarn make build-base gcc postgresql-dev gfortran python python-dev py-pip && \
    pip install --no-cache -r requirements.txt

RUN yarn install && yarn run build

EXPOSE 8000:8000

CMD python manage.py runserver 0.0.0.0:8000
