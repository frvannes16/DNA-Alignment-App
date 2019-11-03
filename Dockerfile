FROM python:3-alpine

COPY . /opt/app/

WORKDIR /opt/app

# install dependencies (yarn, maybe sqllite?)
RUN apk add --update --no-cache yarn make build-base gcc gfortran python python-dev py-pip && \
    pip install --no-cache -r requirements.txt

WORKDIR /opt/app/dna_alignment_manager

RUN yarn install && yarn run build

EXPOSE 8000

CMD python manage.py runserver 0.0.0.0:8000
