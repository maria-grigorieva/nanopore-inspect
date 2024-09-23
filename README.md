# NanoporeInspect
A web tool aimed at evaluating the quality of nanopore sequencing results.

# How it works:
- Input your service sequences (primers, barcodes, adapters).
- Define your desired similarity threshold for searching.
- Choose a smoothing algorithm (optional).
- Upload your FASTQ file.
![screenshot](static/img/input_parameters.png)
- Receive detailed results showing the positional distribution of each specified service sequence across the FASTQ file.
![screenshot](static/img/results.png)
NanoporeInspect empowers users to efficiently evaluate the quality of nanopore sequencing results by providing insightful visualizations and detailed positional information on the distribution of target sequences within the sequencing data.
- The results of each session are saved in /static folder on the server
![screenshot](static/img/sessions.png)

# Installation
- python 3.9
- create virtual environment from the requirements.txt
- install Redis

# Running
1) Run app.py
2) Start redis: brew services start redis (i.e. on Mac OS)
2) Start celery worker: celery -A app.celery_app worker --loglevel INFO
3) Start celery beat: celery -A app.celery_app beat --loglevel INFO
4) Configure settings for your email server