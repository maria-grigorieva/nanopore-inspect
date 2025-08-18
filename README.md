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
- python 3.9 (miniconda)
- install Redis
- create virtual environment from the requirements.txt 
```
pip install -r requirements.txt
```

## Required libraries: 
```
* biopython==1.84
* celery==5.4.0
* confsmooth==1.0.0
* flask==3.0.3
* Bootstrap-Flask==2.4.0
* flask-wtf==1.2.1
* matplotlib==3.9.2
* numpy==2.0.2
* pandas==2.2.3
* plotly==5.24.1
* rapidfuzz==3.10.0
* redis==5.0.8
* scipy==1.13.1
* statsmodels==0.14.3
* Werkzeug==3.0.4
* whittaker-eilers==0.1.3
* wtforms==3.1.2
* python-dotenv
```
# Running
1) Set up environment variables: ```FLASK_SECRET_KEY```
2) Create **sessions** folder in **static**, and don't forget to change permissions to **static** folder to 777
3) Run: **python app.py**
4) Start redis: 
```brew services start redis (i.e. on Mac OS)```
5) Start app: ```./start.sh```