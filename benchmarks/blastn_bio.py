from BioSequenceAligner import BioSequenceAligner
from Bio.Blast import NCBIXML
from Bio import SeqIO
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import tempfile
import time
import subprocess
import pandas as pd

class BlastnBio(BioSequenceAligner):

    def __init__(self, db_file, query_string, similarity_score=0.9):
        super().__init__(db_file, query_string, similarity_score)
        self.db_fasta_file = ""
        self.fastq_to_fasta()
        self.query_to_fasta()
        self.output_dir = f"{self.db_name}"
        self.create_output_dir()
        self.blast_output_xml = f"{self.db_name}.xml"
        self.output_fasta_raw = f"{self.db_name}_RAW.fasta"
        self.output_fasta_blast = f"{self.db_name}_BLAST.fasta"
        self.duration = 0
        self.fuzzy_matches = 0
        self.exact_matches = 0
        self.raw_alignments_list = []
        self.blast_alignments_list = []
        self.occurrences = []
        self.sequences = []

    def query_to_fasta(self):

        record = SeqRecord(
            Seq(self.query_string),
            id="seq1",
            description="Query"
        )
        # Save the record to a FASTA file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as query_fasta_file:
            self.query_fasta_file = query_fasta_file.name
            SeqIO.write(record, self.query_fasta_file, "fasta")
        print(f"Sequence saved to {self.query_fasta_file}")

    def create_blast_db(self):
        subprocess.run([
            "makeblastdb",
            "-in", self.db_fasta_file,
            "-dbtype", "nucl",
            "-out", self.db_name
        ])
        print(f"Created BLAST database '{self.db_name}' from {self.db_fasta_file}")

    def run_blast(self):
        start_time = time.time()
        subprocess.run([
            "blastn",
            "-query", self.query_fasta_file,
            '-task', 'blastn-short',
            "-db", self.db_name,
            "-out", os.path.join(self.output_dir, self.blast_output_xml),
            "-strand", "plus",
            "-word_size", str(int(self.query_string_length * 0.7)), # 70% of the query length
            #"-word_size", str(int(len(self.query_string) - len(self.query_string)*(1-self.similarity_score))),
            "-perc_identity", str(int(self.similarity_score*100)),
            "-evalue", "1",  # Allow higher e-value for partial matches
            "-reward", "1",
            "-penalty", "-2",
            "-gapopen", "2",
            "-gapextend", "1",
            "-max_target_seqs", str(self.get_db_length()),
            "-outfmt", "5"  # Output in XML format
        ])
        self.duration = time.time() - start_time
        print(f"BLAST search completed in {self.duration} seconds. Results saved to {self.blast_output_xml}")

    def parse_results(self):
        with open(os.path.join(self.output_dir, self.blast_output_xml)) as result_handle:
            blast_record = NCBIXML.read(result_handle)

        # Extract and print alignments
        for alignment in blast_record.alignments:
            accession = int(alignment.accession)
            curr_sequence = self.sequences[accession]
            for hsp in alignment.hsps:
                print("\n****Alignment****")
                print(f"Sequence: {alignment.title}")
                print(f"Length: {alignment.length}")
                print(f"Score: {hsp.score}")
                #print(f"Expect: {hsp.expect}")
                #print(f"Identities: {hsp.identities}/{hsp.align_length}")
                print(f"Alignment:\n{hsp.query}\n{hsp.match}\n{hsp.sbjct}")
                start_index = hsp.sbjct_start - hsp.query_start - 1
                end_index = hsp.sbjct_end if hsp.sbjct_end - start_index >= self.query_string_length else \
                            start_index + self.query_string_length + 1
                self.raw_alignments_list.append(curr_sequence[start_index:end_index])
                self.blast_alignments_list.append(hsp.sbjct)
                self.occurrences.append(start_index)

        self.fuzzy_matches = len(blast_record.alignments)
        print(f"Total {len(blast_record.alignments)} alignments have been found out of {self.get_db_length()}.")

    def fasta_from_aligments(self):
        with open(os.path.join(self.output_dir, self.output_fasta_raw), 'w') as f:
            for i, sequence in enumerate(self.raw_alignments_list, start=1):
                identifier = f"seq{i}"
                f.write(f">{identifier}\n{sequence}\n")
        print(f"FASTA file '{self.output_fasta_raw}' generated.")

        with open(os.path.join(self.output_dir, self.output_fasta_blast), 'w') as f:
            for i, sequence in enumerate(self.blast_alignments_list, start=1):
                identifier = f"seq{i}"
                f.write(f">{identifier}\n{sequence}\n")
        print(f"FASTA file '{self.output_fasta_blast}' generated.")

    def unlink_temp_files(self):
        os.unlink(self.db_fasta_file)
        os.unlink(self.query_fasta_file)

    def remove_blast_temp_files(self, db_prefix):
        # List of common BLAST database file extensions
        extensions = ["nhr", "nin", "nsq", "ndb", "not", "ntf", "nto", "njs"]

        # Loop through and delete files with the specified extensions
        for ext in extensions:
            file_path = f"{db_prefix}.{ext}"
            try:
                if os.path.exists(file_path):
                    os.remove(file_path)
                    #print(f"Deleted: {file_path}")
            except Exception as e:
                print(f"Error deleting {file_path}: {e}")

    def calculate_alignments(self):
        self.create_blast_db()
        self.read_fastq_to_array()
        self.run_blast()
        self.parse_results()
        self.fasta_from_aligments()
        self.unlink_temp_files()
        self.remove_blast_temp_files(self.db_name)



long_file = "/Users/maria/МГУ/Nanopores/DataSamples/FL_16S_22032024_VR/barcode03/ANG916_barcode03.fastq"
directory = "/Users/maria/МГУ/Nanopores/DataSamples/FL_16S_22032024_VR/barcode03"
file_pattern = "ANG916_pass"
#file_pattern = "ANG916_barcode03"
#query_pattern = "CCTGGTAACTGGGACACAAGACTC"
query_pattern = "GAGTCTTGTGTCCCAGTTACCAGG"
# blastnbio = BlastnBio("/Users/maria/МГУ/Nanopores/DataSamples/FL_16S_22032024_VR/barcode03/ANG916_barcode03.fastq",
#                       "CCTGGTAACTGGGACACAAGACTC")
# #blastnbio.calculate_alignments()
# print(blastnbio.get_db_length())


# Get a list of files starting with 'ANG916_pass'
files = [f for f in os.listdir(directory) if f.startswith(file_pattern)]
# Limit to first 10 files
files_to_read = files[:1]
results = []
# Read and process the files
for file in files_to_read:
    file_path = os.path.join(directory, file)
    blastnbio = BlastnBio(file_path, query_pattern, similarity_score=0.9)
    blastnbio.calculate_alignments()
    exact_matches = blastnbio.find_exact_matches()
    fuzzy_matches = blastnbio.fuzzy_matches
    results.append({"filename": file,
                    "query_length": len(blastnbio.query_string),
                    "similarity_score": blastnbio.similarity_score,
                    "db_length": blastnbio.get_db_length(),
                    "exact_matches": exact_matches,
                    "fuzzy_matches": fuzzy_matches,
                    "occurrences": blastnbio.occurrences,
                    "duration": blastnbio.duration,
                    "output_fasta_raw": blastnbio.output_fasta_raw,
                    "output_fasta_blast": blastnbio.output_fasta_blast,
                    "output_xml": blastnbio.blast_output_xml})

pd.DataFrame(results).to_csv(f"{file_pattern}_results.csv")
