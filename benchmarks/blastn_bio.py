from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import tempfile

import subprocess

class BlastnBio:

    def __init__(self, db_file, query_string, query_file=None):
        self.db_file = db_file
        self.db_name = os.path.basename(db_file)
        self.query_string = query_string
        self.db_fasta_file = ""
        self.fastq_to_fasta()
        self.query_to_fasta()
        self.blast_output_xml = f"{self.db_name}.xml"

    def read_fastq(self):
        return list(SeqIO.parse(self.db_file, "fastq"))
    def fastq_to_fasta(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as db_fasta_file:
            self.db_fasta_file = db_fasta_file.name
            SeqIO.convert(self.db_file, "fastq", self.db_fasta_file, "fasta")
        print(f"Converted {self.db_file} to {self.db_fasta_file}")

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
        subprocess.run([
            "blastn",
            "-query", self.query_fasta_file,
            '-task', 'blastn-short',
            "-db", self.db_name,
            "-out", self.blast_output_xml,
            "-perc_identity", "90",
            "-outfmt", "5"  # Output in XML format
        ])
        print(f"BLAST search completed. Results saved to {self.blast_output_xml}")

    def parse_blast_results(self):
        with open(self.blast_output_xml) as result_handle:
            blast_record = NCBIXML.read(result_handle)

        # Extract and print alignments
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                print("\n****Alignment****")
                print(f"Sequence: {alignment.title}")
                print(f"Length: {alignment.length}")
                print(f"Score: {hsp.score}")
                print(f"Expect: {hsp.expect}")
                print(f"Identities: {hsp.identities}/{hsp.align_length}")
                print(f"Alignment:\n{hsp.query}\n{hsp.match}\n{hsp.sbjct}")

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
                    print(f"Deleted: {file_path}")
            except Exception as e:
                print(f"Error deleting {file_path}: {e}")
    def calculate_alignments(self):
        self.create_blast_db()
        self.run_blast()
        self.parse_blast_results()
        self.unlink_temp_files()
        self.remove_blast_temp_files(self.db_name)

blastnbio = BlastnBio("../static/sessions/FL_16S_22032024_VR_barcode02/ANG916_pass_barcode02_215459a9_53ca0fe7_0.fastq",
                      "TCGATTCCGTTTGTAGTCGTCTGT")
blastnbio.calculate_alignments()