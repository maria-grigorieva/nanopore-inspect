from abc import ABC, abstractmethod
from Bio import SeqIO
import os
import tempfile

class BioSequenceAligner(ABC):
    """
    Abstract base class for biological sequence alignment algorithms.
    """

    def __init__(self, db_file, query_string, similarity_score=0.9):
        self.db_file = db_file
        self.db_name = os.path.basename(db_file).split('.')[0]
        self.query_string = query_string
        self.query_string_length = len(query_string)
        self.output_dir = f"{self.db_name}"
        self.duration = 0
        self.fuzzy_matches = 0
        self.exact_matches = 0
        self.similarity_score = similarity_score
        self.sequences = []
        self.occurrences = []
        self.raw_alignments_list = []

        # Create output directory
        self.create_output_dir()

    def create_output_dir(self):
        """Create output directory if it doesn't exist."""
        try:
            os.mkdir(self.output_dir)
            print(f"Directory '{self.output_dir}' created successfully.")
        except FileExistsError:
            print(f"Directory '{self.output_dir}' already exists.")
        except PermissionError:
            print(f"Permission denied: Unable to create '{self.output_dir}'.")
        except Exception as e:
            print(f"An error occurred: {e}")

    def read_fastq_to_array(self):
        """Read sequences from FASTQ file into array."""
        self.sequences = [str(seq.seq) for seq in list(SeqIO.parse(self.db_file, "fastq"))]

    def get_db_length(self):
        """Get the number of sequences in the database."""
        return len(list(SeqIO.parse(self.db_file, "fastq")))

    def find_exact_matches(self):
        """Find exact matches of query string in sequences."""
        import re
        return sum([len(re.findall(self.query_string, s))
                    for s in self.sequences])

    def fastq_to_fasta(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as db_fasta_file:
            self.db_fasta_file = db_fasta_file.name
            SeqIO.convert(self.db_file, "fastq", self.db_fasta_file, "fasta")
        print(f"Converted {self.db_file} to {self.db_fasta_file}")

    @abstractmethod
    def calculate_alignments(self):
        """
        Abstract method to calculate sequence alignments.
        Must be implemented by derived classes.
        """
        pass

    @abstractmethod
    def parse_results(self):
        """
        Abstract method to parse alignment results.
        Must be implemented by derived classes.
        """
        pass