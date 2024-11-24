import configparser
from tqdm import tqdm
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np

class SequenceDistribution:
    def __init__(self, config_path):
        self.config = self._load_config(config_path)
        self.searcher = None
        self.sequences = []
        self.n_records = 0
        self.avg_length = 0

    def _load_config(self, config_path):
        config = configparser.ConfigParser()
        config.read(config_path)
        return config

    def _initialize_searcher(self, n_sequences):
        if n_sequences <= 1000:
            self.searcher = LevenshteinSearcher()
        else:
            self.searcher = BlastSearcher()

    def get_all_occurrences(self, reference, type, all_sequences, threshold=0.75):
        occurrences = []
        for s in tqdm(all_sequences, desc=f"Searching for {type}: {reference}", unit="sequence"):
            try:
                matches = self.searcher.find_matches(s, reference, threshold)
                if matches:
                    occurrences.extend(matches)
            except Exception as e:
                print(f"Error processing sequence: {e}")

        if len(occurrences) > 0:
            return self._process_occurrences(occurrences, reference)
        return None

    def _process_occurrences(self, occurrences, reference):
        # Your existing occurrence processing code here
        # (The rest of the processing functions remain the same)
        pass

    def analyze(self):
        file_path = self.config['Parameters']['input_file']
        records = list(SeqIO.parse(file_path, "fastq"))
        self.n_records = len(records)

        limit = int(self.config['Parameters']['limit'])
        self.sequences = [str(rec.seq) for rec in records][:limit if limit > 0 else self.n_records]
        self.avg_length = int(np.mean([len(s) for s in self.sequences]))

        # Initialize appropriate searcher based on sequence count
        self._initialize_searcher(len(self.sequences))

        # Process sequences
        result_data = self._process_sequences()
        return result_data

    def _process_sequences(self):
        # Your existing main processing code here
        # (The rest of the processing logic remains the same)
        pass


# Usage example:
if __name__ == "__main__":
    SESSION = "path/to/session"
    distributor = SequenceDistribution(os.path.join(SESSION, 'config.ini'))
    results = distributor.analyze()