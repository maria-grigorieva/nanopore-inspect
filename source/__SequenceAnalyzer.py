from Bio import SeqIO
import numpy as np
from tqdm import tqdm
import configparser
import os
from scipy.stats import variation
from LevenshteinAligner import LevenshteinBio
from BlastnAligner import BlastnBio
from DataSmoother import DataSmoother
from PeakAnalyzer import PeakAnalyzer

class SequenceAnalyzer:
    def __init__(self, session_path):
        self.session = session_path
        self.config = self._load_config()
        self.sequences = []
        self.result_data = {}
        self.records = []
        self.n_records = 0
        self.avg_length = 0
        self.smoother = DataSmoother()
        self.peak_analyzer = PeakAnalyzer()

    def _load_config(self):
        config = configparser.ConfigParser()
        config.read(os.path.join(self.session, 'config.ini'))
        return config

    @staticmethod
    def noise_level(data, reference):
        return variation(data['reads'], axis=0)

    def get_peak_occurrences(self, ref, smoothing_type):
        """Wrapper method for peak analysis"""
        target = 'proportion'
        smooth_target = ref['occurrences'][target].values
        ref['occurrences']['smoothed'] = self.smoother.smooth_data(smooth_target, ref['sequence'], smoothing_type)
        peaks, distance = self.peak_analyzer.analyze_peaks(ref['type'],
                                                ref['occurrences']['smoothed'].values,
                                                ref['occurrences']['reads'].values,
                                                ref['occurrences']['proportion'].values,
                                                ref['occurrences']['consensus'].values)
        ref['peaks'] = peaks
        ref['average_peaks_distance'] = distance

    def biorecords_to_array(self, limit):
        """Read sequences from FASTQ file into array."""
        return [str(rec.seq) for rec in self.records][:limit if limit > 0 else self.n_records]


    def analyze(self):
        file_path = self.config['Parameters']['input_file']
        threshold = float(self.config['Parameters']['threshold'])
        limit = int(self.config['Parameters']['limit'])
        smoothing_type = self.config['Parameters']['smoothing']
        similarity_search = self.config['Parameters']['similarity_search']

        service_sequences = []
        for key, value in self.config['Sequences'].items():
            param = {'type': key, 'sequence': value, 'occurences': []}
            service_sequences.append(param)

        self.records = list(SeqIO.parse(file_path, "fastq"))
        self.n_records = len(self.records)

        self.sequence_strings = self.biorecords_to_array(limit)
        self.avg_length = int(np.mean([len(s) for s in self.sequence_strings]))

        for ref in service_sequences:
            if similarity_search == 'lev':
                levenshtein_aligner = LevenshteinBio(file_path, ref['sequence'],
                                                     similarity_score=threshold)
                levenshtein_aligner.calculate_alignments()
                ref['occurrences'] = levenshtein_aligner.calculate_proportions_and_motifs(self.n_records,
                                                                                          self.avg_length)
            elif similarity_search == 'blast':
                blast_aligner = BlastnBio(file_path, ref['sequence'],
                                          similarity_score=threshold)
                blast_aligner.calculate_alignments()
                ref['occurrences'] = blast_aligner.calculate_proportions_and_motifs()

            ref['noise_level'] = DataSmoother.get_noise_level(ref['occurrences']['reads'])
            ref['total_reads'] = int(np.sum(ref['occurrences']['reads']))
            ref['total_proportion'] = float(np.round(np.sum(ref['occurrences']['proportion']), 4))
            self.get_peak_occurrences(ref, smoothing_type)
            ref['value_counts'] = ref['occurrences'].to_dict('records') if ref['occurrences'].shape[0] > 0 else []

        self.result_data['sequences'] = service_sequences
        self.result_data['parameters'] = {
            'n_records': self.n_records,
            'similarity_search': similarity_search,
            'smoothing': smoothing_type,
            'limit': limit,
            'threshold': threshold,
            'file_path': file_path,
            'avg_noise_level': np.round(np.mean([item['noise_level']
                                                 for item in service_sequences
                                                 if 'noise_level' in item]), 4)
        }

        return self.result_data