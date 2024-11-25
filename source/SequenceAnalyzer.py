from dataclasses import dataclass
from typing import List, Dict, Any, Optional
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


@dataclass
class AnalysisParameters:
    """Data class to store analysis parameters"""
    file_path: str
    threshold: float
    limit: int
    smoothing_type: str
    similarity_search: str


@dataclass
class SequenceData:
    """Data class to store sequence information"""
    type: str
    sequence: str
    occurrences: Dict[str, Any]
    noise_level: float = 0.0
    total_reads: int = 0
    total_proportion: float = 0.0
    peaks: List[Dict] = None
    average_peaks_distance: float = 0.0
    value_counts: List[Dict] = None


class SequenceAnalyzer:
    """Main class for analyzing biological sequences"""

    def __init__(self, session_path: str):
        """Initialize the sequence analyzer with session path"""
        self.session = session_path
        self.parameters = self._load_parameters()
        self.sequences: List[SequenceData] = []
        self.result_data: Dict[str, Any] = {}
        self.records = []
        self.n_records: int = 0
        self.avg_length: int = 0
        self.sequence_strings: List[str] = []

        # Initialize analysis components
        self.smoother = DataSmoother()
        self.peak_analyzer = PeakAnalyzer()

    def _load_parameters(self) -> AnalysisParameters:
        """Load and parse configuration parameters"""
        config = configparser.ConfigParser()
        config.read(os.path.join(self.session, 'config.ini'))

        return AnalysisParameters(
            file_path=config['Parameters']['input_file'],
            threshold=float(config['Parameters']['threshold']),
            limit=int(config['Parameters']['limit']),
            smoothing_type=config['Parameters']['smoothing'],
            similarity_search=config['Parameters']['similarity_search']
        )

    def _initialize_sequences(self) -> List[SequenceData]:
        """Initialize sequence data from configuration"""
        config = configparser.ConfigParser()
        config.read(os.path.join(self.session, 'config.ini'))

        return [
            SequenceData(
                type=key,
                sequence=value,
                occurrences={},
                peaks=[],
                value_counts=[]
            )
            for key, value in config['Sequences'].items()
        ]

    def _load_bio_records(self) -> None:
        """Load biological records from file"""
        self.records = list(SeqIO.parse(self.parameters.file_path, "fastq"))
        self.n_records = len(self.records)
        self.sequence_strings = self._biorecords_to_array()
        self.avg_length = int(np.mean([len(s) for s in self.sequence_strings]))

    def _biorecords_to_array(self) -> List[str]:
        """Convert bio records to string array"""
        limit = self.parameters.limit
        return [str(rec.seq) for rec in self.records][:limit if limit > 0 else self.n_records]

    def _process_sequence(self, sequence: SequenceData) -> None:
        """Process individual sequence data"""
        sequence.occurrences = self._get_alignments(sequence)
        sequence.noise_level = DataSmoother.get_noise_level(sequence.occurrences['reads'])
        sequence.total_reads = int(np.sum(sequence.occurrences['reads']))
        sequence.total_proportion = float(np.round(np.sum(sequence.occurrences['proportion']), 4))

        self._analyze_peaks(sequence)
        sequence.value_counts = (sequence.occurrences.to_dict('records')
                                 if sequence.occurrences.shape[0] > 0 else [])

    def _get_alignments(self, sequence: SequenceData) -> Dict[str, Any]:
        """Get sequence alignments based on similarity search method"""
        if self.parameters.similarity_search == 'lev':
            aligner = LevenshteinBio(
                self.parameters.file_path,
                sequence.sequence,
                similarity_score=self.parameters.threshold
            )
        else:  # blast
            aligner = BlastnBio(
                self.parameters.file_path,
                sequence.sequence,
                similarity_score=self.parameters.threshold
            )

        aligner.calculate_alignments()
        return aligner.calculate_proportions_and_motifs(self.n_records, self.avg_length)

    def _analyze_peaks(self, sequence: SequenceData) -> None:
        """Analyze peaks in sequence data"""
        target = 'proportion'
        smooth_target = sequence.occurrences[target].values
        sequence.occurrences['smoothed'] = self.smoother.smooth_data(
            smooth_target,
            sequence.sequence,
            self.parameters.smoothing_type
        )

        peaks, distance = self.peak_analyzer.analyze_peaks(
            sequence.type,
            sequence.occurrences['smoothed'].values,
            sequence.occurrences['reads'].values,
            sequence.occurrences['proportion'].values,
            sequence.occurrences['consensus'].values
        )

        sequence.peaks = peaks
        sequence.average_peaks_distance = distance

    def _prepare_result_data(self) -> Dict[str, Any]:
        """Prepare final result data"""
        return {
            'sequences': self.sequences,
            'parameters': {
                'n_records': self.n_records,
                'similarity_search': self.parameters.similarity_search,
                'smoothing': self.parameters.smoothing_type,
                'limit': self.parameters.limit,
                'threshold': self.parameters.threshold,
                'file_path': self.parameters.file_path,
                'avg_noise_level': np.round(
                    np.mean([seq.noise_level for seq in self.sequences if seq.noise_level]),
                    4
                )
            }
        }

    def analyze(self) -> Dict[str, Any]:
        """Main analysis method"""
        self._load_bio_records()
        self.sequences = self._initialize_sequences()

        for sequence in self.sequences:
            self._process_sequence(sequence)

        self.result_data = self._prepare_result_data()
        return self.result_data