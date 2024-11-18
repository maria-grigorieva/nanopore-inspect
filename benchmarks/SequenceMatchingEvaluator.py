# Part 1: Imports and setup
import time
import csv
from Bio import SeqIO
from Bio import pairwise2
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import subprocess
import os
from collections import defaultdict
import regex
from rapidfuzz.distance import Levenshtein
import numpy as np
from itertools import combinations
import tempfile
from Bio.pairwise2 import format_alignment
from rapidfuzz import fuzz
from scipy.signal import find_peaks


class SequenceMatchingEvaluator:
    def __init__(self, fastq_file, service_sequence, error_rate=0.1):
        self.fastq_file = fastq_file
        self.service_sequence = service_sequence
        self.error_rate = error_rate
        self.max_distance = int(len(service_sequence) * error_rate)
        self.results = []
        self.window_size = len(service_sequence)

    def load_sequences(self):
        """Load sequences from FASTQ file"""
        return list(SeqIO.parse(self.fastq_file, "fastq"))

    def sliding_window(self, sequence, window_size):
        """Generate sliding windows from a sequence"""
        return [sequence[i:i + window_size]
                for i in range(len(sequence) - window_size + 1)]

    # Part 2: Implementation of individual algorithms
    # def levenshtein_search(self, sequences):
    #     """Levenshtein distance using rapidfuzz"""
    #     start_time = time.time()
    #     positions = []
    #
    #     for seq in sequences:
    #         seq_str = str(seq.seq)
    #         windows = self.sliding_window(seq_str, self.window_size)
    #
    #         for i, window in enumerate(windows):
    #             distance = Levenshtein.distance(window, self.service_sequence)
    #             if distance <= self.max_distance:
    #                 positions.append(i)
    #
    #     duration = time.time() - start_time
    #     return {
    #         "algorithm": "Levenshtein (rapidfuzz)",
    #         "duration": duration,
    #         "positions": positions
    #     }

    def levenshtein_search(self, sequences):
        positions = []
        peaks = []
        start_time = time.time()
        for seq in sequences:
            s = str(seq.seq)
            length = len(self.service_sequence)
            substrings = [s[i:i + length] for i in range(len(s) - length + 1)]
            arr = [(s.find(i), fuzz.ratio(i, self.service_sequence)) for i in substrings if
                   fuzz.ratio(i, self.service_sequence) >= 0.9 * 100]
            positions.extend([i for i, j in arr])

            if len(arr) > 0:
                idx, values = zip(*arr)
                peak_indices, _ = find_peaks(values, distance=length)
                if len(peak_indices) > 0:
                    peaks.extend([idx[item] for item in peak_indices.tolist()])
        duration = time.time() - start_time

        return {
            "algorithm": "Levenshtein",
            "duration": duration,
            "positions": peaks
        }
    def smith_waterman(self, sequences):
        """Smith-Waterman local alignment"""
        start_time = time.time()
        positions = []
        for seq in sequences:
            # Perform Smith-Waterman alignment
            alignments = pairwise2.align.localms(str(seq.seq), self.service_sequence, match=2, mismatch=-1, open=-2, extend=-0.5)

            if len(alignments) > 0:
                for a in alignments:
                    matches = sum(1 for a, b in zip(str(seq.seq[a.start:a.end]), a.seqB.strip('-')) if a == b)
                    alignment_length = a.end - a.start
                    accuracy = matches / alignment_length
                    if accuracy >= (1 - self.error_rate):
                        positions.append(a.start)
            # best_alignment = alignments[0]
            # aligned_seq1, aligned_seq2, score, start, end = best_alignment
            # aligned_seq2 = aligned_seq2.strip('-')
            #
            # # Calculate accuracy
            # matches = sum(1 for a, b in zip(str(seq.seq[start:end]), aligned_seq2) if a == b)
            # alignment_length = len(str(seq.seq[start:end]))
            # accuracy = matches / alignment_length
            #
            # # # print(format_alignment(*best_alignment))
            # if len(alignments) > 0 and accuracy >= (1 - self.error_rate):
            #     positions.append(start)
            #     # print(format_alignment(*best_alignment))
            #     # Extract and print the position of occurrence in sequence A
            #     # print(f"Service sequence B found in A at position: {start}")
            #     # print(f"Aligned segment: {A[start:end]}")
        duration = time.time() - start_time
        return {
            "algorithm": "Smith-Waterman",
            "duration": duration,
            "positions": positions
        }
    # def smith_waterman(self, sequences):
    #     """Smith-Waterman local alignment"""
    #     start_time = time.time()
    #     positions = []
    #
    #     for seq in sequences:
    #         alignments = pairwise2.align.localms(
    #             str(seq.seq),
    #             self.service_sequence,
    #             2,  # match score
    #             -1,  # mismatch penalty
    #             -2,  # gap opening penalty
    #             -0.5  # gap extension penalty
    #         )
    #
    #         if alignments:
    #             for alignment in alignments:
    #                 if alignment[2] >= len(self.service_sequence) * (1 - self.error_rate):
    #                     positions.append(alignment[3])  # start position
    #
    #     duration = time.time() - start_time
    #     return {
    #         "algorithm": "Smith-Waterman",
    #         "duration": duration,
    #         "positions": positions
    #     }

    # def bitap_search(self, sequences):
    #     """Bitap algorithm implementation"""
    #     start_time = time.time()
    #     positions = []
    #
    #     pattern_length = len(self.service_sequence)
    #     pattern_mask = {}
    #
    #     # Preprocess pattern
    #     for i in range(pattern_length):
    #         pattern_mask[self.service_sequence[i]] = ~(1 << i)
    #
    #     for seq in sequences:
    #         seq_str = str(seq.seq)
    #         R = ~0
    #         for i, c in enumerate(seq_str):
    #             # Update state
    #             if c in pattern_mask:
    #                 R = (R << 1) | pattern_mask[c]
    #             else:
    #                 R = (R << 1) | ~0
    #
    #             # Check for match
    #             if R & (1 << (pattern_length - 1)) == 0:
    #                 positions.append(i - pattern_length + 1)
    #
    #     duration = time.time() - start_time
    #     return {
    #         "algorithm": "Bitap",
    #         "duration": duration,
    #         "positions": positions
    #     }

    def bwa_search(self, sequences):
        """BWA alignment"""
        start_time = time.time()
        positions = []

        # Create temporary files
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as ref_file, \
                tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as query_file:

            # Write reference sequence
            ref_file.write(f">reference\n{self.service_sequence}\n")

            # Write query sequences
            for i, seq in enumerate(sequences):
                query_file.write(f">seq_{i}\n{str(seq.seq)}\n")

        try:
            # Index reference
            subprocess.run(['bwa', 'index', ref_file.name], check=True)

            # Perform alignment
            result = subprocess.run(
                ['bwa', 'mem', ref_file.name, query_file.name],
                capture_output=True,
                text=True,
                check=True
            )

            # Parse BWA output
            for line in result.stdout.split('\n'):
                if not line.startswith('@') and line.strip():
                    fields = line.split('\t')
                    if len(fields) >= 4:
                        positions.append(int(fields[3]))

        finally:
            # Cleanup
            os.unlink(ref_file.name)
            os.unlink(query_file.name)

        duration = time.time() - start_time
        return {
            "algorithm": "BWA",
            "duration": duration,
            "positions": positions
        }

    def blast_search(self, sequences):
        """BLAST search using local BLAST"""
        start_time = time.time()
        positions = []

        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as db_file, \
                tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as query_file:

            # Write database sequence
            for i, seq in enumerate(sequences):
                db_file.write(f">seq_{i}\n{str(seq.seq)}\n")

            # Write query sequence
            query_file.write(f">query\n{self.service_sequence}\n")

        try:
            # Make BLAST database
            subprocess.run([
                '/usr/local/Cellar/blast/2.16.0/bin/makeblastdb',
                '-in', db_file.name,
                '-dbtype', 'nucl'
            ], check=True)

            # Run BLAST
            start_time = time.time()
            result = subprocess.run([
                '/usr/local/Cellar/blast/2.16.0/bin/blastn',
                '-task', 'blastn-short',  # Optimized for short sequences
                '-query', query_file.name,
                '-db', db_file.name,
                '-strand', 'plus',
                '-perc_identity', '90',
                #'-word_size', str(len(self.service_sequence)-10),  # Reduce word size to increase sensitivity
                '-evalue', '1',  # Allow higher e-value for partial matches
                '-reward', '1',
                '-penalty', '-2',
                '-gapopen', '2',
                '-gapextend', '1',
                '-outfmt', '6',
                '-max_target_seqs', '57654'
            ], capture_output=True, text=True, check=True)
            #
            # # Parse BLAST output
            positions = [
                int(fields[8])
                for line in result.stdout.splitlines()
                if len((fields := line.split('\t'))) >= 9
            ]
            duration = time.time() - start_time
            #
            # for line in result.stdout.split('\n'):
            #     if line.strip():
            #         fields = line.split('\t')
            #         if len(fields) >= 9:
            #             positions.append(int(fields[8]))
        finally:
            # Cleanup
            os.unlink(db_file.name)
            os.unlink(query_file.name)

        return {
            "algorithm": "BLAST",
            "duration": duration,
            "positions": positions
        }

    def jaccard_similarity_search(self, sequences, k=5):
        """Jaccard similarity search"""
        start_time = time.time()
        positions = []

        def calculate_jaccard(set1, set2):
            intersection = len(set1 & set2)
            union = len(set1 | set2)
            return intersection / union if union > 0 else 0

        # Create k-mer set for service sequence
        service_kmers = set()
        for i in range(len(self.service_sequence) - k + 1):
            service_kmers.add(self.service_sequence[i:i + k])

        for seq in sequences:
            seq_str = str(seq.seq)
            for i in range(len(seq_str) - len(self.service_sequence) + 1):
                window = seq_str[i:i + len(self.service_sequence)]
                window_kmers = set()

                for j in range(len(window) - k + 1):
                    window_kmers.add(window[j:j + k])

                similarity = calculate_jaccard(service_kmers, window_kmers)
                if similarity >= (1 - self.error_rate):
                    positions.append(i)

        duration = time.time() - start_time
        return {
            "algorithm": "Jaccard Similarity",
            "duration": duration,
            "positions": positions
        }

    def fuzzy_regex_search(self, sequences):
        """Fuzzy regular expression search"""
        start_time = time.time()
        positions = []

        pattern = f"({self.service_sequence}){{e<={self.max_distance}}}"

        for seq in sequences:
            seq_str = str(seq.seq)
            matches = regex.finditer(pattern, seq_str, regex.BESTMATCH)
            positions.extend(match.start() for match in matches)

        duration = time.time() - start_time
        return {
            "algorithm": "Fuzzy Regex",
            "duration": duration,
            "positions": positions
        }

    # def fuzzy_ahocorasick_search(self, sequences):
    #     """Fuzzy Aho-Corasick search"""
    #     start_time = time.time()
    #     positions = []
    #
    #     # Create Aho-Corasick automaton
    #     automaton = ahocorasick.Automaton()
    #
    #     # Add the pattern to the automaton
    #     automaton.add_word(self.service_sequence, (0, self.service_sequence))
    #
    #     # Make automaton immutable (required before searching)
    #     automaton.make_automaton()
    #
    #     # Process each sequence
    #     for seq in sequences:
    #         seq_str = str(seq.seq)
    #
    #         # Sliding window approach for fuzzy matching
    #         window_size = len(self.service_sequence) + self.max_distance
    #
    #         # Store found positions to avoid duplicates
    #         found_positions = set()
    #
    #         # Iterate through the sequence with overlapping windows
    #         for i in range(len(seq_str) - len(self.service_sequence) + 1):
    #             window = seq_str[i:i + window_size]
    #
    #             # Find exact matches in the current window
    #             for end_index, (pattern_index, pattern) in automaton.iter(window):
    #                 start_index = end_index - len(pattern) + 1
    #
    #                 # Calculate Levenshtein distance for the potential match
    #                 candidate = window[start_index:end_index + 1]
    #                 edit_dist = distance(self.service_sequence, candidate)
    #
    #                 # If within allowed distance and not already found
    #                 if edit_dist <= self.max_distance and i + start_index not in found_positions:
    #                     positions.append(i + start_index)
    #                     found_positions.add(i + start_index)
    #
    #     duration = time.time() - start_time
    #     return {
    #         "algorithm": "Fuzzy Aho-Corasick",
    #         "duration": duration,
    #         "positions": sorted(positions)
    #     }

    # Part 4: Evaluation and results handling
    def evaluate_all(self):
        """Evaluate all implemented algorithms"""
        sequences = self.load_sequences()

        algorithms = [
            self.levenshtein_search,
            #self.smith_waterman,
            # self.bwa_search,
            #self.blast_search,
            # self.ukkonen_search,
            # self.kmer_based_search,
            # self.jaccard_similarity_search,
            # self.fuzzy_regex_search
        ]

        for algo in algorithms:
            try:
                print(f"Running {algo.__name__}...")
                result = algo(sequences)
                self.results.append(result)
                print(f"Completed {algo.__name__}")
            except Exception as e:
                print(f"Error in {algo.__name__}: {str(e)}")

    def save_results(self, output_file):
        """Save results to CSV file"""
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([
                'Algorithm',
                'Duration (seconds)',
                'Matches found'
            ])

            for result in self.results:
                writer.writerow([
                    result['algorithm'],
                    f"{result['duration']:.4f}",
                    len(result['positions'])
                ])

    def save_detailed_results(self, output_file):
        """Save detailed results including positions"""
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([
                'Algorithm',
                'Duration (seconds)',
                'Matches found',
                'Positions'
            ])

            for result in self.results:
                writer.writerow([
                    result['algorithm'],
                    f"{result['duration']:.4f}",
                    len(result['positions']),
                    ','.join(map(str, result['positions']))
                ])


def main():
    # Example usage
    # fastq_file = "/Users/maria/МГУ/Nanopores/DataSamples/FL_16S_22032024_VR/barcode03/merged_FL_16S_22032024_VR.fastq"
    fastq_file = "bigfile.fastq"
    # fastq_file = "merged_PAU55281_barcode12.fastq"
    service_sequence = """
    CCCACTGAGCCGAGACAAGATTCTGCTGTAGTCAGTGCTGCCTGGGAATCTATTTTCACAAAGTTCTCCA
AAAAATGTGATGATCAAAACTAGGAATTAGTGTTCTGTGTCTTAGGCCCTAAAATCTTCCTGTGAATTCC
ATTTTTAAGGTAGTCGAGGTGAACCGCGTCTGGTCTGCAGAGGATAGAAAAAAGGCCCTCTGATACCTCA
AGTTAGTTTCACCTTTAAAGAAGGTCGGAAGTAAAGACGCAAAGCCTTTCCCGGACGTGCGGAAGGGCAA
CGTCCTTCCTCATGGCCGGAAATGGAACTTTAATTTCCCGTTCCCCCCAACCAGCCCGCCCGAGAGAGTG
ACTCTCACGAGAGCCGCGAGAGTCAGCTTGGCCAATCCGTGCGGTCGGCGGCCGCTCCCTTTATAAGCCG
ACTCGCCCGGCAGCGCACCGGGTTGCGGAGGGTGGGCCTGGGAGGGGTGGTGGCCATTTTTTGTCTAACC
CTAACTGAGAAGGGCGTAGGCGCCGTGCTTTTGCTCCCCGCGCGCTGTTTTTCTCGCTGACTTTCAGCGG
GCGGAAAAGCCTCGGCCTGCCGCCTTCCACCGTTCATTCTAGAGCAAACAAAAAATGTCAGCTGCTGGCC
CGTTCGCCCCTCCCGGGGACCTGCGGCGGGTCGCCTGCCCAGCCCCCGAACCCCGCCTGGAGGCCGCGGT
CGGCCCGGGGCTTCTCCGGAGGCACCCACTGCCACCGCGAAGAGTTGGGCTCTGTCAGCCGCGGGTCTCT
CGGGGGCGAGGGCGAGGTTCAGGCCTTTCAGGCCGCAGGAAGAGGAACGGAGCGAGTCCCCGCGCGCGGC
GCGATTCCCTGAGCTGTGGGACGTGCACCCAGGACTCGGCTCACACATGCAGTTCGCTTTCCTGTTGGTG
GGGGGAACGCCGATCGTGCGCATCCGTCACCCCTCGCCGGCAATGGGGGCTTGTGAACCCCCAAACCTGA
CTGACTGGGCCAGTGTGCTGCAAATTGGCAGGAGACGTGAAGGCACCTCCAAAGTCGGCCAAAATGAATG
GGCAGTGAGCCGGGGTTGCCTGGAGCCGTTCCTGCGTGGGTTCTCCCGTCTTCCGCTTTTTGTTGCCTTT
TATGGTTGTATTACAACTTAGTTCCTGCTCTGCAGATTTTGTTGAGGTTTTTGCTTCTCCCAAGGTAGAT
CTCGACCAGTCCCCTCAACGGGGTGTGGGAGAACAGTCATTTTTTTTTGAGAGATCATTTAACATTTAAT
GAATATTTAATTAGAAGATCTAAATGAACATTGGAAATTGTGTTCCTTTAATGGTCATCGGTTTATGCCA
GAGGTTAGAAGTTTCTTTTTTGAAAAATTAGACCTTGGCGATGACCTTGAGCAGTAGGATATAACCCCCA
CAAGCTTAGCGTTCCAATAACGGAACACTAGGCATAATGAAAGACGGAAAAGAAATTTATTCCACACCCC
CCACCCCCAACCCTCCCAGCCGGCAGTCTCCCACAAGAATTGGCTCTGATTTCCTTTAAGAAAATAAGCT
AATGGCCCACTAGGTTGAAAGAAACAAACAGGCCATTCTGCCCACTTATCACGACAAGGTAATTCCGTCC
CAAATCCATACCTTCAATTCCTTAGGATCATCTGGGGGTAGTTGCCAAGAGAGAGAGCAGCCAGGGGTGT
GTATATAAAGGCCCACCTTTAGGGTGATAGCCTGAATCCTGATGATTGAAAGTCAGAAGTCAAGACCATG
TCAGTAAATTTATGAACATGTATATTTGTAGCTTTTTAACCTATTACCTAAGTAGGTCCCCTGGAATTGT
CTAGCAGATACATTTCTTAGCACTATTAGAATTAAGAAATGTAAAAAAACCTCTAGAGTCCACTGTTGCT
GGTAATGTTCTCTAAATAAGGAAAAATATTTTTCCTATCAGCATGGTTTTGTGGAAAAGTAAGGAATGAT
TTTGCCAAGAACTTGTCTAGAGTCTTGAGGAAGTAGCTCAGCTTCAGTAAGCCTCAGTTTACTCAAGGTG
ATCTAAGATCGCTTTTTCTTCTCTTTCTTTTGAGACGGAGTCCCACTCTGTCACCCAGAGCTGGAGTGCA
GTGTCCCCATCTTGGCTCACTGCAACCTCTGCCTCCTCGGTTTCTAGCGATTTTCTCTCAGCCTCCCAAG
TAGCTGGGTTTACAGGCACACACCACCATACCCGGCTGATTTTTTTGTATTTTCAGTAAAGTTGGGCAGG
CTGGCCTCGAACTCCTGACCTCAGGTGATCCGCCCGCTTCTGCCTCCCAAAGTGCTGGGATTACAGGCGT
GAGCCACCGTGCGGGGACTAAGATCCCTTTTGCAAGGGCGGGGACTCCTTGAACTACAGACAAGTGCATG
TCTTTTGTTCTTACTCCATCTAGTGGTACTTTTGCTTCTCCACATTTACAACATTTGTGGTGGTGCAGGG
CCGTGAAGCAAGAATGGCCACAGACTGGGAATGTTTTTCTTTCATTTTTAAAAAGAAAACAGTGATTTAA
TACATTCAATTTAAAAATAAGTCCATAGTTGTTTTCTAGCTCATATTTGAAGTCGTCCCTCTGACATGTC
TATTACAAACAGTCCAAAGGACAATTCTTAGTATTATAAACATTAAATATATTGTTTCTTTTTGGTAAAA
TGATTAAGTTATATAATTCTAGAGTAAAGACGGCCAGTCATCCTGGATGACCATACAGGTCATCCAGTGA
GATAATTTATGTGAAAACTTTCTGTACATAATAATTATTGCAAAAGTGACTTCTTACTGCCACCGATTTG
AGTGACTCTAGCAGATGAAGAAACTGCTGCGGAGAGGGGTTTTAGTTACTTGCCTAAGTTATTTACTTGC
CTAAGTTAATTTTTTTGCAGAGTTGTGACCAGAAGCTAGGCTTCCCAAGTTCTAGTCCAGTGTTCTGTTC
ACTGTATTGCACCCTAAGTCTTTGTTGGTTAACCTGGGGCAGATGGTTAGCATTAACTCCATCGTCACTT
TCTGGTTTAAAATAAAATAACCTCCTTAATGGTTGTAACCTTGGTTTCCACCCCACATCCACTCCCGCCC
TAATACTCTTCGATGTGATTTTGTGCTGTGTGTGTTTCACTTTTTTACAAGTACAGTTGATTCTAGTTAT
TTGCAGTAGTTACATTCTGTGAAGCCTCTGTGAACACTGGATTAGTGAATACTGAACCATTGCTCCCAGA
GGAATTACGGGATTAGGTTCCTGGGAACATTCGGTCACAACATTTTTGTCAAAGAATCGATACATAATCT
TGTGTGTGTTTCTGTTAAAGACATCCTAAAGTATGTTGTTGATTCATTAGCACTGAACTCATATTCAACA
ACACTATCACTTAAGCCTGAATGAAGGCTAACACAAGTATTTTTTTTCTGTAAGACACATCACAGCTTTC
TTGTGCTTAGGAACGCTAGACAGCACTTAAGCACTATGATTGGGGGTCATTTTAAACAGCACAGTCACCA
ACTAAAAGCATACAAATGTGAAAAACTTGCCATTAAATAGACCAGAAATAGGACACTTGTTTATAGTGTC
CTTTTTGGACTACAAAGGCAAAAGAGAACCTATGTGTTGGGGAAATCAAACTTTTCAGGGCTCTGCACAT
GCCCATGAATGACTATGAAAGCACTGCAGGTATTGCTTTTGAGATTATAAACAATTTCGGTGACTAGGTG
AATCCACAAATACAGAATCTATGAATATTGAGGGTTGACTGTACTTCAAGGTCTATGTGACAGCTGTATC
TGAGAGAGACTACAGATGTCCACATGGTTCCTGAGGTACTTAGGGCTTGTTTCATGCCTTTTATATTATT
TAAACCAGTGATTCTCAAAGGGGTATAGTGAGAGAAAGGGTATATGAAGTTTGCAAATGACGATACTTGT
TCTCTTTATTTTGTTACTTAGGTTGTATTTCACATGAAGGATGTGAGAATTTTATGCGTATGAAGGAAGC
ATTAGCTTTAAAAATACTGAGAAACATGCTTTTTCTTTTGTCCTTCAGAATGAGGGGCTTCCTTTGTAAG
GTCTGGAGTTGGTTATGGGTTTTGACTCTGAGGATGACATTGAGTTTTAATTCAGCTTGTCCTGGTTCTG
CAGTCAAAAATTCCATGTTAGCCTGGCTCGAGACTTAAGAGACAAATAGATGGCACAGGCACAAGCGTGT
GTACACATGCACCCACCGTCATACAAAGTGTAAGGGAGGAGTGGTTTTAAATAAGATGAGAATGAAGTGA
CAAAAGGTAACTGTCATTCCTTCATGTTCATTCTTTTAAGGTAATAATATGTGATGGAAGTAAAAACAGC
GAAAATACCCAAGAATTAGCCCAGGATATGCAAATTTAAATGTCTGCAGTGGACAGAAAGGAAACAAAGA
AGTGAATCACACTGGGAATAAGAAAATGGAGACTATAGAGATGGTAGGGAGAATGAGCTTCCCTTTTTAA
GGAGGCAGCAATCTTACAGAAATTGAGCCTAAAATTGTCAGTATTGCAAGAGAAGACAAAAATTCAGATT
TTTATGTCAACACCACTGATTATAGAATGTTGGTAGCTAAATTGAAATAAAATTTAAATATTATATAAGT
GCATGGTAAAGAATTTGGATGGCTTTTGTTCCTGAACCCTGGGAAGTAGCCTCTAAACCCTTGGAATTCA
CCAAGAATGGTGGGAGTGTCTTTTCTTGGTAGGCCTCTTGGACTAGACCCAATAGTTTACACTAACGAGA
TGACTCATGACTGGGGGTTGGTCAGCCATGTGATTAGAGTGTTGGGGCTTTGGGCCCCATGATACTAGCC
CGACCTCTTGAGAGACTTGGAGCCACATGTCTGAGCCTTATGACTAGTGATTCATCCAATCATGCCTGTA
ATGAAGCCCCAATAAAAACTCTGGACACACTGAAGCTTGGATGAACGTCCGTGGTTGGCAGTACTCAGTG
TGCTAGGAGGAAAATGTGTCCCTGAGGACATGAGGAGAGGGAC
    """
    service_sequence = service_sequence.replace("\n", "").replace(" ", "")
    print(service_sequence)
    #service_sequence = "TCCGATTCTGCTTCTTTCTACCTG"
    # service_sequence = "GAGTCTTGTGTCCCAGTTACCAGG"

    print("Starting algorithm evaluation...")

    evaluator = SequenceMatchingEvaluator(fastq_file, service_sequence)
    evaluator.evaluate_all()

    # Save results
    # evaluator.save_results("sequence_matching_comparison50000.csv")
    # evaluator.save_detailed_results("sequence_matching_detailed_comparison50000.csv")
    # evaluator.save_results("long_sequence.csv")
    evaluator.save_detailed_results("long_sequence_detailed.csv")
    print("Evaluation complete. Results saved to CSV files.")


if __name__ == "__main__":
    main()