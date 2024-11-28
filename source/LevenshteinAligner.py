from .BioSequenceAligner import BioSequenceAligner
import os
import time
from rapidfuzz import fuzz
from scipy.signal import find_peaks
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import tempfile


class LevenshteinBio(BioSequenceAligner):
    def __init__(self, db_file, query_string, similarity_score=0.9):
        super().__init__(db_file, query_string, similarity_score)
        self.output_fasta = f"{self.db_name}_LEVENSHTEIN.fasta"
        self.similarity_threshold = similarity_score * 100
        self.peak_distance = len(query_string)
        self.line_matches = []  # store tuples of (line_idx, positions, scores)
        self.matches_df = None

    def find_matches_in_sequence(self, sequence, line_idx):
        """
        Find start positions of similar substrings in a sequence.
        Returns positions and scores where similarity exceeds threshold.
        """
        positions = []
        scores = []

        # Find all positions where similarity exceeds threshold
        for i in range(len(sequence) - self.query_string_length + 1):
            substring = sequence[i:i + self.query_string_length]
            score = fuzz.ratio(self.query_string, substring)
            if score >= self.similarity_threshold:
                positions.append(i)
                scores.append(score)

        # Convert to numpy arrays for peak finding
        if positions:
            positions = np.array(positions)
            scores = np.array(scores)

            # Find peaks in the scores array
            peaks, properties = find_peaks(
                scores,
                distance=self.peak_distance
            )

            # Return line index and the positions/scores at peaks
            return line_idx, positions[peaks], scores[peaks]

        return line_idx, np.array([]), np.array([])

    def calculate_alignments(self):
        """
        Process all sequences and find match positions using peak detection.
        """
        start_time = time.time()
        self.read_fastq_to_array()

        matches_data = []

        # Process each sequence
        for line_idx, sequence in enumerate(self.sequences):
            # Find matches in current sequence
            curr_line_idx, positions, scores = self.find_matches_in_sequence(sequence, line_idx)

            if len(positions) > 0:
                # Store matches for this line
                self.line_matches.append((curr_line_idx, positions, scores))

                # Add each match to the results
                for pos, score in zip(positions, scores):
                    match_seq = sequence[pos:pos + self.query_string_length]
                    matches_data.append({
                        'line_idx': curr_line_idx,
                        'position': pos,
                        'score': score,
                        'match': match_seq
                    })

                self.fuzzy_matches += len(positions)
                self.occurrences.extend(positions.tolist())

        self.duration = time.time() - start_time

        # Convert matches to DataFrame for easy analysis
        self.matches_df = pd.DataFrame(matches_data)

        # if save:
        #     self.parse_results()

        # Print summary
        print(f"Analysis completed in {self.duration:.2f} seconds")
        print(f"Found {self.fuzzy_matches} total matches above {self.similarity_threshold}% similarity")
        print(f"Lines with matches: {len(self.line_matches)}/{len(self.sequences)}")

    def parse_results(self):
        """
        Save results to FASTA and positions to CSV with line indices.
        """
        # Save matched sequences to FASTA
        if not self.matches_df.empty:
            with open(os.path.join(self.output_dir, self.output_fasta), 'w') as f:
                for i, row in enumerate(self.matches_df.itertuples(), 1):
                    identifier = f"line{row.line_idx}_pos{row.position}_score{row.score:.1f}"
                    f.write(f">{identifier}\n{row.match}\n")

            # Save positions and scores with line indices to CSV
            positions_data = []
            for line_idx, positions, scores in self.line_matches:
                positions_data.append({
                    'line_idx': line_idx,
                    'positions': ','.join(map(str, positions)),
                    'scores': ','.join(map(lambda x: f'{x:.1f}', scores))
                })

            positions_df = pd.DataFrame(positions_data)
            positions_csv = os.path.join(self.output_dir, f"{self.db_name}_positions.csv")
            positions_df.to_csv(positions_csv, index=False)

            print(f"Results saved to:")
            print(f"FASTA: {self.output_fasta}")
            print(f"Positions CSV: {positions_csv}")
        else:
            print("No matches found to save.")

    def get_line_matches(self, line_idx):
        """
        Get all matches for a specific line index.
        """
        for stored_line_idx, positions, scores in self.line_matches:
            if stored_line_idx == line_idx:
                return positions, scores
        return np.array([]), np.array([])

    def plot_sequence_matches(self, sequence_idx):
        """
        Plot similarity scores and detected peaks for a specific sequence.
        """

        # Find the matches for this sequence from line_matches
        sequence_match = None
        for line_idx, positions, scores in self.line_matches:
            if line_idx == sequence_idx:
                sequence_match = (positions, scores)
                break

        sequence = self.sequences[sequence_idx]

        # Calculate similarity scores for all positions
        similarity_scores = []
        for i in range(len(sequence) - self.query_string_length + 1):
            substring = sequence[i:i + self.query_string_length]
            score = fuzz.ratio(self.query_string, substring)
            similarity_scores.append(score)

        plt.style.use('seaborn-v0_8-paper')
        plt.figure(figsize=(10, 5))
        plt.plot(similarity_scores, label='Similarity Scores', alpha=0.6)

        # Plot detected peaks if matches were found
        if sequence_match is not None:
            positions, scores = sequence_match
            plt.scatter(positions, scores, color='red', marker='x', s=100,
                        label='Detected Peaks')

            # Add position annotations to peaks
            for pos, score in zip(positions, scores):
                plt.annotate(f'pos:{pos}\nscore:{score:.1f}',
                             xy=(pos, score),
                             xytext=(0, 10),
                             textcoords='offset points',
                             ha='center',
                             va='bottom')

        plt.axhline(y=self.similarity_threshold, color='green', linestyle='--',
                    label=f'Threshold ({self.similarity_threshold}%)')

        # plt.title(f'Sequence {sequence_idx} Matches')
        plt.xlabel('Position', fontsize=16)
        plt.ylabel('Similarity Score (%)', fontsize=16)
        # Set ticks for better readability
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.legend(fontsize=14)
        plt.grid(True, alpha=0.3)

        # Set y-axis limits to show full range of scores
        plt.ylim(0, 105)

        plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)

        # Tight layout for better spacing
        plt.tight_layout()
        plt.savefig('academic_plot.png', dpi=300)

        plt.show()


    def get_matches_summary(self):
        """
        Return summary of matches with line indices.
        """
        if self.matches_df is not None and not self.matches_df.empty:
            summary = {
                'total_matches': self.fuzzy_matches,
                'lines_with_matches': len(self.line_matches),
                'total_lines': len(self.sequences),
                'matches_per_line': self.matches_df['line_idx'].value_counts().to_dict(),
                'average_score': self.matches_df['score'].mean(),
                'execution_time': self.duration
            }
        else:
            summary = {
                'total_matches': 0,
                'lines_with_matches': 0,
                'total_lines': len(self.sequences),
                'matches_per_line': {},
                'average_score': 0,
                'execution_time': self.duration
            }

        return summary
