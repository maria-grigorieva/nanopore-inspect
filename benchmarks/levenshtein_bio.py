from BioSequenceAligner import BioSequenceAligner
import os
import time
from rapidfuzz import fuzz
from scipy.signal import find_peaks
import numpy as np
import pandas as pd

class LevenshteinBio(BioSequenceAligner):
    def __init__(self, db_file, query_string, similarity_score=0.9):
        super().__init__(db_file, query_string, similarity_score)
        self.output_fasta = f"{self.db_name}_LEVENSHTEIN.fasta"
        self.similarity_threshold = similarity_score * 100
        self.peak_distance = len(query_string)
        self.peak_prominence = 10
        # Store line indices along with positions and scores
        self.line_matches = []  # Will store tuples of (line_idx, positions, scores)
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
                distance=self.peak_distance,
                #prominence=self.peak_prominence
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
                self.occurrences.extend(positions)

        self.duration = time.time() - start_time

        # Convert matches to DataFrame for easy analysis
        self.matches_df = pd.DataFrame(matches_data)

        self.parse_results()

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

    def plot_line_positions(self, line_idx):
        """
        Plot positions of matches in a specific line.
        """
        import matplotlib.pyplot as plt

        positions, scores = self.get_line_matches(line_idx)

        if len(positions) == 0:
            print(f"No matches found in line {line_idx}")
            return

        sequence = self.sequences[line_idx]

        # Create a binary array showing match positions
        match_indicator = np.zeros(len(sequence))
        for pos in positions:
            match_indicator[pos:pos + self.query_string_length] = 1

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 6), height_ratios=[1, 1])

        # Plot match positions
        ax1.plot(match_indicator, drawstyle='steps-post')
        ax1.set_title(f'Match Positions in Line {line_idx}')
        ax1.set_xlabel('Position')
        ax1.set_ylabel('Match Present')
        ax1.grid(True, alpha=0.3)

        # Plot scores
        ax2.bar(positions, scores, alpha=0.6)
        ax2.set_title('Similarity Scores at Match Positions')
        ax2.set_xlabel('Position')
        ax2.set_ylabel('Score')
        ax2.grid(True, alpha=0.3)

        # Add position annotations
        for pos, score in zip(positions, scores):
            ax1.annotate(f'{pos}',
                         xy=(pos, 1.1),
                         xytext=(pos, 1.1),
                         ha='center',
                         rotation=45)
            ax2.annotate(f'{score:.1f}',
                         xy=(pos, score),
                         xytext=(pos, score + 2),
                         ha='center')

        ax1.set_ylim(-0.1, 1.5)
        plt.tight_layout()
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

long_file = "/Users/maria/МГУ/Nanopores/DataSamples/FL_16S_22032024_VR/barcode03/ANG916_barcode03.fastq"
directory = "/Users/maria/МГУ/Nanopores/DataSamples/FL_16S_22032024_VR/barcode03"
file_pattern = "ANG916_pass"
#file_pattern = "ANG916_barcode03"
#query_pattern = "CCTGGTAACTGGGACACAAGACTC"
query_pattern = "GAGTCTTGTGTCCCAGTTACCAGG"

# Get a list of files starting with 'ANG916_pass'
files = [f for f in os.listdir(directory) if f.startswith(file_pattern)]
# Limit to first 10 files
files_to_read = files[:1]
results = []
# Read and process the files
for file in files_to_read:
    file_path = os.path.join(directory, file)
    # Initialize and run alignment
    levenshtein_aligner = LevenshteinBio(file_path, query_pattern, 0.9)
    levenshtein_aligner.calculate_alignments()

    summary = levenshtein_aligner.get_matches_summary()
    print("Matches summary:", summary)

    # Plot matches for a specific line
    levenshtein_aligner.plot_line_positions(0)  # Plot first line

    # Get matches for a specific line
    positions, scores = levenshtein_aligner.get_line_matches(0)
    print(f"Line 0 matches:")
    print(f"Positions: {positions}")
    print(f"Scores: {scores}")

    # Access all matches through DataFrame
    print("\nAll matches:")
    print(levenshtein_aligner.matches_df)

#     results.append({"filename": file,
#                     "query_length": len(blastnbio.query_string),
#                     "similarity_score": blastnbio.similarity_score,
#                     "db_length": blastnbio.get_db_length(),
#                     "exact_matches": exact_matches,
#                     "fuzzy_matches": fuzzy_matches,
#                     "occurrences": blastnbio.occurrences,
#                     "duration": blastnbio.duration,
#                     "output_fasta_raw": blastnbio.output_fasta_raw,
#                     "output_fasta_blast": blastnbio.output_fasta_blast,
#                     "output_xml": blastnbio.blast_output_xml})
#
# pd.DataFrame(results).to_csv(f"{file_pattern}_results.csv")
