from scipy.signal import find_peaks, peak_prominences, peak_widths
import numpy as np
class PeakAnalyzer:
    """Class for handling peak detection and analysis operations"""

    def __init__(self):
        pass

    def analyze_peaks(self, reference_type,
                            smoothed_array,
                            reads_array,
                            proportion_array,
                      consensus_array,
                      ):
        """Main method for peak analysis"""
        list_of_peaks = []
        peaks, initial_bases = find_peaks(smoothed_array)
        prominences = peak_prominences(smoothed_array, peaks)[0]
        avg_prominence = np.mean(prominences)
        widths = peak_widths(smoothed_array, peaks, rel_height=1)[0]
        peak_indices, bases = find_peaks(smoothed_array, prominence=avg_prominence)

        if len(peak_indices) == 0:
            peak_indices = peaks
            bases = initial_bases

        if bases:
            extremums = self._process_peaks(reads_array,
                                            proportion_array,
                                            consensus_array,
                                            peak_indices,
                                            bases)
            list_of_peaks = extremums if len(extremums) > 0 else []
            average_peaks_distance = self.calculate_average_peaks_distance(list_of_peaks)

        return list_of_peaks, average_peaks_distance

    def _process_peaks(self, reads_array,
                             proportion_array,
                             consensus_array,
                             peak_indices,
                             bases):
        """Process all peaks and calculate their properties"""
        extremums = []
        for i in range(0, len(peak_indices)):
            peak_index = peak_indices[i]
            extremums.append(self.aggregate_peak_values(i, reads_array,
                             proportion_array, consensus_array, peak_index, bases))

        for i in range(1, len(extremums)):
            extremums[i]['peak_dist'] = extremums[i]['peak_index'] - extremums[i - 1]['peak_index']

        return extremums

    @staticmethod
    def aggregate_peak_values(step, reads_array,
                             proportion_array, consensus_array, peak_index, bases):
        """Calculate aggregate values for a single peak"""
        left_bases = bases['left_bases'][step] if 'left_bases' in bases else 0
        right_bases = bases['right_bases'][step] if 'right_bases' in bases else 0
        total_proportion = np.round(np.sum(proportion_array[left_bases:right_bases]), 4)
        total_occurrences = np.round(np.sum(reads_array[left_bases:right_bases]), 4)
        consensus = consensus_array[int(peak_index)]

        return {
            'peak_index': int(peak_index),
            'left_bases': int(left_bases),
            'right_bases': int(right_bases),
            'total_proportion': float(total_proportion),
            'total_reads': int(total_occurrences),
            'motif_consensus': str(consensus)
        }

    @staticmethod
    def calculate_average_peaks_distance(peaks):
        """Calculate the average distance between peaks"""
        indexes = [p['peak_index'] for p in peaks]
        if len(indexes) > 2:
            distances = [indexes[i + 1] - indexes[i] for i in range(len(indexes) - 1)]
            return sum(distances) / len(distances)
        else:
            return 0
