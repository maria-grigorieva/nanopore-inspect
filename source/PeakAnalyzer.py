from scipy.signal import find_peaks, peak_prominences, peak_widths
import numpy as np


class PeakAnalyzer:
    """Class for handling peak detection and analysis operations"""

    def __init__(self, window_size=20):
        self.window_size = window_size

    def find_peak_boundaries(self, smoothed_array, peak_indices):
        """
        Find more accurate peak boundaries by analyzing the slope changes
        around each peak.
        """
        peak_info = []

        for peak_idx in peak_indices:
            # Define initial search ranges
            left_limit = max(0, peak_idx - self.window_size)
            right_limit = min(len(smoothed_array), peak_idx + self.window_size)

            # Find left base by looking for slope change
            left_segment = smoothed_array[left_limit:peak_idx]
            left_base = left_limit
            if len(left_segment) > 0:
                slopes = np.diff(left_segment)
                for i in range(len(slopes) - 1):
                    if slopes[i] < slopes[i + 1]:
                        left_base = left_limit + i
                        break

            # Find right base by looking for slope change
            right_segment = smoothed_array[peak_idx:right_limit]
            right_base = right_limit
            if len(right_segment) > 0:
                slopes = np.diff(right_segment)
                for i in range(len(slopes) - 1):
                    if slopes[i] < 0 and slopes[i + 1] > 0:
                        right_base = peak_idx + i
                        break

            peak_info.append({
                'peak_index': peak_idx,
                'left_base': left_base,
                'right_base': right_base
            })

        return peak_info

    def analyze_peaks(self, reference_type,
                      reference_value,
                      smoothed_array,
                      reads_array,
                      proportion_array,
                      consensus_array):
        """Main method for peak analysis"""
        # Initialize default return values
        list_of_peaks = []
        average_peaks_distance = 0

        try:
            # Find initial peaks
            peaks, initial_bases = find_peaks(smoothed_array)

            # If no initial peaks found, return default values
            if len(peaks) == 0:
                return list_of_peaks, average_peaks_distance

            # Calculate prominences
            prominences = peak_prominences(smoothed_array, peaks)[0]
            avg_prominence = np.mean(prominences) if len(prominences) > 0 else 0

            # Find peaks with adjusted parameters
            peak_indices, properties = find_peaks(smoothed_array,
                                                  prominence=avg_prominence,
                                                  distance=len(reference_value),
                                                  width=5,
                                                  rel_height=0.7)

            # If no peaks found with adjusted parameters, use initial peaks
            if len(peak_indices) == 0:
                peak_indices = peaks

            # Get peak boundaries
            peak_info = self.find_peak_boundaries(smoothed_array, peak_indices)

            # If peak boundaries found, process peaks
            if peak_info:
                extremums = self._process_peaks(reads_array,
                                                proportion_array,
                                                consensus_array,
                                                peak_indices,
                                                peak_info)
                list_of_peaks = extremums if len(extremums) > 0 else []
                if list_of_peaks:
                    average_peaks_distance = self.calculate_average_peaks_distance(list_of_peaks)

        except Exception as e:
            print(f"Error during peak analysis: {e}")
            # Log the error if you have logging configured
            # logger.error(f"Peak analysis failed: {e}")

        # Always return a tuple of list and number
        return list_of_peaks, average_peaks_distance

    def _process_peaks(self, reads_array,
                       proportion_array,
                       consensus_array,
                       peak_indices,
                       peak_info):
        """Process all peaks and calculate their properties"""
        extremums = []
        for i in range(len(peak_indices)):
            peak_index = peak_indices[i]
            extremums.append(self.aggregate_peak_values(i, reads_array,
                                                        proportion_array,
                                                        consensus_array,
                                                        peak_index,
                                                        peak_info))

        for i in range(1, len(extremums)):
            extremums[i]['peak_dist'] = extremums[i]['peak_index'] - extremums[i - 1]['peak_index']

        return extremums

    @staticmethod
    def aggregate_peak_values(step, reads_array,
                              proportion_array, consensus_array, peak_index, peak_info):
        """Calculate aggregate values for a single peak"""
        left_bases = peak_info[step]['left_base']
        right_bases = peak_info[step]['right_base']
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
