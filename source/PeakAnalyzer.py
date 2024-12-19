from scipy.signal import find_peaks, peak_prominences, peak_widths
import numpy as np
from typing import List, Tuple, Dict, Any


class PeakAnalyzer:
    def __init__(self, window_size: int = 20, min_distance_factor: float = 0.8):
        self.window_size = window_size
        self.min_distance_factor = min_distance_factor

    def _find_initial_peaks(self,
                            smoothed_array: np.ndarray,
                            reference_length: int) -> Tuple[np.ndarray, Dict]:
        """
        Find initial peaks with minimal restrictions
        """
        # Calculate basic statistics
        signal_mean = np.mean(smoothed_array)
        signal_std = np.std(smoothed_array)
        min_distance = int(reference_length * self.min_distance_factor)

        # Try different height thresholds until we find peaks
        height_thresholds = [
            None,  # No height restriction
            signal_mean * 0.5,  # Very permissive
            signal_mean * 0.7,  # Moderate
            signal_mean  # Stricter
        ]

        for height in height_thresholds:
            peaks, properties = find_peaks(
                smoothed_array,
                height=height,
                distance=min_distance  # This is the only strict requirement
            )

            if len(peaks) > 0:
                # Calculate prominences for found peaks
                prominences = peak_prominences(smoothed_array, peaks)[0]

                # Keep only peaks with prominence above minimal threshold
                min_prominence = np.mean(prominences) * 0.3  # Very permissive threshold
                good_peaks_mask = prominences > min_prominence

                if np.any(good_peaks_mask):
                    return peaks[good_peaks_mask], properties

        # If no peaks found, return empty arrays
        return np.array([]), {}

    def find_peak_boundaries(self, smoothed_array: np.ndarray, peak_indices: np.ndarray) -> List[Dict[str, int]]:
        """
        Find peak boundaries using a simplified approach
        """
        peak_info = []

        for peak_idx in peak_indices:
            left_limit = max(0, peak_idx - self.window_size)
            right_limit = min(len(smoothed_array), peak_idx + self.window_size)

            # Find left boundary - look for minimum value in left window
            left_segment = smoothed_array[left_limit:peak_idx]
            left_base = left_limit + np.argmin(left_segment) if len(left_segment) > 0 else left_limit

            # Find right boundary - look for minimum value in right window
            right_segment = smoothed_array[peak_idx:right_limit]
            right_base = peak_idx + np.argmin(right_segment) if len(right_segment) > 0 else right_limit

            peak_info.append({
                'peak_index': peak_idx,
                'left_base': left_base,
                'right_base': right_base
            })

        return peak_info

    def analyze_peaks(self,
                      reference_type: str,
                      reference_value: str,
                      smoothed_array: np.ndarray,
                      reads_array: np.ndarray,
                      proportion_array: np.ndarray,
                      consensus_array: np.ndarray) -> Tuple[List[Any], float]:
        """
        Main method for peak analysis with simplified approach
        """
        list_of_peaks = []
        average_peaks_distance = 0

        try:
            # Find peaks with minimal restrictions
            peak_indices, properties = self._find_initial_peaks(
                smoothed_array,
                len(reference_value)
            )

            if len(peak_indices) > 0:
                # Get peak boundaries
                peak_info = self.find_peak_boundaries(smoothed_array, peak_indices)

                # Process the peaks
                extremums = self._process_peaks(
                    reads_array,
                    proportion_array,
                    consensus_array,
                    peak_indices,
                    peak_info
                )

                list_of_peaks = extremums if len(extremums) > 0 else []
                if list_of_peaks:
                    average_peaks_distance = self.calculate_average_peaks_distance(list_of_peaks)

        except Exception as e:
            print(f"Error during peak analysis: {e}")
            print(f"Error details: {str(e)}")
            # You might want to add proper logging here

        return list_of_peaks, average_peaks_distance
# from scipy.signal import find_peaks, peak_prominences, peak_widths
# import numpy as np
#
#
# class PeakAnalyzer:
#     """Class for handling peak detection and analysis operations"""
#
#     def __init__(self, window_size=20):
#         self.window_size = window_size
#
#     def find_peak_boundaries(self, smoothed_array, peak_indices):
#         """
#         Find more accurate peak boundaries by analyzing the slope changes
#         around each peak.
#         """
#         peak_info = []
#
#         for peak_idx in peak_indices:
#             # Define initial search ranges
#             left_limit = max(0, peak_idx - self.window_size)
#             right_limit = min(len(smoothed_array), peak_idx + self.window_size)
#
#             # Find left base by looking for slope change
#             left_segment = smoothed_array[left_limit:peak_idx]
#             left_base = left_limit
#             if len(left_segment) > 0:
#                 slopes = np.diff(left_segment)
#                 for i in range(len(slopes) - 1):
#                     if slopes[i] < slopes[i + 1]:
#                         left_base = left_limit + i
#                         break
#
#             # Find right base by looking for slope change
#             right_segment = smoothed_array[peak_idx:right_limit]
#             right_base = right_limit
#             if len(right_segment) > 0:
#                 slopes = np.diff(right_segment)
#                 for i in range(len(slopes) - 1):
#                     if slopes[i] < 0 and slopes[i + 1] > 0:
#                         right_base = peak_idx + i
#                         break
#
#             peak_info.append({
#                 'peak_index': peak_idx,
#                 'left_base': left_base,
#                 'right_base': right_base
#             })
#
#         return peak_info
#
#     def analyze_peaks(self, reference_type,
#                       reference_value,
#                       smoothed_array,
#                       reads_array,
#                       proportion_array,
#                       consensus_array):
#         """Main method for peak analysis"""
#         # Initialize default return values
#         list_of_peaks = []
#         average_peaks_distance = 0
#
#         try:
#             # Find initial peaks
#             peaks, initial_bases = find_peaks(smoothed_array)
#
#             # If no initial peaks found, return default values
#             if len(peaks) == 0:
#                 return list_of_peaks, average_peaks_distance
#
#             # Calculate prominences
#             prominences = peak_prominences(smoothed_array, peaks)[0]
#             avg_prominence = np.mean(prominences) if len(prominences) > 0 else 0
#
#             # Find peaks with adjusted parameters
#             peak_indices, properties = find_peaks(smoothed_array,
#                                                   prominence=avg_prominence,
#                                                   distance=len(reference_value),
#                                                   width=5,
#                                                   rel_height=0.7)
#
#             # If no peaks found with adjusted parameters, use initial peaks
#             if len(peak_indices) == 0:
#                 peak_indices = peaks
#
#             # Get peak boundaries
#             peak_info = self.find_peak_boundaries(smoothed_array, peak_indices)
#
#             # If peak boundaries found, process peaks
#             if peak_info:
#                 extremums = self._process_peaks(reads_array,
#                                                 proportion_array,
#                                                 consensus_array,
#                                                 peak_indices,
#                                                 peak_info)
#                 list_of_peaks = extremums if len(extremums) > 0 else []
#                 if list_of_peaks:
#                     average_peaks_distance = self.calculate_average_peaks_distance(list_of_peaks)
#
#         except Exception as e:
#             print(f"Error during peak analysis: {e}")
#             # Log the error if you have logging configured
#             # logger.error(f"Peak analysis failed: {e}")
#
#         # Always return a tuple of list and number
#         return list_of_peaks, average_peaks_distance

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
