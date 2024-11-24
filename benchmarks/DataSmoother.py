from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.signal import savgol_filter
import numpy as np
from whittaker_eilers import WhittakerSmoother
from confsmooth import confsmooth
from scipy.stats import variation

class DataSmoother:
    """A class for smoothing numerical arrays using various smoothing algorithms."""

    SMOOTHING_TYPES = ['whittaker', 'lowess', 'savgol', 'confsmooth']

    def __init__(self):
        """Initialize the DataSmoother class."""
        self.result = {}

    def validate_input(self, data):
        """
        Validate input data.

        Args:
            data: Input array-like object

        Raises:
            ValueError: If data is empty or contains non-numeric values
            TypeError: If data is not array-like
        """
        try:
            if len(data) == 0:
                raise ValueError("Input data cannot be empty")

            # Convert to numpy array if not already
            data = np.array(data, dtype=float)

            if np.isnan(data).any():
                raise ValueError("Input data contains NaN values")

            return data

        except TypeError:
            raise TypeError("Input must be array-like")

    def whittaker_smooth(self, data):
        """Apply Whittaker smoothing."""
        whittaker_smoother = WhittakerSmoother(
            lmbda=len(data),
            order=2,
            data_length=len(data)
        )
        return whittaker_smoother.smooth(data)

    def lowess_smooth(self, data):
        """Apply LOWESS smoothing."""
        return lowess(
            endog=data,
            exog=range(len(data)),
            frac=0.15,
            it=3,
            return_sorted=False
        )

    def savgol_smooth(self, data):
        """Apply Savitzky-Golay smoothing."""
        window_length = 23  # Previously hardcoded as len("TCGATTCCGTTTGTAGTCGTCTGT")
        return savgol_filter(
            data,
            window_length=window_length,
            polyorder=2
        )

    def confsmooth_smooth(self, data):
        """Apply Confidence Interval smoothing."""
        return confsmooth(
            data,
            variation(data),
            confidence=0.9,
            deg=2
        )

    def smooth(self, data, sequence_type, smoothing_types=None):
        """
        Smooth the input data using specified smoothing methods.

        Args:
            data: Array-like object containing numerical data
            sequence_type: String identifier for the type of sequence
            smoothing_types: List of smoothing methods to apply (default: all methods)

        Returns:
            dict: Dictionary containing smoothed data for each method

        Raises:
            ValueError: If invalid smoothing type is specified
        """
        # Validate and prepare input data
        data = self.validate_input(data)

        # Use all smoothing types if none specified
        if smoothing_types is None:
            smoothing_types = self.SMOOTHING_TYPES

        # Validate smoothing types
        invalid_types = [t for t in smoothing_types if t.lower() not in self.SMOOTHING_TYPES]
        if invalid_types:
            raise ValueError(f"Invalid smoothing type(s): {invalid_types}")

        # Dictionary mapping smoothing types to their corresponding methods
        smoothing_methods = {
            'whittaker': self.whittaker_smooth,
            'lowess': self.lowess_smooth,
            'savgol': self.savgol_smooth,
            'confsmooth': self.confsmooth_smooth
        }

        # Apply each smoothing method
        result = {}
        for smoothing_type in smoothing_types:
            try:
                smooth_method = smoothing_methods[smoothing_type.lower()]
                result[f'{smoothing_type}_{sequence_type}_smooth'] = smooth_method(data)
            except Exception as e:
                print(f"Warning: Failed to apply {smoothing_type} smoothing: {str(e)}")

        return result
