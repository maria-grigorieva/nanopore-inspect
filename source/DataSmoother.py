import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess
from whittaker_eilers import WhittakerSmoother
from statsmodels.nonparametric.kernel_regression import KernelReg
from confsmooth import confsmooth
from scipy.stats import variation
from scipy.signal import savgol_filter

class DataSmoother:
    """Class handling all data smoothing operations"""

    @staticmethod
    def moving_average(arr, window_size):
        """Calculate moving average of an array."""
        return np.convolve(arr, np.ones(window_size) / window_size, mode='valid')

    @staticmethod
    def calculate_noise(arr, window_size=3):
        """Calculate noise in the data."""
        smoothed = DataSmoother.moving_average(arr, window_size)
        extended_smoothed = np.concatenate((arr[:window_size - 1], smoothed))
        return arr - extended_smoothed

    @staticmethod
    def signaltonoise(a, axis=0, ddof=0):
        """Calculate signal-to-noise ratio."""
        a = np.asanyarray(a)
        m = a.mean(axis)
        sd = a.std(axis=axis, ddof=ddof)
        return np.where(sd == 0, 0, m / sd)

    @staticmethod
    def noise_level(data, reference):
        """Calculate noise level using coefficient of variation."""
        return variation(data['reads'], axis=0)

    @staticmethod
    def get_noise_level(array):
        return float(np.round(DataSmoother.signaltonoise(array), 4))

    def smooth_data(self, array, reference, smooth_type='whittaker'):
        """Apply various smoothing methods to the data."""
        # target = 'proportion'
        self.array = array
        self.reference = reference
        if smooth_type.lower() != 'none':
            # target = 'smoothed'
            if smooth_type.lower() == 'whittaker':
                return self._apply_whittaker_smoothing()
            elif smooth_type.lower() == 'lowess':
                return self._apply_lowess_smoothing()
            elif smooth_type.lower() == 'savgol':
                return self._apply_savgol_smoothing()
            elif smooth_type.lower() == 'confsmooth':
                return self._apply_confsmooth()
        else:
            return array

    def _apply_whittaker_smoothing(self):
        """Apply Whittaker smoothing."""
        whittaker_smoother = WhittakerSmoother(
            lmbda=len(self.array), order=2, data_length=len(self.array)
        )
        return whittaker_smoother.smooth(self.array)

    def _apply_lowess_smoothing(self):
        """Apply LOWESS smoothing."""
        return lowess(self.array, range(len(self.array)), frac=0.1)[:, 1]

    def _apply_savgol_smoothing(self):
        """Apply Savitzky-Golay smoothing."""
        return savgol_filter(self.array,
                                   window_length=len(self.reference),
                                   polyorder=2)

    def _apply_confsmooth(self):
        """Apply confidence smoothing."""
        print('Estimated noise level =', self.noise_level(self.array, self.reference))
        return confsmooth(self.array, self.noise_level, confidence=0.9, deg=2)