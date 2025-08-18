# app/forms/constants.py
from enum import Enum
from dataclasses import dataclass
from typing import FrozenSet

class SimilarityAlgorithm(Enum):
    LEVENSHTEIN = 'lev'
    BLAST = 'blast'

class SmoothingType(Enum):
    NONE = 'None'
    LOWESS = 'lowess'
    WHITTAKER = 'whittaker'
    SAVGOL = 'savgol'
    CONFSMOOTH = 'confsmooth'

@dataclass(frozen=True)
class FormConfig:
    MIN_SEQUENCES: int = 1
    MAX_SEQUENCES: int = 10
    MIN_THRESHOLD: float = 0.1
    MAX_THRESHOLD: float = 1.0
    DEFAULT_THRESHOLD: float = 0.9
    DEFAULT_LIMIT: int = 0
    ALLOWED_FILE_EXTENSIONS: FrozenSet[str] = frozenset({'fastq'})
    MAX_SESSION_NAME_LENGTH: int = 50
    SESSION_NAME_PATTERN: str = r'^[a-zA-Z0-9_]+$'