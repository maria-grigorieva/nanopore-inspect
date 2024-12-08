from typing import Dict, Any, List
import logging
import json
import uuid
from pathlib import Path

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class DataProcessingError(Exception):
    """Custom exception for data processing errors"""
    pass
def save_json(file_path: str, data: Dict[str, Any]) -> None:
    """Save data to JSON file"""
    try:
        with open(file_path, 'w') as f:
            json.dump(data, f, default=str)
    except Exception as e:
        logger.error(f"Failed to save JSON file: {e}")
        raise DataProcessingError(f"JSON save failed: {e}")


def save_csv(file_path: str, df) -> None:
    """Save DataFrame to CSV file"""
    try:
        df.to_csv(file_path)
    except Exception as e:
        logger.error(f"Failed to save CSV file: {e}")
        raise DataProcessingError(f"CSV save failed: {e}")


def save_plot(fig, file_path: str) -> None:
    """Save plot to file"""
    try:
        with open(file_path, "wb") as f:
            fig.write_image(f)
    except Exception as e:
        logger.error(f"Failed to save plot: {e}")
        raise DataProcessingError(f"Plot save failed: {e}")

def generate_session_id():
    # Generate a random UUID (version 4)
    session_id = uuid.uuid4()
    return str(session_id)

def load_output_data(file_path: Path) -> Dict:
    """Load and parse output data from JSON file"""
    try:
        with open(file_path, 'r') as f:
            return json.load(f)
    except Exception as e:
        logger.error(f"Error loading output data: {e}")
        raise DataProcessingError(f"Failed to load output data: {e}")

def ensure_directory_exists(directory: Path) -> None:
    """Ensure directory exists, create if it doesn't"""
    directory.mkdir(parents=True, exist_ok=True)
