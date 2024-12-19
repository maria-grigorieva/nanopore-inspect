# app/forms/validators.py
class FileValidator:
    """Custom file validator class"""

    @staticmethod
    def validate_fastq(filename: str) -> bool:
        return filename.lower().endswith('.fastq')

    @staticmethod
    def validate_file_size(file_data, max_size_mb: int = 1000) -> bool:
        if file_data:
            file_size = len(file_data.read())
            file_data.seek(0)
            return file_size <= max_size_mb * 1024 * 1024
        return False