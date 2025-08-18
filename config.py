# config.py
import os
import secrets
from pathlib import Path
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

class BaseConfig:
    """Base configuration class"""
    SECRET_KEY = os.getenv('FLASK_SECRET_KEY', secrets.token_urlsafe(16))
    redis_host = os.getenv("REDIS_HOST", "127.0.0.1")
    redis_port = os.getenv("REDIS_PORT", 6379)
    BROKER_URL = f"redis://{redis_host}:{redis_port}/0"
    BACKEND_URL = f"redis://{redis_host}:{redis_port}/0"

    # File upload configuration
    #UPLOAD_FOLDER = 'static/sessions/'
    ALLOWED_EXTENSIONS = {'fastq', 'fq'}
    MAX_CONTENT_LENGTH = 5 * 1024 ** 3  # 1GB max file size
    # Increase request timeouts
    SEND_FILE_MAX_AGE_DEFAULT = 0

    # Celery configuration
    CELERY = {
        "broker_url": BROKER_URL,
        "result_backend": BACKEND_URL,
        "task_ignore_result": True,
    }

    UPLOAD_FOLDER = os.getenv('RESULTS_PATH')

    @staticmethod
    def init_app(app):
        """Initialize application configuration"""
        # Create upload directory if it doesn't exist
        upload_path = Path(app.config['UPLOAD_FOLDER'])
        upload_path.mkdir(parents=True, exist_ok=True)


class DevelopmentConfig(BaseConfig):
    """Development configuration"""
    DEBUG = True
    TESTING = False
    # Add development-specific settings here


class ProductionConfig(BaseConfig):
    """Production configuration"""
    DEBUG = False
    TESTING = False

    # Override with more secure production settings
    SECRET_KEY = os.environ.get('PRODUCTION_SECRET_KEY')

    redis_host = os.getenv("REDIS_HOST", "127.0.0.1")
    redis_url = f"redis://{redis_host}:6379/0"

    # Production Celery settings
    CELERY = {
        "broker_url": os.environ.get('CELERY_BROKER_URL', redis_url),
        "result_backend": os.environ.get('CELERY_RESULT_BACKEND', redis_url),
        "task_ignore_result": True,
    }


class TestingConfig(BaseConfig):
    """Testing configuration"""
    TESTING = True
    DEBUG = True
    # Add testing-specific settings here


# Configuration dictionary
config = {
    'development': DevelopmentConfig,
    'production': ProductionConfig,
    'testing': TestingConfig,
    'default': DevelopmentConfig
}