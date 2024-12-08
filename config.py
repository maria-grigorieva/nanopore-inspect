# config.py
import os
import secrets
from pathlib import Path

class BaseConfig:
    """Base configuration class"""
    # Basic Flask configuration
    SECRET_KEY = os.environ.get('FLASK_SECRET_KEY', secrets.token_urlsafe(16))

    # File upload configuration
    UPLOAD_FOLDER = 'static/sessions/'
    ALLOWED_EXTENSIONS = {'fastq', 'fq'}
    MAX_CONTENT_LENGTH = 500 * 1024 * 1024  # 500MB max file size

    # Celery configuration
    CELERY = {
        "broker_url": "redis://localhost:6379/0",
        "result_backend": "redis://localhost:6379/0",
        "task_ignore_result": True,
    }

    # Mail configuration
    MAIL_SERVER = 'smtp.yandex.ru'
    MAIL_PORT = 465
    MAIL_USERNAME = os.environ.get('MAIL_USERNAME')
    MAIL_PASSWORD = os.environ.get('MAIL_PASSWORD')
    MAIL_USE_TLS = False
    MAIL_USE_SSL = True

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

    # Production Celery settings
    CELERY = {
        "broker_url": os.environ.get('CELERY_BROKER_URL', "redis://localhost:6379/0"),
        "result_backend": os.environ.get('CELERY_RESULT_BACKEND', "redis://localhost:6379/0"),
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