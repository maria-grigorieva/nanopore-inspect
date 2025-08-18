# factory.py
from flask import Flask
from celery import Celery, Task
from flask_bootstrap import Bootstrap5
from flask_wtf import CSRFProtect
from flask_mail import Mail


def celery_init_app(app: Flask) -> Celery:
    """Initialize Celery instance"""

    class FlaskTask(Task):
        def __call__(self, *args: object, **kwargs: object) -> object:
            with app.app_context():
                return self.run(*args, **kwargs)

    celery_app = Celery(app.name, task_cls=FlaskTask)
    celery_app.config_from_object(app.config["CELERY"])
    celery_app.set_default()
    app.extensions["celery"] = celery_app
    return celery_app


def create_app(config_name='default'):
    """Application factory function"""
    app = Flask(__name__)

    # Load configuration
    from config import config
    app.config.from_object(config[config_name])
    config[config_name].init_app(app)

    # Initialize extensions
    bootstrap = Bootstrap5(app)
    csrf = CSRFProtect(app)
    celery_app = celery_init_app(app)
    mail = Mail(app)

    # Store extensions in app.extensions
    app.extensions['bootstrap'] = bootstrap
    app.extensions['csrf'] = csrf
    app.extensions['mail'] = mail

    return app, celery_app