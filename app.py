from datetime import datetime
from flask import Flask, render_template, request, url_for, redirect, session, jsonify
from flask_bootstrap import Bootstrap5
from flask_wtf import FlaskForm, CSRFProtect
from werkzeug.utils import secure_filename
from flask_wtf.file import FileField, FileRequired, FileAllowed
from wtforms.validators import DataRequired, Email, NumberRange, ValidationError, Regexp, Length
from wtforms import StringField, IntegerField, SelectField, FloatField, FileField, SubmitField, FormField, FieldList
import configparser
from source import sequence_distribution, visualization
from source.SequenceAnalyzer import SequenceAnalyzer
import json
import plotly
import shutil
from celery import Celery, Task
from celery import shared_task
from celery.result import AsyncResult
import uuid
from pathlib import Path
from flask_mail import Mail, Message
import os
from concurrent.futures import ThreadPoolExecutor
from typing import Dict, Any, List
import logging
import pandas as pd
import secrets
from typing import Optional
import re
from enum import Enum
from dataclasses import dataclass
from pathlib import Path
# from functools import lru_cache

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class DataProcessingError(Exception):
    """Custom exception for data processing errors"""
    pass

class ConfigurationError(Exception):
    """Custom exception for configuration errors"""
    pass

# Configuration
class Config:
    UPLOAD_FOLDER = 'uploads'
    ALLOWED_EXTENSIONS = {'fastq', 'fq'}
    MAX_CONTENT_LENGTH = 500 * 1024 * 1024  # 500MB max file size

def celery_init_app(app: Flask) -> Celery:
    class FlaskTask(Task):
        def __call__(self, *args: object, **kwargs: object) -> object:
            with app.app_context():
                return self.run(*args, **kwargs)

    celery_app = Celery(app.name, task_cls=FlaskTask)
    celery_app.config_from_object(app.config["CELERY"])
    celery_app.set_default()
    app.extensions["celery"] = celery_app
    return celery_app

app = Flask(__name__)
app.secret_key = os.environ.get('FLASK_SECRET_KEY')
app.config['UPLOAD_FOLDER'] = 'static/sessions/'
app.config.from_mapping(
    CELERY=dict(
        broker_url="redis://localhost:6379/0",
        result_backend="redis://localhost:6379/0",
        task_ignore_result=True,
    ),
)
# Bootstrap-Flask requires this line
bootstrap = Bootstrap5(app)
# Flask-WTF requires this line
csrf = CSRFProtect(app)
celery_app = celery_init_app(app)


# Configuration for Flask-Mail
app.config['MAIL_SERVER'] = 'smtp.yandex.ru'
app.config['MAIL_PORT'] = 465
app.config['MAIL_USERNAME'] = os.environ.get('MAIL_USERNAME')
app.config['MAIL_PASSWORD'] = os.environ.get('MAIL_PASSWORD')
app.config['MAIL_USE_TLS'] = False
app.config['MAIL_USE_SSL'] = True

foo = secrets.token_urlsafe(16)
app.secret_key = foo

# Initialize Flask-Mail
mail = Mail(app)

# Utility functions
def allowed_file(filename: str) -> bool:
    """Check if file extension is allowed"""
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in Config.ALLOWED_EXTENSIONS

def ensure_directory_exists(directory: Path) -> None:
    """Ensure directory exists, create if it doesn't"""
    directory.mkdir(parents=True, exist_ok=True)

def process_sequences(form_data: List) -> List[Dict]:
    """Process sequence data from form"""
    return [{
        'type': field['type'],
        'sequence': field['sequence'].replace("\n", "").replace(" ", "")
    } for field in form_data]

def create_parameters_dict(form, filename: str, new_dir: str) -> Dict:
    """Create parameters dictionary from form data"""
    return {
        'fuzzy_similarity': dict(form.fuzzy_similarity.choices).get(form.fuzzy_similarity.data),
        'limit': form.limit.data,
        'threshold': form.threshold.data,
        'filename': filename,
        'new_dir': new_dir,
        'smoothing': dict(form.smoothing.choices).get(form.smoothing.data),
        'email': form.email.data,
        'datetime': str(datetime.now())
    }

def load_output_data(file_path: Path) -> Dict:
    """Load and parse output data from JSON file"""
    try:
        with open(file_path, 'r') as f:
            return json.load(f)
    except Exception as e:
        logger.error(f"Error loading output data: {e}")
        raise DataProcessingError(f"Failed to load output data: {e}")

def process_sequence_data(sequence_data: Dict) -> Dict:
    """Process sequence data for template rendering"""
    return {
        'type': sequence_data['type'],
        'sequence': sequence_data['sequence'],
        'peaks': sequence_data.get('peaks', []),
        'noise_level': sequence_data.get('noise_level', 0),
        'total_reads': sequence_data.get('total_reads', 0),
        'total_proportion': sequence_data.get('total_proportion', 0),
        'value_counts': sequence_data.get('value_counts', [])
    }


# Constants
class SimilarityAlgorithm(Enum):
    LEVENSHTEIN = 'lev'
    BLAST = 'blast'


class SmoothingType(Enum):
    NONE = 'None'
    LOWESS = 'lowess'
    WHITTAKER = 'whittaker'
    SAVGOL = 'savgol'
    CONFSMOOTH = 'confsmooth'


@dataclass
class FormConfig:
    MIN_SEQUENCES: int = 1
    MAX_SEQUENCES: int = 10
    MIN_THRESHOLD: float = 0.1
    MAX_THRESHOLD: float = 1.0
    DEFAULT_THRESHOLD: float = 0.9
    DEFAULT_LIMIT: int = 0
    ALLOWED_FILE_EXTENSIONS: set = frozenset({'fastq'})
    MAX_SESSION_NAME_LENGTH: int = 50
    SESSION_NAME_PATTERN: str = r'^[a-zA-Z0-9_]+$'


class FileValidator:
    """Custom file validator class"""

    @staticmethod
    def validate_fastq(filename: str) -> bool:
        return filename.lower().endswith('.fastq')

    @staticmethod
    def validate_file_size(file_data, max_size_mb: int = 100) -> bool:
        # Check if file size is within limits
        if file_data:
            file_size = len(file_data.read())
            file_data.seek(0)  # Reset file pointer
            return file_size <= max_size_mb * 1024 * 1024
        return False


class SequenceItem(FlaskForm):
    """Form for individual sequence items"""
    type = StringField(
        'Type',
        validators=[
            DataRequired(),
            Length(max=50),
            Regexp(r'^[a-zA-Z0-9_-]+$', message="Type must contain only letters, numbers, hyphen and underscore")
        ]
    )

    sequence = StringField(
        'Sequence',
        validators=[
            DataRequired(),
            Regexp(r'^[ATCG]+$', message="Sequence must contain only valid DNA bases (A, T, C, G)")
        ]
    )

    def validate_sequence(self, field):
        """Validate sequence length and content"""
        if len(field.data) < 10:
            raise ValidationError("Sequence must be at least 10 bases long")
        if len(field.data) > 1000:
            raise ValidationError("Sequence must not exceed 1000 bases")


class InputForm(FlaskForm):
    # @lru_cache(maxsize=128)
    """Main input form with enhanced validation and features"""
    session_name = StringField(
        'Session Name',
        validators=[
            DataRequired(),
            Length(max=FormConfig.MAX_SESSION_NAME_LENGTH),
            Regexp(
                FormConfig.SESSION_NAME_PATTERN,
                message="Session name must contain only letters, numbers, and underscore (_)"
            )
        ]
    )

    items = FieldList(
        FormField(SequenceItem),
        min_entries=FormConfig.MIN_SEQUENCES,
        max_entries=FormConfig.MAX_SEQUENCES
    )

    limit = IntegerField(
        'Limit',
        default=FormConfig.DEFAULT_LIMIT,
        validators=[NumberRange(min=0, message="Limit must be non-negative")]
    )

    fuzzy_similarity = SelectField(
        'Sequence Similarity Algorithm',
        choices=[(algo.value, algo.name.title()) for algo in SimilarityAlgorithm],
        default=SimilarityAlgorithm.LEVENSHTEIN.value
    )

    threshold = FloatField(
        'Threshold',
        validators=[
            DataRequired(),
            NumberRange(
                min=FormConfig.MIN_THRESHOLD,
                max=FormConfig.MAX_THRESHOLD,
                message=f"Threshold must be between {FormConfig.MIN_THRESHOLD} and {FormConfig.MAX_THRESHOLD}"
            )
        ],
        default=FormConfig.DEFAULT_THRESHOLD
    )

    smoothing = SelectField(
        'Smoothing',
        choices=[(smooth.value, smooth.name.title()) for smooth in SmoothingType],
        default=SmoothingType.NONE.value
    )

    email = StringField(
        'Email',
        validators=[
            DataRequired(),
            Email(message="Please enter a valid email address")
        ]
    )

    file = FileField(
        'fastq_file',
        validators=[
            FileRequired(message="Please select a FASTQ file"),
            FileAllowed(
                FormConfig.ALLOWED_FILE_EXTENSIONS,
                message='Only .fastq files are allowed!'
            )
        ]
    )

    submit = SubmitField('Submit')

    def __init__(self, *args, **kwargs):
        """Initialize form with custom validation"""
        super(InputForm, self).__init__(*args, **kwargs)
        self.file_validator = FileValidator()

    def validate_session_name(self, field) -> None:
        """Validate session name and check for existing sessions"""
        # Check if session already exists
        session_path = Path(f"uploads/{field.data}")
        if session_path.exists():
            raise ValidationError("Session name already exists")

    def validate_file(self, field) -> None:
        """Validate file type and size"""
        if field.data:
            if not self.file_validator.validate_fastq(field.data.filename):
                raise ValidationError('File must be in .fastq format')

            if not self.file_validator.validate_file_size(field.data):
                raise ValidationError('File size exceeds maximum limit (100MB)')

    def validate_items(self, field) -> None:
        """Validate sequence items"""
        if len(field.data) < FormConfig.MIN_SEQUENCES:
            raise ValidationError(f"At least {FormConfig.MIN_SEQUENCES} sequence is required")

        # Check for duplicate sequence types
        types = [item['type'] for item in field.data]
        if len(types) != len(set(types)):
            raise ValidationError("Duplicate sequence types are not allowed")

    def clean_data(self) -> dict:
        """Clean and format form data"""
        return {
            'session_name': self.session_name.data,
            'sequences': [{
                'type': item['type'],
                'sequence': item['sequence'].upper().strip()
            } for item in self.items.data],
            'parameters': {
                'limit': self.limit.data,
                'fuzzy_similarity': self.fuzzy_similarity.data,
                'threshold': self.threshold.data,
                'smoothing': self.smoothing.data,
                'email': self.email.data.lower().strip()
            }
        }
#
# class SequenceItem(Form):
#     type = StringField('Type')
#     sequence = StringField('Sequence')
#
# class InputForm(FlaskForm):
#     session_name = StringField('Session Name', validators=[DataRequired()])
#     items = FieldList(FormField(SequenceItem), min_entries=1, max_entries=10)
#     limit = IntegerField('Limit', default=0)
#     fuzzy_similarity = SelectField('Sequence Similarity Algorithm',
#                                     choices=[('lev', 'Levenshtein'), ('blast', 'BLASTn')], default='lev')
#     threshold = FloatField('Threshold', validators=[DataRequired(), NumberRange(min=0.1, max=1.0)], default=0.9)
#     smoothing = SelectField('Smoothing',
#                 choices=[('None','None'),('LOWESS','lowess'),('Whittaker Smoother', 'whittaker'),('savgol','savgol'),('confsmooth','confsmooth')])
#     email = StringField('Email', validators=[DataRequired(), Email()])
#     file = FileField('fastq_file', validators=[FileRequired(), FileAllowed(['fastq'], 'Only .fastq files are allowed!')])
#     submit = SubmitField('Submit')
#
#     # Custom validator for session_name
#     def validate_session_name(self, field):
#         # Only allow latin letters, numbers, and underscore (_)
#         if not re.match(r'^[a-zA-Z0-9_]+$', field.data):
#             raise ValidationError("Session name must contain only letters, numbers, and underscore (_)")
#
#     # Custom validator for file type
#     def validate_file(self, field):
#         # Check that the file has a .fastq extension
#         if field.data:
#             filename = field.data.filename
#             if not filename.lower().endswith('.fastq'):
#                 raise ValidationError('File must be in .fastq format')

def read_config(directory_path):
    config_file = os.path.join(directory_path, 'config.ini')
    config = configparser.ConfigParser()
    config.read(config_file)
    parameters = dict(config['Parameters'])
    sequences = dict(config['Sequences'])
    return parameters, sequences

def process_directories(base_dir):
    all_sessions = []
    for directory in os.listdir(base_dir):
        if os.path.isdir(os.path.join(base_dir, directory)):
            if os.path.isfile(os.path.join(base_dir, directory, 'config.ini')):
                session_id = directory
                parameters, sequences = read_config(os.path.join(base_dir, directory))
                imageURL = os.path.join(base_dir, directory, 'distribution_proportional.png') if \
                    os.path.isfile(os.path.join(base_dir, directory, 'distribution_proportional.png')) \
                    else None
                session_dict = {'sessionID': session_id,
                                'parameters': parameters,
                                'sequences': sequences,
                                'imageURL': imageURL}
                all_sessions.append(session_dict)
    return all_sessions

@app.route('/check')
def check():
    current_directory = os.getcwd()  # Get the current working directory
    return f"Current working directory: {current_directory}"

def generate_config(sequences,
                    parameters):
    config = configparser.ConfigParser()
    config['Sequences'] = {item['type']:item['sequence'] for item in sequences}
    config['Parameters'] = {'fuzzy_similarity': parameters['fuzzy_similarity'],
                            'limit': parameters['limit'],
                            'threshold': parameters['threshold'],
                            'input_file': os.path.join(parameters['new_dir'], parameters['filename']),
                            'smoothing': parameters['smoothing']
                            }
    with open(os.path.join(parameters['new_dir'], 'config.ini'), 'w') as configfile:
        config.write(configfile)


@app.route('/', methods=['GET', 'POST'])
def index():
    """Handle main page requests"""
    form = InputForm()

    if not form.validate_on_submit():
        return render_template('index.html', form=form, page='index')

    try:
        # Process sequences
        sequences = process_sequences(form.items.data)

        # Handle file upload
        file = form.file.data
        if not file or not allowed_file(file.filename):
            raise ValueError("Invalid file type")

        filename = secure_filename(file.filename)
        new_dir = Path(app.config['UPLOAD_FOLDER']) / str(form.session_name.data)

        # Ensure directory exists
        ensure_directory_exists(new_dir)

        # Save file
        file.save(new_dir / filename)

        # Create parameters
        parameters = create_parameters_dict(form, filename, str(new_dir))

        # Generate configuration
        generate_config(sequences, parameters)

        # Store session data
        session['session'] = form.session_name.data
        session['input_data'] = {
            'sequences': {item['type']: item['sequence'] for item in sequences},
            'parameters': parameters,
            'session_name': form.session_name.data
        }

        return redirect(url_for('results'))

    except Exception as e:
        logger.error(f"Error processing form: {e}")
        #flash(f"An error occurred: {str(e)}", 'error')
        return render_template('index.html', form=form, page='index')


@app.route('/contacts')
def contacts():
    return render_template('contacts.html')

@app.route('/sessions')
def sessions():
    base_directory = app.config['UPLOAD_FOLDER']
    all_sessions_list = process_directories(base_directory)

    return render_template('sessions.html', sessions=all_sessions_list, page='sessions')
#
@app.route('/experiment/<sessionID>')
def experiment(sessionID):
    base_directory = app.config['UPLOAD_FOLDER']
    # app.config['UPLOAD_FOLDER']
    try:
        directory_path = os.path.join(base_directory, sessionID)
        parameters, sequences = read_config(directory_path)
        parameters['filename'] = os.path.basename(parameters['input_file'])
        parameters['new_dir'] = directory_path
        data = {'sequences': sequences,
                'parameters': parameters,
                'session_name': sessionID}
        session['input_data'] = data
        return redirect(url_for('results'))
    except Exception as e:
        return redirect(url_for('no_results'))

@app.route('/delete/<sessionID>')
def delete(sessionID):
    base_directory = app.config['UPLOAD_FOLDER']
    directory_path = os.path.join(base_directory, sessionID)
    try:
        shutil.rmtree(directory_path)
        print(f"Directory '{directory_path}' and all its contents have been successfully removed.")
    except Exception as e:
        print(f"An error occurred while removing '{directory_path}': {e}")
    session['input_data'] = None
    return redirect(url_for('sessions'))

def generate_session_id():
    # Generate a random UUID (version 4)
    session_id = uuid.uuid4()
    return str(session_id)


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


def create_merged_dataframe(sequences: list) -> pd.DataFrame:
    """Create merged DataFrame from sequences"""
    try:
        # Get first sequence's occurrences for initialization
        merged_df = sequences[0]['occurrences'][['index']].copy()

        # Dictionary comprehension for sequence mapping
        sequence_dict = {seq['type']: seq['occurrences'] for seq in sequences}

        # Merge all sequences
        for name, df in sequence_dict.items():
            column_mapping = {
                'reads': f'{name}_reads',
                'proportion': f'{name}_proportion',
                'consensus': f'{name}_consensus'
            }
            merged_df = merged_df.merge(
                df.rename(columns=column_mapping),
                on='index',
                how='outer'
            )

        return merged_df
    except Exception as e:
        logger.error(f"Failed to create merged DataFrame: {e}")
        raise DataProcessingError(f"DataFrame merge failed: {e}")


@shared_task(ignore_result=False)
def data_processing(data: Dict[str, Any]) -> Dict[str, str]:
    """Main data processing function"""
    output = {
        'session_id': Path(data['parameters']['new_dir']).name,
        'email': data['parameters']['email']
    }

    try:
        # Initialize analyzer and get output data
        analyzer = SequenceAnalyzer(data['parameters']['new_dir'])
        output_data = analyzer.analyze()

        # Create paths
        base_dir = data['parameters']['new_dir']
        file_paths = {
            'json': os.path.join(base_dir, 'sequences.json'),
            'csv': os.path.join(base_dir, 'occurrences.csv'),
            'prop_plot': os.path.join(base_dir, 'distribution_proportional.png'),
            'abs_plot': os.path.join(base_dir, 'distribution_absolute.png')
        }

        # Generate plots
        plots = {
            'prop_plot': visualization.plot_distribution(
                output_data['sequences'],
                data['parameters']['smoothing'],
                mode='proportion'
            ),
            'abs_plot': visualization.plot_distribution(
                output_data['sequences'],
                data['parameters']['smoothing'],
                mode='reads'
            )
        }

        # Create merged DataFrame
        merged_df = create_merged_dataframe(output_data['sequences'])

        # Use ThreadPoolExecutor for parallel I/O operations
        with ThreadPoolExecutor() as executor:
            # Submit all save tasks
            futures = [
                executor.submit(save_json, file_paths['json'], output_data),
                executor.submit(save_csv, file_paths['csv'], merged_df),
                executor.submit(save_plot, plots['prop_plot'], file_paths['prop_plot']),
                executor.submit(save_plot, plots['abs_plot'], file_paths['abs_plot'])
            ]

            # Wait for all tasks to complete and check for exceptions
            for future in futures:
                future.result()

        logger.info(f"Successfully processed data for session {output['session_id']}")

    except Exception as e:
        logger.error(f"Error processing data for session {output['session_id']}: {e}")

    finally:
        return output

@app.route('/no_results')
def no_results():
    return render_template('no_results.html')

@app.route('/results')
def results():
    """Handle results page requests"""
    if 'input_data' not in session or session['input_data'] is None:
        return render_template('no_results.html', page='results')

    try:
        data = session['input_data']
        sequences_file = Path(data['parameters']['new_dir']) / 'sequences.json'

        if not sequences_file.exists():
            # Start async processing
            result = data_processing.delay(data)
            return render_template('async_result.html',
                                result_id=result.id,
                                parameters=data['parameters'],
                                page='results')

        # Load and process output data
        output_data = load_output_data(sequences_file)

        # Generate plots
        plots = {
            'hist1': json.dumps(
                visualization.plot_distribution(output_data['sequences'],
                                             data['parameters']['smoothing'],
                                             mode='proportion'),
                cls=plotly.utils.PlotlyJSONEncoder
            ),
            'hist2': json.dumps(
                visualization.plot_distribution(output_data['sequences'],
                                             data['parameters']['smoothing'],
                                             mode='reads'),
                cls=plotly.utils.PlotlyJSONEncoder
            )
        }

        # Process sequences
        sequences = [process_sequence_data(seq) for seq in output_data['sequences']]

        # Extract fastq parameters
        fastq_parameters = {
            'n_records': output_data['parameters']['n_records'],
            'avg_noise_level': output_data['parameters']['avg_noise_level']
        }

        return render_template('results.html',
                             plots=plots,
                             data=data,
                             sequences=sequences,
                             fastq_parameters=fastq_parameters,
                             page='results')

    except Exception as e:
        logger.error(f"Error processing results: {e}")
        #flash("An error occurred while processing results", 'error')
        return redirect(url_for('index'))

def send_email(email, sessionID, path):
    recipient = email
    subject = f'NanoporeInspect: session {sessionID} results are ready'
    message_body = f'The results are available in the web application by the link: {path}'
    msg = Message(
        subject=subject,
        sender=app.config['MAIL_USERNAME'],
        recipients=[recipient]
    )
    msg.body = message_body  # Plain text email body

    try:
        mail.send(msg)
        return f"Email sent to {recipient}!"
    except Exception as e:
        return f"Failed to send email. Error: {str(e)}"

@app.route('/send_email')
def send_email_test():
    recipient = 'magsend@gmail.com'
    subject = 'TEST from NanoporeInspect'
    message_body = 'Test from NanoporeInspect'
    sender = app.config['MAIL_USERNAME']
    print(sender)

    msg = Message(
        subject=subject,
        sender=app.config['MAIL_USERNAME'],
        recipients=[recipient]
    )
    msg.body = message_body  # Plain text email body

    try:
        mail.send(msg)
        return f"Email sent to {recipient}!"
    except Exception as e:
        return f"Failed to send email. Error: {str(e)}"

@app.route("/result/<id>", methods=['GET','POST'])
def task_result(id: str) -> dict[str, object]:
    result = AsyncResult(id)
    if result.ready():
        session_id = result.result['session_id']
        email = result.result['email']
        path = request.url_root + 'experiment/' + session_id
        send_email(email, session_id, path)
        return redirect(url_for('experiment', sessionID=session_id))
    else:
        return render_template('in_progress.html', result_id = id, page='results')

    # return {
    #     "ready": result.ready(),
    #     "successful": result.successful(),
    #     "value": result.result if result.ready() else None,
    # }

# Error handlers
@app.errorhandler(404)
def not_found_error(error):
    return render_template('404.html'), 404

@app.errorhandler(500)
def internal_error(error):
    return render_template('500.html'), 500


if __name__ == '__main__':
    app.run(debug=True)