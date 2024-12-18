# Standard library imports
import configparser
import json
import logging
import os
import secrets
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Dict, Any

# Third-party imports
import pandas as pd
import plotly
from celery import Celery, Task, shared_task
from celery.result import AsyncResult
from flask import Flask, render_template, request, url_for, redirect, session, jsonify
from flask_bootstrap import Bootstrap5
from flask_mail import Mail, Message
from flask_wtf import CSRFProtect
from werkzeug.utils import secure_filename

# Local application imports
from config import config
from forms import InputForm
from source import visualization
from source.SequenceAnalyzer import SequenceAnalyzer
from utils import (
    save_csv,
    save_json,
    save_plot,
    load_output_data,
    ensure_directory_exists,
    remove_session_dir,
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class DataProcessingError(Exception):
    """Custom exception for data processing errors"""
    pass


class ConfigurationError(Exception):
    """Custom exception for configuration errors"""
    pass


# Get environment configuration
# env = os.environ.get('FLASK_ENV', 'default')
# app, celery_app = create_app(env)

# Make the app and celery instances available at module level
# celery = celery_app
#
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
app.config.from_object(config['default'])
config['default'].init_app(app)
# Bootstrap-Flask requires this line
bootstrap = Bootstrap5(app)
# Flask-WTF requires this line
csrf = CSRFProtect(app)
celery_app = celery_init_app(app)

foo = secrets.token_urlsafe(16)
app.secret_key = foo
#
# Initialize Flask-Mail
mail = Mail(app)


# Utility functions
def allowed_file(filename: str) -> bool:
    """Check if file extension is allowed"""
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in app.config['ALLOWED_EXTENSIONS']


def process_sequence_data(sequence_data: Dict) -> Dict:
    """Process sequence data for template rendering"""
    return {
        'type': sequence_data.get('type', 'unknown'),
        'sequence': sequence_data.get('sequence', ''),
        'peaks': sequence_data.get('peaks', []),
        'noise_level': sequence_data.get('noise_level', 0.0),
        'total_reads': sequence_data.get('total_reads', 0),
        'total_proportion': sequence_data.get('total_proportion', 0.0),
        'value_counts': sequence_data.get('value_counts', []),
    }


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
    config['Sequences'] = {item['type']: item['sequence'] for item in sequences}
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
        sequences = form.process_sequences()

        # Handle file upload
        file = form.file.data
        if not file or not allowed_file(file.filename):
            raise ValueError("Invalid file type")

        filename = secure_filename(file.filename)
        # Create session directory
        new_dir = Path(app.config['UPLOAD_FOLDER']) / str(form.session_name.data)

        # Ensure directory exists
        ensure_directory_exists(new_dir)

        # Save file
        file.save(new_dir / filename)

        # Create parameters
        parameters = form.create_parameters_dict(filename, str(new_dir))

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
        # flash(f"An error occurred: {str(e)}", 'error')
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
    remove_session_dir(base_directory, sessionID)
    session['input_data'] = None
    return redirect(url_for('sessions'))


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
        # flash("An error occurred while processing results", 'error')
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


@app.route("/result/<id>", methods=['GET', 'POST'])
def task_result(id: str) -> object:
    """
    Handle the result of an asynchronous task.

    Args:
        id (str): Task ID.

    Returns:
        Response: Redirect to the experiment page if ready, or renders an 'in progress' template.
    """
    try:
        # Get task result
        result = AsyncResult(id)

        if result.ready():
            # Extract task results
            task_data = result.result or {}
            session_id = task_data.get('session_id')
            email = task_data.get('email')

            if not session_id or not email:
                logging.error(f"Task {id} result is incomplete: {task_data}")
                return jsonify({"error": "Task result is incomplete."}), 500

            # Construct path and send email
            experiment_path = request.url_root + 'experiment/' + session_id
            try:
                send_email(email, session_id, experiment_path)
                logging.info(f"Email sent to {email} for session {session_id}.")
            except Exception as e:
                logging.error(f"Failed to send email for session {session_id}: {e}")
                return jsonify({"error": "Failed to send email."}), 500

            # Redirect to experiment page
            return redirect(url_for('experiment', sessionID=session_id))
        else:
            # Task is still in progress
            return render_template('in_progress.html', result_id=id, page='results')

    except Exception as e:
        logging.error(f"Error processing task {id}: {e}")
        return jsonify({"error": "An error occurred while processing the task."}), 500


# Error handlers
@app.errorhandler(404)
def not_found_error(error):
    return render_template('404.html'), 404


@app.errorhandler(500)
def internal_error(error):
    return render_template('500.html'), 500


if __name__ == '__main__':
    app.run(debug=True)
