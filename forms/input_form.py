# app/forms/input_form.py
from flask_wtf import FlaskForm
from flask_wtf.file import FileRequired, FileAllowed
from wtforms import (
    StringField, IntegerField, SelectField,
    FloatField, FileField, SubmitField,
    FormField, FieldList
)
from wtforms.validators import (
    DataRequired, Email, NumberRange,
    ValidationError, Regexp, Length
)
from pathlib import Path

from .constants import SimilarityAlgorithm, SmoothingType, FormConfig
from .validators import FileValidator


class SequenceItem(FlaskForm):
    """Form for individual sequence items"""
    type = StringField(
        'Type',
        validators=[
            DataRequired(),
            Length(max=50),
            Regexp(r'^[a-zA-Z0-9_-]+$',
                   message="Type must contain only letters, numbers, hyphen and underscore")
        ]
    )

    sequence = StringField(
        'Sequence',
        validators=[
            DataRequired(),
            Regexp(r'^[ATCG]+$',
                   message="Sequence must contain only valid DNA bases (A, T, C, G)")
        ]
    )

    def validate_sequence(self, field):
        if len(field.data) < 10:
            raise ValidationError("Sequence must be at least 10 bases long")
        if len(field.data) > 100:
            raise ValidationError("Sequence must not exceed 100 bases")


class InputForm(FlaskForm):
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
        super(InputForm, self).__init__(*args, **kwargs)
        self.file_validator = FileValidator()

    def validate_session_name(self, field):
        session_path = Path(f"uploads/{field.data}")
        if session_path.exists():
            raise ValidationError("Session name already exists")

    def validate_file(self, field):
        if field.data:
            if not self.file_validator.validate_fastq(field.data.filename):
                raise ValidationError('File must be in .fastq format')

            if not self.file_validator.validate_file_size(field.data):
                raise ValidationError('File size exceeds maximum limit (100MB)')

    def validate_items(self, field):
        if len(field.data) < FormConfig.MIN_SEQUENCES:
            raise ValidationError(f"At least {FormConfig.MIN_SEQUENCES} sequence is required")

        types = [item['type'] for item in field.data]
        if len(types) != len(set(types)):
            raise ValidationError("Duplicate sequence types are not allowed")

    def clean_data(self):
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