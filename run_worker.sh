#!/bin/bash

# Set default concurrency if not provided
export CELERY_CONCURRENCY="${CELERY_CONCURRENCY:-4}"

# Run the Celery worker
celery -A app.celery_app worker --loglevel=INFO --concurrency=$CELERY_CONCURRENCY