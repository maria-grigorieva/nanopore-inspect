#!/bin/bash

# Run the Celery beat scheduler
celery -A app.celery_app beat --loglevel=INFO