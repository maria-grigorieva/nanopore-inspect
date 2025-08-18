#!/bin/bash

# Set default port if not provided
export FLASK_RUN_PORT="${FLASK_RUN_PORT:-5000}"

# Run the Flask app
python app.py