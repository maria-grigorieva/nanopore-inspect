#!/bin/bash
set -e
#
# : "${REDIS_HOST:=redis}"
# : "${REDIS_PORT:=6379}"
#
# # Function to check Redis readiness
# wait_for_redis() {
#   echo "Waiting for Redis at $REDIS_HOST:$REDIS_PORT..."
#
#   while ! nc -z "$REDIS_HOST" "$REDIS_PORT"; do
#     echo "Redis is not available yet - sleeping"
#     sleep 1
#   done
#
#   echo "Redis is up!"
# }
#
# # Wait for Redis before starting services
# wait_for_redis

# Start services
# ./run_web.sh &
./run_worker.sh &
./run_beat.sh &

# Wait for all background processes
wait