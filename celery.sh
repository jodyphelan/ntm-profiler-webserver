cd bin
celery -A celery_queue worker --loglevel=INFO --concurrency=1 --time-limit 3600

