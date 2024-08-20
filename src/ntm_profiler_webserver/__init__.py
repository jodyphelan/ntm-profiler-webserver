"""Package to create ntm-profiler webserver Flask app"""

__version__ = '0.1.0'

import os
from flask import Flask
from celery import Celery, Task

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

def create_app(test_config=None):
    # create and configure the app
    app = Flask(__name__, instance_relative_config=True)
    app.config.from_mapping(
        SECRET_KEY='dev',
        UPLOAD_FOLDER="/tmp/ntm-profiler",
        THREADS=4,
        RUN_SUBMISSION="local",
        APP_ROOT=os.path.dirname(os.path.abspath(__file__)),
        RESULTS_DIR=os.path.dirname(os.path.abspath(__file__))+"/static/results",
        REFERENCE_DIR=os.path.dirname(os.path.abspath(__file__))+"/static/reference_files",
    )

    app.config.from_mapping(
          CELERY=dict(
              broker_url= os.environ['REDIS_URL'] if 'REDIS_URL' in os.environ else 'redis://localhost:6379',
              result_backend=os.environ['REDIS_URL'] if 'REDIS_URL' in os.environ else 'redis://localhost:6379',
          ),
    )

    if test_config is None:
        # load the instance config, if it exists, when not testing
        app.config.from_pyfile('config.py', silent=True)
    else:
        # load the test config if passed in
        app.config.from_mapping(test_config)

    # ensure the instance folder exists
    try:
        os.makedirs(app.instance_path)
    except OSError:
        pass

    celery_init_app(app)

   
    from . import main
    app.register_blueprint(main.bp)

    return app
