# MPapp

## Install env
```
conda env create -f MPapp_env.yml
```

## Install malaria-profiler
```
pip install git+https://github.com/jodyphelan/malaria-profiler.git
pip install git+https://github.com/jodyphelan/pathogen-profiler.git
```

## Run app
```
redis-server
sh celery.sh

cd MPapp
export FLASK_APP=MPapp
export FLASK_ENV=development
flask run --host=0.0.0.0


SCRIPT_NAME=/test-app gunicorn --worker-class gevent  --bind 127.0.0.1:5001  "mypackage:create_app()"
```