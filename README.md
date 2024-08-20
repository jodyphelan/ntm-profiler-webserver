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
```