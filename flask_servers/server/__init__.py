from flask import Flask
from flask_cors import CORS
from logging.config import dictConfig

from server.config import SQLALCHEMY_DATABASE_URI, db
from server.models import Base
from server.views.litvar_views import litvar_views
from server.views.sample_views import sample_views
from server.views.vus_views import vus_views


def create_app():
    # TODO: write logging to file
    # logging configuration
    dictConfig({
        'version': 1,
        'formatters': {'default': {
            'format': '[%(asctime)s] %(levelname)s in %(module)s: %(message)s',
        }},
        'handlers': {'wsgi': {
            'class': 'logging.StreamHandler',
            'stream': 'ext://flask.logging.wsgi_errors_stream',
            'formatter': 'default'
        }},
        'root': {
            'level': 'INFO',
            'handlers': ['wsgi']
        }
    })

    app = Flask(__name__)
    app.config.from_pyfile('config.py')

    # Initialize SQLAlchemy with the Flask app
    db.init_app(app)

    CORS(app)

    app.register_blueprint(litvar_views, url_prefix='/litvar')
    app.register_blueprint(vus_views, url_prefix='/vus')
    app.register_blueprint(sample_views, url_prefix='/sample')

    return app
