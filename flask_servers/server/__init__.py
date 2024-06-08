from datetime import datetime

from flask import Flask
from flask_cors import CORS
from logging.config import dictConfig

from flask_login import LoginManager

from server.config import SQLALCHEMY_DATABASE_URI, db
from server.models import Base, ScientificMembers
from server.services.clinvar_service import scheduled_clinvar_updates
from server.services.publications_service import check_for_new_litvar_publications
from server.views.auth_views import auth_views
from server.views.profile_views import profile_views
from server.views.publication_views import publication_views
from server.views.review_views import review_views
from server.views.sample_views import sample_views
from server.views.vus_views import vus_views

import atexit

from apscheduler.schedulers.background import BackgroundScheduler


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

    app.register_blueprint(publication_views, url_prefix='/publication')
    app.register_blueprint(vus_views, url_prefix='/vus')
    app.register_blueprint(sample_views, url_prefix='/sample')
    app.register_blueprint(auth_views, url_prefix='/auth')
    app.register_blueprint(profile_views)
    app.register_blueprint(review_views, url_prefix='/review')

    # specify user loader: tells Flask-Login how to find a specific user from the ID that is stored in their session cookie
    login_manager = LoginManager()
    login_manager.login_view = 'auth_views.login'
    login_manager.init_app(app)

    @login_manager.user_loader
    def load_scientific_member(scientific_member_id):
        # the user_id is just the primary key of our scientific members table
        return ScientificMembers.query.get(int(scientific_member_id))

    # schedule clinvar updates
    def scheduled_clinvar_updates_():
        with app.app_context():
            scheduled_clinvar_updates()

    # schedule litvar publications updates
    def scheduled_litvar_updates_():
        with app.app_context():
            check_for_new_litvar_publications()

    scheduler = BackgroundScheduler()
    scheduler.add_job(func=scheduled_clinvar_updates_, run_date=datetime.now())
    # 1 week
    scheduler.add_job(func=scheduled_clinvar_updates_, trigger="interval", seconds=604800)

    scheduler.add_job(func=scheduled_litvar_updates_, run_date=datetime.now())
    # 10 days
    scheduler.add_job(func=scheduled_litvar_updates_, trigger="interval", seconds=864000)
    scheduler.start()

    # Shut down the scheduler when exiting the app
    atexit.register(lambda: scheduler.shutdown())

    return app
