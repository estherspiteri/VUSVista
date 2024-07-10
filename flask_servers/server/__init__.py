from datetime import datetime

import pytz
from flask import Flask
from flask_cors import CORS
from logging.config import dictConfig

from flask_login import LoginManager

from server.config import *
from server.db_setup.populate_gene_annotations_table import store_gtf_file_in_db
from server.error_handlers import register_error_handlers
from server.models import Base, ScientificMembers, GeneAnnotations
from server.services.clinvar_service import scheduled_clinvar_updates
from server.services.publications_service import check_for_new_litvar_publications
from server.services.vus_preprocess_service import scheduled_file_upload_events
from server.views.auth_views import auth_views
from server.views.profile_views import profile_views
from server.views.publication_views import publication_views
from server.views.review_views import review_views
from server.views.sample_views import sample_views
from server.views.vus_views import vus_views

import atexit

from apscheduler.schedulers.background import BackgroundScheduler


def create_app():
    # logging configuration
    dictConfig({
        'version': 1,
        'formatters': {
            'default': {
                'format': '[%(asctime)s] %(levelname)s in %(module)s: %(message)s',
            },
            'file_formatter': {
                'format': '[%(asctime)s] %(levelname)s in %(module)s: %(message)s',
            },
        },
        'handlers': {
            'wsgi': {
                'class': 'logging.StreamHandler',
                'stream': 'ext://flask.logging.wsgi_errors_stream',
                'formatter': 'default',
            },
            'file': {
                'class': 'logging.handlers.RotatingFileHandler',
                'filename': 'logs/app.log',
                'maxBytes': 10240,
                'backupCount': 10,
                'formatter': 'file_formatter',
            },
        },
        'root': {
            'level': 'INFO',
            'handlers': ['wsgi', 'file'],
        },
    })

    app = Flask(__name__)
    app.config.from_pyfile('config.py')

    # Initialize Flask-Mail
    mail.init_app(app)

    # Initialize SQLAlchemy with the Flask app
    db.init_app(app)

    with app.app_context():
        gene_annotations = db.session.query(GeneAnnotations).all()
        if len(gene_annotations) == 0:
            # Populate gene annotations table
            store_gtf_file_in_db()

    CORS(app)

    app.register_blueprint(publication_views, url_prefix='/publication')
    app.register_blueprint(vus_views, url_prefix='/vus')
    app.register_blueprint(sample_views, url_prefix='/sample')
    app.register_blueprint(auth_views, url_prefix='/auth')
    app.register_blueprint(profile_views, url_prefix='/user')
    app.register_blueprint(review_views, url_prefix='/review')

    # Register error handlers
    register_error_handlers(app)

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

    # schedule file upload event check
    def scheduled_file_upload_events_():
        with app.app_context():
            scheduled_file_upload_events()

    with app.app_context():
        scheduler = BackgroundScheduler()

        # Get the current time with timezone information
        timezone = pytz.timezone('Europe/Amsterdam')
        run_date = datetime.now(timezone)

        scheduler.add_job(func=scheduled_clinvar_updates_, run_date=run_date)
        # 1 week
        scheduler.add_job(func=scheduled_clinvar_updates_, trigger="interval", seconds=604800)

        scheduler.add_job(func=scheduled_litvar_updates_, run_date=run_date)
        # 10 days
        scheduler.add_job(func=scheduled_litvar_updates_, trigger="interval", seconds=864000)

        scheduler.add_job(func=scheduled_file_upload_events_, run_date=run_date)
        # 1 min
        scheduler.add_job(func=scheduled_file_upload_events_, trigger="interval", seconds=180, max_instances=2)

        scheduler.start()

        # Shut down the scheduler when exiting the app
        atexit.register(lambda: scheduler.shutdown())

    return app
