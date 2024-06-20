import os
from dotenv import load_dotenv

from flask_mail import Mail
from flask_sqlalchemy import SQLAlchemy

load_dotenv()

# Create an instance of SQLAlchemy
db = SQLAlchemy()

DEBUG = True
SECRET_KEY = os.getenv('SECRET_KEY')
CORS_HEADERS = 'Content-Type'
SQLALCHEMY_DATABASE_URI = os.getenv('SQLALCHEMY_DATABASE_URI')

MAIL_SERVER = "smtp.gmail.com"
MAIL_PORT = 587
MAIL_USERNAME = 'vus.curation.system@gmail.com'
MAIL_PASSWORD = os.getenv('MAIL_PASSWORD')
MAIL_DEFAULT_SENDER = 'vus.curation.system@gmail.com'
MAIL_USE_TLS = True
MAIL_USE_SSL = False

mail = Mail()
