# Create an instance of SQLAlchemy
from flask_sqlalchemy import SQLAlchemy

db = SQLAlchemy()

DEBUG = True
SECRET_KEY = 'h227GW-MI.5k}@H+Ppi"NOXO2#c)_z'
CORS_HEADERS = 'Content-Type'
SQLALCHEMY_DATABASE_URI = 'postgresql://postgres:21641@localhost:5432/vus-app-db'
