from flask import Flask


def create_app():
    app = Flask(__name__)
    app.config.from_pyfile('config.py')

    from .litvar_views import litvar_views
    from .entrez_views import entrez_views
    from .views import views

    app.register_blueprint(litvar_views, url_prefix='/litvar')
    app.register_blueprint(entrez_views, url_prefix='/entrez')
    app.register_blueprint(views, url_prefix='/')

    return app
